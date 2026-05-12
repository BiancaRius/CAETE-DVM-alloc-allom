module carbon_allocation_offline_kernel

  !==========================================================================
  ! Carbon allocation offline solver - Step 1
  !==========================================================================
  ! Purpose
  ! -------
  ! This module implements a self-contained woody-plant carbon allocation
  ! kernel based on the LPJ-style allometric allocation logic.
  !
  ! This file is intentionally verbose. The comments are part of the model
  ! documentation and are meant to help a future reader understand why each
  ! equation appears in the code.
  !
  ! Scope of this first version
  ! ---------------------------
  ! This kernel solves allocation for ONE average woody individual over ONE
  ! allocation period. The allocation period can be annual, monthly, seasonal,
  ! or any other interval. The equations do not know the calendar frequency.
  ! The caller only needs to provide the carbon available over that period.
  !
  ! Key idea
  ! --------
  ! The model has an available carbon increment, here called c_available:
  !
  !     c_available = dL + dR + dS
  !
  ! where:
  !
  !     dL = increment to leaf carbon mass
  !     dR = increment to fine-root carbon mass
  !     dS = increment to sapwood carbon mass
  !
  ! For woody plants, the allocation cannot be arbitrary. The final plant
  ! state after allocation must satisfy allometric constraints linking leaf
  ! mass, fine-root mass, sapwood mass, heartwood mass, height, stem diameter,
  ! and sapwood cross-sectional area.
  !
  ! The normal-allocation problem is reduced to one unknown:
  !
  !     x = dL
  !
  ! Once x is known, dR follows from the leaf-to-root allometry, and dS follows
  ! from carbon conservation. The bisection method is then used to find the x
  ! that makes the final structure consistent with both stem geometry and the
  ! pipe-model constraint.
  !
  ! Important modeling decision
  ! ---------------------------
  ! This module should not be called automatically every day unless the model
  ! has a storage/buffer pool. With a daily carbon increment, c_available may be
  ! too small for structural growth, causing frequent abnormal allocation.
  ! A safer first implementation is to accumulate daily carbon fluxes outside
  ! this kernel and call this kernel monthly, seasonally, or annually.
  !==========================================================================

   use, intrinsic :: iso_fortran_env, only: real64

   implicit none

   private   

   public :: Parameters
   public :: PlantCarbonState
   public :: AllocationOutput
   public :: allocate
   public :: allocation_residual
   public :: leaf_requirement


  !--------------------------------------------------------------------------
  ! Numerical constants
  !--------------------------------------------------------------------------

   real(real64), parameter :: pi = 3.1415926535897932384626433832795_real64

   ! Default bisection settings.
   ! x_tolerance controls convergence in terms of the unknown dL.
   ! f_tolerance controls convergence in terms of the residual f(dL).

   ! FIRST TRIAL VALUES (MIGHT CHANGE IN THE FUTURE)
   real(real64), parameter :: default_x_tolerance = 1.0e-8_real64 !! ATTENTION: must be in the units you use for carbon stocks
   real(real64), parameter :: default_f_tolerance = 1.0e-10_real64

   integer, parameter :: default_max_iterations = 200
   integer, parameter :: default_scan_segments  = 100

  !--------------------------------------------------------------------------
  ! Input parameters
  !--------------------------------------------------------------------------

   type :: Parameters

      !!! PLS-specific parameters !!!
      !-----------------------------!

      ! Specific leaf area (PLS specific)
      ! Units must be consistent with leaf_mass and area.
      ! Example: if leaf_mass is gC individual-1 and SLA is m2 gC-1,
      ! then leaf_area = leaf_mass * sla is m2 individual-1.
      real(real64) :: sla

      ! Wood density (PLS specific)
      ! Units must be consistent with carbon mass and volume.
      ! Example: gC m-3 if carbon pools are in gC.
      real(real64) :: wood_density

      ! Leaf area to sapwood cross-sectional area ratio (global value)
      ! This is the pipe-model coefficient.
      ! Pipe model:
      !     leaf_area = latosa * sapwood_cross_sectional_area
      real(real64) :: latosa


     !!! Global parameters !!!
     !-----------------------!
     
     ! Leaf-to-fine-root mass ratio for the allocation period.
     ! Functional balance:
     !     leaf_mass = leaf_to_root_ratio * root_mass
     ! In LPJ-style code this is often computed as:
     !     leaf_to_root_ratio = max(lm2rm_max * water_scalar, 0.1)
     real(real64) :: leaf_to_root_ratio

     ! Height-diameter allometry coefficient.
     ! Stem mechanics:
     !     height = allom2 * stem_diameter**allom3
     real(real64) :: allom2

     ! Height-diameter allometry exponent.
     ! Stem mechanics:
     !     height = allom2 * stem_diameter**allom3
     real(real64) :: allom3

   end type Parameters



  !--------------------------------------------------------------------------
  ! Carbon state of one average woody individual
  !
  ! These variables represent the current carbon pools at the moment when the
  ! allocation routine is called.  
  !--------------------------------------------------------------------------

   type :: PlantCarbonState

     ! Leaf carbon mass 
     real(real64) :: leaf_mass

     ! Fine-root carbon mass 
     real(real64) :: root_mass

     ! Sapwood carbon mass 
     real(real64) :: sapwood_mass

     ! Heartwood carbon mass 
     real(real64) :: heartwood_mass

     ! Plant height.
     ! The variable comes from the current state, before allocation
     ! Required for detecting the minimum leaf mass needed to maintain the
     ! currently existing sapwood mass before normal allocation.
     real(real64) :: height

   end type PlantCarbonState



  !--------------------------------------------------------------------------
  ! Output from the allocation routine
  !--------------------------------------------------------------------------

   type :: AllocationOutput

      !! Allocation pathway used by the solver.
      !
      ! .true.  = normal allocation was used.
      !           In this pathway, the plant had enough available carbon to
      !           attempt a fully allometric growth solution with positive
      !           increments to the living tissues:
      !
      !               delta_leaf    > 0
      !               delta_root    > 0
      !               delta_sapwood > 0
      !
      !           The unknown leaf increment, delta_leaf, was solved by the
      !           bisection method so that the final plant state satisfies the
      !           leaf-root functional balance, the pipe model, and the stem
      !           height-diameter/volume constraints simultaneously.
      !
      ! .false. = abnormal allocation was used.
      !           In this pathway, the normal allometric problem was not feasible
      !           with the available carbon. This can happen when the carbon
      !           available over the allocation period is too small to maintain
      !           or increase leaf, root, and sapwood pools while satisfying the
      !           allometric constraints. The solver then applies a corrective
      !           allocation that may reduce one or more pools to restore the
      !           leaf-root and pipe-model relationships.

      ! True if the normal allometric allocation problem was solved by bisection.
      ! False if abnormal allocation was used.
      logical :: normal_allocation = .false.

      !! Convergence status of the numerical solver.
      ! For normal allocation:
      !     .true.  = the bisection method found a value of delta_leaf that
      !               satisfies the nonlinear allometric residual equation within
      !               the prescribed numerical tolerances.
      !
      !     .false. = the bisection method did not converge within the maximum
      !               number of iterations, or no valid sign-changing bracket was
      !               found before bisection. In that case, the diagnostic message
      !               should be inspected before trusting the allocation result.
      !
      ! For abnormal allocation:
      !     This flag is set to .true. after the fallback algebraic correction is
      !     completed, because no bisection root-finding is required. In this case,
      !     normal_allocation should be used together with converged to interpret
      !     what happened.
      logical :: converged = .false.

      ! Number of bisection iterations used in the normal-allocation case.
      ! Used for numerical diagnose. If iterations is too high it my indicate some issue
      integer :: iterations = 0

      ! Carbon increments of the average individual over the allocation period.
      !! _real64 is used for precision safety 
      real(real64) :: delta_leaf    = 0.0_real64 
      real(real64) :: delta_root    = 0.0_real64
      real(real64) :: delta_sapwood = 0.0_real64

      ! New carbon pools after allocation.
      real(real64) :: leaf_mass_new      = 0.0_real64
      real(real64) :: root_mass_new      = 0.0_real64
      real(real64) :: sapwood_mass_new   = 0.0_real64
      real(real64) :: heartwood_mass_new = 0.0_real64

      ! New structural variables inferred from the final allometric state.
      real(real64) :: sapwood_area_new = 0.0_real64
      real(real64) :: height_new       = 0.0_real64
      real(real64) :: stem_diameter_new = 0.0_real64

      !! Diagnostic residuals and balance checks.
      ! These variables are numerical diagnostics used to verify
      ! whether the allocation result satisfies the equations that it is supposed
      ! to satisfy.
      !
      ! In an ideal mathematical solution, the residuals would be exactly zero.
      ! In floating-point computation, they are expected to be very close to zero,
      ! but not necessarily exactly zero, because the bisection method stops when
      ! numerical tolerances are reached.
      !
      ! carbon_balance_error checks conservation of the available carbon over the
      ! living allocation increments:
      !     carbon_balance_error = c_available
      !                            - (delta_leaf + delta_root + delta_sapwood)
      ! For normal allocation, this should be approximately zero.
      !
      ! leaf_root_residual checks whether the final leaf and fine-root pools obey
      ! the functional balance relationship:
      !     leaf_mass_new = leaf_to_root_ratio * root_mass_new
      !
      ! Therefore:
      !     leaf_root_residual = leaf_mass_new
      !                          - leaf_to_root_ratio * root_mass_new
      ! For normal allocation, this should be approximately zero.
      !
      ! pipe_model_residual checks whether the final leaf area and sapwood area
      ! obey the pipe-model relationship:
      !     leaf_area_new = latosa * sapwood_area_new
      !
      ! Therefore:
      !     pipe_model_residual = leaf_area_new
      !                           - latosa * sapwood_area_new
      ! For normal allocation, this should be approximately zero.
      !
      ! final_root_residual is reserved for additional root-related diagnostics
      ! during integration with the full model. It can be removed if it remains
      ! redundant with leaf_root_residual.

      real(real64) :: carbon_balance_error = 0.0_real64
      real(real64) :: leaf_root_residual   = 0.0_real64
      real(real64) :: pipe_model_residual  = 0.0_real64
      real(real64) :: final_root_residual  = 0.0_real64

      ! Residual of the nonlinear equation f(dL) = 0 at the selected solution.
      real(real64) :: allocation_residual_final = 0.0_real64

      !! Bounds used by the normal-allocation solver 
      ! Inside this interval, the solver searches for the
      ! delta_leaf value that also makes the stem geometry and pipe-model
      ! constraints mutually consistent (more detailed explanations below)
      real(real64) :: lower_bound_delta_leaf = 0.0_real64
      real(real64) :: upper_bound_delta_leaf = 0.0_real64

      ! Diagnostic message
      character(len=160) :: message = ""

   end type AllocationOutput


   contains

   !==========================================================================
   !> Compute the leaf mass required to maintain the existing sapwood mass.
   !==========================================================================
      function leaf_requirement(state, params) result(leaf_required)

         type(PlantCarbonState), intent(in) :: state
         type(Parameters), intent(in) :: params

         real(real64) :: leaf_required

         !-----------------------------------------------------------------------
         ! Derivation
         ! ----------
         ! The pipe model states:
         !     leaf_area = latosa * sapwood_area
         !
         ! Leaf area is:
         !     leaf_area = leaf_mass * SLA
         !
         ! Sapwood volume is:
         !     sapwood_volume = height * sapwood_area
         !
         ! Sapwood mass is:
         !     sapwood_mass = wood_density * sapwood_volume
         !                  = wood_density * height * sapwood_area
         !
         ! Therefore:
         !     sapwood_area = sapwood_mass / (wood_density * height)
         !
         ! Substitute this into the pipe model:
         !     leaf_mass * SLA = latosa * sapwood_mass / (wood_density * height)
         !
         ! Solve for leaf_mass:
         !     leaf_mass_required = latosa * sapwood_mass /
         !                          (wood_density * height * SLA)
         !
         ! This is the minimum leaf mass needed to keep the current sapwood mass
         ! consistent with the pipe model, assuming no new sapwood is produced.
         !-----------------------------------------------------------------------

         if (state%height <= 0.0_real64) then
            leaf_required = 0.0_real64
         else
            leaf_required = params%latosa * state%sapwood_mass / &
                           (params%wood_density * state%height * params%sla)
         end if

      end function leaf_requirement


   !==========================================================================
   !> Nonlinear residual f(dL) used by the normal-allocation bisection solver.
   !==========================================================================
      function allocation_residual(state, params, c_available, delta_leaf) result(f)

         type(PlantCarbonState),    intent(in) :: state
         type(Parameters), intent(in) :: params
         real(real64),              intent(in) :: c_available
         real(real64),              intent(in) :: delta_leaf

         real(real64) :: f

         real(real64) :: leaf_new
         real(real64) :: root_new
         real(real64) :: delta_root
         real(real64) :: delta_sapwood
         real(real64) :: sapwood_new
         real(real64) :: stem_carbon_total_new
         real(real64) :: a1
         real(real64) :: a2
         real(real64) :: a3
         real(real64) :: pi_over_four
         real(real64) :: height_power_pipe_model
         real(real64) :: height_power_total_stem

         !-----------------------------------------------------------------------
         ! Unknown
         ! -------
         ! The bisection method solves for:
         !
         !     x = delta_leaf = dL
         !
         ! Once x is chosen, all other increments are determined.
         !-----------------------------------------------------------------------

         leaf_new = state%leaf_mass + delta_leaf

         !-----------------------------------------------------------------------
         ! Fine-root increment from functional balance
         ! -------------------------------------------
         ! Functional balance after allocation:
         !     leaf_new = leaf_to_root_ratio * root_new
         !
         ! Therefore:
         !     root_new = leaf_new / leaf_to_root_ratio
         !
         ! The increment is:
         !     delta_root = root_new - root_old
         !
         ! which gives:
         !     delta_root = (leaf_old + delta_leaf) / leaf_to_root_ratio
         !                  - root_old
         !-----------------------------------------------------------------------
         root_new   = leaf_new / params%leaf_to_root_ratio
         delta_root = root_new - state%root_mass

         !-----------------------------------------------------------------------
         ! Sapwood increment from carbon conservation
         ! ------------------------------------------
         ! Carbon conservation over the allocation period:
         !     c_available = delta_leaf + delta_root + delta_sapwood
         !
         ! Therefore:
         !     delta_sapwood = c_available - delta_leaf - delta_root
         !
         ! Substitute delta_root:
         !     delta_sapwood = c_available - delta_leaf
         !                     - ((leaf_old + delta_leaf) / leaf_to_root_ratio
         !                        - root_old)
         !
         ! Final sapwood mass:
         !     sapwood_new = sapwood_old + delta_sapwood
         !
         ! which expands to:
         !     sapwood_new = sapwood_old + c_available - delta_leaf
         !                   - ((leaf_old + delta_leaf) / leaf_to_root_ratio)
         !                   + root_old
         !
         ! It is the final sapwood mass written as a function of delta_leaf.
         !-----------------------------------------------------------------------

         delta_sapwood = c_available - delta_leaf - delta_root
         sapwood_new   = state%sapwood_mass + delta_sapwood

         ! Guard against non-physical values during root search.
         ! The bisection bounds should normally avoid these cases, but this protects
         ! the residual function if it is called outside the valid interval.
         if (leaf_new <= 0.0_real64 .or. sapwood_new <= 0.0_real64) then
            f = huge(1.0_real64)
            return
         end if

         ! ------------------------------------------------
         ! The normal-allocation solution requires consistency between:
         !   (1) total stem geometry, using sapwood + heartwood, and
         !   (2) the pipe model, using leaf mass + sapwood mass.
         !
         ! Both routes can be rearranged to compute the same quantity (plant height raised to the exponent (1 + 2 / allom3)):
         !     height_new**(1 + 2 / allom3)
         !
         ! The residual is:
         !     f(delta_leaf) = height_power_total_stem
         !                     - height_power_pipe_model
         !
         ! The correct allocation satisfies:
         !     f(delta_leaf) = 0
         !-----------------------------------------------------------------------

         a1 = 2.0_real64 / params%allom3
         a2 = 1.0_real64 + a1
         a3 = params%allom2**a1
         pi_over_four = pi / 4.0_real64

         !-----------------------------------------------------------------------
         ! Expression 1: total stem geometry
         ! ---------------------------------
         ! Height-diameter allometry:
         !     height = allom2 * diameter**allom3
         !
         ! Stem volume:
         !     stem_volume = height * pi * diameter**2 / 4
         !
         ! Total stem carbon includes living sapwood and non-living heartwood:
         !     stem_carbon_total_new = sapwood_new + heartwood_old
         !
         ! Wood density:
         !     wood_density = stem_carbon_total_new / stem_volume
         !
         ! Combining these equations and eliminating diameter gives:
         !     height_new**(1 + 2 / allom3)
         !       = allom2**(2 / allom3)
         !         * ((sapwood_new + heartwood_old) / wood_density)
         !         / (pi / 4)
         !
         ! This is height_power_total_stem.
         !-----------------------------------------------------------------------

         stem_carbon_total_new = sapwood_new + state%heartwood_mass

         ! Note this is not the actual height, but height raised to the power (1 + 2 / allom3), which is the form that allows direct comparison with the pipe-model expression.
         height_power_total_stem = a3 * (stem_carbon_total_new / params%wood_density) &
                                    / pi_over_four

         !-----------------------------------------------------------------------
         ! Expression 2: pipe model plus sapwood volume
         ! --------------------------------------------
         ! Pipe model:
         !     leaf_area = latosa * sapwood_area
         !
         ! Leaf area:
         !     leaf_area = leaf_new * SLA
         !
         ! Therefore:
         !     sapwood_area = leaf_new * SLA / latosa
         !
         ! Sapwood volume:
         !     sapwood_volume = height * sapwood_area
         !
         ! Sapwood mass:
         !     sapwood_new = wood_density * sapwood_volume
         !                 = wood_density * height * sapwood_area
         !
         ! Solve for height:
         !     height = sapwood_new / (wood_density * sapwood_area)
         !
         ! Substitute sapwood_area:
         !     height = sapwood_new /
         !              (wood_density * leaf_new * SLA / latosa)
         !
         ! Raise both sides to:
         !     1 + 2 / allom3
         ! to match the total-stem expression:
         !     height_new**(1 + 2 / allom3)
         !       = [ sapwood_new /
         !           (wood_density * leaf_new * SLA / latosa) ]
         !         **(1 + 2 / allom3)
         !
         ! This is height_power_pipe_model.
         !-----------------------------------------------------------------------

         ! Note this is not actual height, but height raised to the power (1 + 2 / allom3), which is the form that allows direct comparison with the total-stem expression.
         height_power_pipe_model = &
            (sapwood_new / (leaf_new * params%sla * params%wood_density / params%latosa))**a2

         ! The allocation residual is zero only when both structural routes agree.
         f = height_power_total_stem - height_power_pipe_model

      end function allocation_residual


   !==========================================================================
   !> Main allocation routine for one woody individual and one allocation period.
   !==========================================================================
      subroutine allocate(state, params, c_available, result)

         type(PlantCarbonState),    intent(in)  :: state
         type(Parameters), intent(in)  :: params
         real(real64),              intent(in)  :: c_available
         type(AllocationOutput),    intent(out) :: result

         real(real64) :: leaf_required
         real(real64) :: delta_leaf_min
         real(real64) :: delta_root_min
         real(real64) :: lower
         real(real64) :: upper
         real(real64) :: f_lower
         real(real64) :: f_upper
         real(real64) :: left
         real(real64) :: right
         real(real64) :: mid
         real(real64) :: f_left
         real(real64) :: f_mid
         real(real64) :: scan_step
         real(real64) :: previous_x
         real(real64) :: previous_f
         real(real64) :: scan_interval_width
         real(real64) :: n_scan_segments_real
         logical      :: bracket_found
         integer      :: i

         ! Initialize output.
         result = AllocationOutput()

         ! Basic checks.
         if (params%sla <= 0.0_real64 .or. params%latosa <= 0.0_real64 .or. &
               params%wood_density <= 0.0_real64 .or. params%leaf_to_root_ratio <= 0.0_real64 .or. &
               params%allom2 <= 0.0_real64 .or. params%allom3 <= 0.0_real64) then
               result%message = "Invalid parameter value. All core allometric parameters must be positive."
               return
         end if

         if (state%height <= 0.0_real64) then
               result%message = "Invalid initial state. Height must be positive for woody    allocation."
               return
         end if

         !-----------------------------------------------------------------------
         ! Step 1: compute minimum leaf and root increments needed for normal
         ! allocation.
         !-----------------------------------------------------------------------
         leaf_required  = leaf_requirement(state, params)
         delta_leaf_min = leaf_required - state%leaf_mass

            ! Root mass needed to support leaf_required under functional balance:
            !     root_required = leaf_required / leaf_to_root_ratio
            !
            ! Minimum root increment:
            !     delta_root_min = root_required - root_old
         delta_root_min = leaf_required / params%leaf_to_root_ratio - state%root_mass
            
         !-----------------------------------------------------------------------
            ! Step 2: decide whether the normal-allocation problem is feasible.
            !-----------------------------------------------------------------------
            ! Normal allocation requires a feasible interval for delta_leaf.
            ! The lower bound is defined by three constraints:
            ! 1. Leaf increment should not be negative:
            !    NOTE: it can be 0 if the current leaf mass is already sufficient to maintain the existing sapwood under the pipe model, but it cannot be negative 
            !        delta_leaf >= 0
            !
            ! 2. Leaf mass after allocation must be sufficient to maintain the
            !    already existing sapwood under the pipe model:
            !        delta_leaf >= delta_leaf_min
            !
            ! 3. Root increment should not be negative:
            !     NOTE: it can be 0 if the current root mass is already sufficient to maintain the existing leaf mass under functional balance, but it cannot be negative
            !        delta_root >= 0
            !
            !    Since:
            !        root_new = leaf_new / leaf_to_root_ratio
            !        delta_root = root_new - root_old
            !
            !    then:
            !        delta_root >= 0
            !        (leaf_old + delta_leaf) / leaf_to_root_ratio - root_old >= 0
            !
            !    therefore:
            !        delta_leaf >= root_old * leaf_to_root_ratio - leaf_old
            !
            ! The upper bound is the maximum delta_leaf that still leaves
            ! non-negative carbon for new sapwood.
            !-----------------------------------------------------------------------

         lower = max( &
         0.0_real64, &
         delta_leaf_min, &
         state%root_mass * params%leaf_to_root_ratio - state%leaf_mass)

            !--------------------------------------------------------------------
            ! Upper bound for dL
            ! ------------------
            ! At the upper bound, no carbon is left for new sapwood:
            !     c_available = dL + dR
            ! with:
            !     dR = (leaf_old + dL) / leaf_to_root_ratio - root_old
            !
            ! Substitute dR:
            !     c_available = dL + (leaf_old + dL) / leaf_to_root_ratio - root_old
            !
            ! Rearrange:
            !     c_available + root_old - leaf_old / leaf_to_root_ratio
            !       = dL * (1 + 1 / leaf_to_root_ratio)
            !
            ! Therefore:
            !     dL_max = (c_available - leaf_old / leaf_to_root_ratio + root_old)
            !              / (1 + 1 / leaf_to_root_ratio)
            !
            ! This is the maximum possible dL before delta_sapwood becomes zero.
            !--------------------------------------------------------------------

            upper = (c_available - state%leaf_mass / params%leaf_to_root_ratio + &
               state%root_mass) / &
               (1.0_real64 + 1.0_real64 / params%leaf_to_root_ratio)

            result%lower_bound_delta_leaf = lower
            result%upper_bound_delta_leaf = upper

         if (upper > lower) then

            result%normal_allocation = .true.



            ! if (upper <= lower) then
            !    result%message = "Normal allocation requested, but the bisection interval is invalid because the upper bound is *lower/equal to* lower bound."
            !    call abnormal_allocation(state, params, c_available, result)
            !    return
            ! end if


            !--------------------------------------------------------------------
            ! Step 3: find a sign-changing bracket for f(delta_leaf).
            !--------------------------------------------------------------------
            ! The bisection method needs two values of delta_leaf, called left and
            ! right, such that f(left) and f(right) have opposite signs. This sign
            ! change indicates that f(delta_leaf) crosses zero between them.
            !
            ! The full interval [lower, upper] may not show a sign change at its
            ! endpoints, even if a root exists somewhere inside it. Therefore, we
            ! scan the interval in smaller segments and look for a subinterval where
            ! the sign changes. That subinterval is then used as the initial bracket
            ! for bisection.
            !--------------------------------------------------------------------
            f_lower = allocation_residual(state, params, c_available, lower)
            f_upper = allocation_residual(state, params, c_available, upper)

            bracket_found = .false.

            !! Start the interval scan at the lower bound.
            ! previous_x stores the last tested delta_leaf value, and previous_f
            ! stores the residual evaluated at that value. These are compared with
            ! the next scanned point to detect a sign change.
            previous_x = lower
            previous_f = f_lower

            ! Compute the total width of the delta_leaf interval to be scanned.
            scan_interval_width = upper - lower
            ! Convert the integer number of scan segments to real64 before division.
            n_scan_segments_real = real(default_scan_segments, real64)
            ! Each scan step is one equal subdivision of the full interval.
            scan_step = scan_interval_width / n_scan_segments_real

            if (abs(f_lower) <= 0.0_real64) then
               left = lower
               right = lower
               bracket_found = .true.
            else if (f_lower * f_upper <= 0.0_real64) then
               left = lower
               right = upper
               bracket_found = .true.
            else
               do i = 1, default_scan_segments
                  mid   = lower + scan_step * real(i, real64)
                  f_mid = allocation_residual(state, params, c_available, mid)

                  if (previous_f * f_mid <= 0.0_real64) then
                     left = previous_x
                     right = mid
                     bracket_found = .true.
                     exit
                  end if

                  previous_x = mid
                  previous_f = f_mid
               end do
            end if

            if (.not. bracket_found) then
               result%message = "No sign-changing bracket found for normal allocation; using abnormal allocation."
               call abnormal_allocation(state, params, c_available, result)
               return
            end if

            !--------------------------------------------------------------------
            ! Step 4: solve f(dL) = 0 by bisection.
            !--------------------------------------------------------------------
            ! Bisection repeatedly cuts the current bracket in half:
            !
            !     mid = 0.5 * (left + right)
            !
            ! Then it keeps the half-interval where the sign change remains.
            ! This is robust because it does not require derivatives.
            !--------------------------------------------------------------------

            if (abs(left - right) <= default_x_tolerance) then
               result%delta_leaf = left
               result%converged = .true.
               result%iterations = 0
            else
               f_left = allocation_residual(state, params, c_available, left)

                do i = 1, default_max_iterations
                    mid   = 0.5_real64 * (left + right)
                      f_mid = allocation_residual(state, params, c_available, mid)

                    if (abs(f_mid) <= default_f_tolerance .or. &
                        abs(right - left) <= default_x_tolerance) then
                        result%delta_leaf = mid
                        result%converged = .true.
                        result%iterations = i
                        exit
                    end if

                      if (f_left * f_mid <= 0.0_real64) then
                        right = mid
                      else
                         left = mid
                         f_left = f_mid
                      end if
                   end do

                if (.not. result%converged) then
                    result%delta_leaf = 0.5_real64 * (left + right)
                    result%iterations = default_max_iterations
                    result%message = "Bisection reached maximum iterations; returning best midpoint estimate."
                end if
            end if

            ! Use the solved dL to compute the other increments and final state.
            call final_allocation(state, params, c_available, result%delta_leaf, result)

            if (result%message == "") then
               result%message = "Normal allocation solved."
            end if

         else

            !--------------------------------------------------------------------
            ! Abnormal allocation
            ! -------------------
            ! Normal allocation is not feasible. The plant does not have enough
            ! carbon to increase leaf, root, and sapwood while also maintaining all
            ! allometric constraints. The model then reallocates/reduces some pools
            ! to restore allometry.
            !--------------------------------------------------------------------
            call abnormal_allocation(state, params, c_available, result)

         end if

      end subroutine allocate

   !==========================================================================
   !> Compute final increments and diagnostic residuals from a selected delta Leaf
   !==========================================================================
      subroutine final_allocation(state, params, c_available, delta_leaf, result)

         type(PlantCarbonState),    intent(in)    :: state
         type(Parameters), intent(in)    :: params
         real(real64),              intent(in)    :: c_available
         real(real64),              intent(in)    :: delta_leaf
         type(AllocationOutput),    intent(inout) :: result

         real(real64) :: leaf_area_new
         real(real64) :: sapwood_area_from_mass

         result%delta_leaf = delta_leaf

         ! Fine-root increment from functional balance.
         result%delta_root = (state%leaf_mass + result%delta_leaf) / &
                           params%leaf_to_root_ratio - state%root_mass

         ! Sapwood increment from carbon conservation.
         result%delta_sapwood = c_available - result%delta_leaf - result%delta_root

         ! Final pools.
         result%leaf_mass_new    = state%leaf_mass    + result%delta_leaf
         result%root_mass_new    = state%root_mass    + result%delta_root
         result%sapwood_mass_new = state%sapwood_mass + result%delta_sapwood

         ! In normal allocation, heartwood does not receive new carbon directly.
         ! It changes elsewhere through sapwood-to-heartwood conversion or turnover.
         result%heartwood_mass_new = state%heartwood_mass

         ! Pipe-model sapwood cross-sectional area.
         result%sapwood_area_new = result%leaf_mass_new * params%sla / params%latosa

         ! Height from sapwood mass and sapwood cross-sectional area:
         !     sapwood_mass = wood_density * height * sapwood_area
         !
         ! Therefore:
         !     height = sapwood_mass / (wood_density * sapwood_area)
         if (result%sapwood_area_new > 0.0_real64) then
            result%height_new = result%sapwood_mass_new / &
                                 (params%wood_density * result%sapwood_area_new)
         else
            result%height_new = 0.0_real64
         end if

         ! Stem diameter from height-diameter allometry:
         !     height = allom2 * diameter**allom3
         !
         ! Therefore:
         !     diameter = (height / allom2)**(1 / allom3)
         if (result%height_new > 0.0_real64) then
            result%stem_diameter_new = (result%height_new / params%allom2)**(1.0_real64 / params%allom3)
         else
            result%stem_diameter_new = 0.0_real64
         end if

         ! Diagnostics.
         result%carbon_balance_error = c_available - &
            (result%delta_leaf + result%delta_root + result%delta_sapwood)

         result%leaf_root_residual = result%leaf_mass_new - &
            params%leaf_to_root_ratio * result%root_mass_new

         leaf_area_new = result%leaf_mass_new * params%sla
         sapwood_area_from_mass = 0.0_real64
         if (result%height_new > 0.0_real64) then
            sapwood_area_from_mass = result%sapwood_mass_new / &
                                    (params%wood_density * result%height_new)
         end if

         result%pipe_model_residual = leaf_area_new - params%latosa * sapwood_area_from_mass

         result%allocation_residual_final = allocation_residual(state, params, c_available, result%delta_leaf)

      end subroutine final_allocation

   !==========================================================================
   !> Abnormal woody allocation.
   !==========================================================================
      subroutine abnormal_allocation(state, params, c_available, result)

         type(PlantCarbonState),    intent(in)    :: state
         type(Parameters), intent(in)    :: params
         real(real64),              intent(in)    :: c_available
         type(AllocationOutput),    intent(inout) :: result

         real(real64) :: delta_leaf
         real(real64) :: delta_root
         real(real64) :: sapwood_required

         result%normal_allocation = .false.
         result%converged = .true.

         !-----------------------------------------------------------------------
         ! Abnormal allocation derivation
         ! ------------------------------
         ! Normal allocation cannot be solved with positive increments to leaf,
         ! fine root, and sapwood. The fallback logic first tries to use the
         ! available carbon to restore the leaf-root functional balance only:
         !     c_available = dL + dR
         ! and:
         !     leaf_old + dL = leaf_to_root_ratio * (root_old + dR)
         !
         ! From the second equation:
         !     dR = (leaf_old + dL) / leaf_to_root_ratio - root_old
         !
         ! Substitute into c_available = dL + dR:
         !     c_available = dL + (leaf_old + dL) / leaf_to_root_ratio - root_old
         !
         ! Rearrange:
         !     c_available + root_old - leaf_old / leaf_to_root_ratio
         !       = dL * (1 + 1 / leaf_to_root_ratio)
         !
         ! Therefore:
         !     dL = (c_available - leaf_old / leaf_to_root_ratio + root_old)
         !          / (1 + 1 / leaf_to_root_ratio)
         !-----------------------------------------------------------------------

         delta_leaf = (c_available - state%leaf_mass / params%leaf_to_root_ratio + &
                        state%root_mass) / (1.0_real64 + 1.0_real64 / params%leaf_to_root_ratio)

         if (delta_leaf > 0.0_real64) then
            ! There is still positive allocation to leaves. The remaining carbon is
            ! assigned to fine roots.
            delta_root = c_available - delta_leaf

            if (delta_root < 0.0_real64) then
               ! If restoring the leaf-root ratio would require negative root
               ! allocation, allocate all available carbon to leaves and reduce roots
               ! to the value implied by the final leaf mass.
               delta_leaf = c_available
               delta_root = (state%leaf_mass + delta_leaf) / params%leaf_to_root_ratio - &
                           state%root_mass
            end if

         else

            ! Leaf increment is negative. Allocate available carbon to roots, then
            ! reduce leaves to restore the leaf-root ratio.
            delta_root = c_available
            delta_leaf = params%leaf_to_root_ratio * (state%root_mass + delta_root) - &
                        state%leaf_mass

         end if

         result%delta_leaf = delta_leaf
         result%delta_root = delta_root

         result%leaf_mass_new = state%leaf_mass + result%delta_leaf
         result%root_mass_new = state%root_mass + result%delta_root

         !-----------------------------------------------------------------------
         ! Sapwood adjustment under abnormal allocation
         ! --------------------------------------------
         ! After leaf and root pools are adjusted, compute the sapwood mass required
         ! by the pipe model at the CURRENT height:
         !
         !     sapwood_area_required = leaf_new * SLA / latosa
         !
         !     sapwood_required = wood_density * height_old * sapwood_area_required
         !
         ! Therefore:
         !
         !     sapwood_required = leaf_new * SLA / latosa * wood_density * height_old
         !
         ! The sapwood increment is then:
         !
         !     dS = sapwood_required - sapwood_old
         !
         ! In abnormal allocation this value is expected to be negative, meaning
         ! excess sapwood is converted into heartwood.
         !-----------------------------------------------------------------------

         sapwood_required = result%leaf_mass_new * params%sla / params%latosa * &
                        params%wood_density * state%height

         result%delta_sapwood = sapwood_required - state%sapwood_mass

         result%sapwood_mass_new = state%sapwood_mass + result%delta_sapwood

         ! Convert reduced sapwood to heartwood if delta_sapwood is negative.
         ! If numerical conditions produce positive delta_sapwood here, we do not
         ! add it to heartwood; the positive sapwood increment is already included
         ! in sapwood_mass_new.
         result%heartwood_mass_new = state%heartwood_mass + max(-result%delta_sapwood, 0.0_real64)

         ! Structural variables are recomputed from the final leaf/sapwood state.
         result%sapwood_area_new = result%leaf_mass_new * params%sla / params%latosa

         if (result%sapwood_area_new > 0.0_real64) then
            result%height_new = result%sapwood_mass_new / &
                                 (params%wood_density * result%sapwood_area_new)
         else
            result%height_new = 0.0_real64
         end if

         if (result%height_new > 0.0_real64) then
            result%stem_diameter_new = (result%height_new / params%allom2)**(1.0_real64 / params%allom3)
         else
            result%stem_diameter_new = 0.0_real64
         end if

         result%carbon_balance_error = c_available - &
            (result%delta_leaf + result%delta_root + result%delta_sapwood)

         ! Note: in abnormal allocation, carbon balance over living tissues alone
         ! can look different because negative sapwood is converted to heartwood.
         ! For a full model integration, track litter and sapwood-to-heartwood
         ! fluxes explicitly outside this minimal kernel.

         result%leaf_root_residual = result%leaf_mass_new - &
            params%leaf_to_root_ratio * result%root_mass_new

         result%allocation_residual_final = allocation_residual(state, params, c_available, result%delta_leaf)

         result%message = "Abnormal allocation used. Track litter and sapwood-to-heartwood fluxes explicitly in integration."

      end subroutine abnormal_allocation

end module carbon_allocation_offline_kernel
