module carbon_allocation_offline_kernel

!==========================================================================
! Carbon allocation
!==========================================================================
! This file is intentionally verbose. The comments are part of the model
! documentation and are meant to help a future reader understand why each
! equation appears in the code. 
!
! Purpose
! -------
! This module implements a self-contained woody-plant carbon allocation
! kernel based on the LPJ-style allometric allocation logic. BUT important:
! Not using "abnormal allocation" as LPJ-style. If there is no feasible solution 
! to the normal allometric problem, the carbon goes to storage instead of being 
! forced into a non-allometric solution. This is a more
! conservative approach that avoids unrealistic jumps in plant structure
!
! Scope of this first version
! ---------------------------
! This kernel solves allocation for ONE average woody individual over ONE
! allocation period. The allocation period can be annual, monthly, seasonal,
! or any other interval. The equations do not know the calendar frequency.
! The caller only needs to provide the carbon available over that period (here
! we are using daily because CAETE works on a daily basis).
!
! Key idea for the gradual storage-based allocation scheme
! -------------------------------------------------------
! Daily carbon input is first added to a labile storage pool. Structural
! growth is then paid from storage only when there is positive allocation
! demand and enough available carbon.
!
! In contrast to the legacy rigid allocation scheme (strictly based in LPJ), this routine does not
! force the plant to satisfy allometric constraints exactly at each daily
! timestep. Instead, the leaf-root relationship and the pipe-model relationship
! are used to compute structural demand. This demand gradually moves the plant
! toward allometric consistency over the timescale defined by
! allometric_adjustment_days.
!
! The structural increment is:
!
!     structural_growth = dL + dR + dS
!
! where:
!
!     dL = increment to leaf carbon mass
!     dR = increment to fine-root carbon mass
!     dS = increment to sapwood carbon mass
!
! If there is no structural demand, or if daily allocation is limited by the
! maximum allocation fraction, carbon remains in storage instead of being
! forced into an abnormal allocation pathway.
!
! Height is updated from total stem carbon rather than by forcing the
! pipe-model residual to zero at each daily timestep. The pipe model therefore
! acts as a gradual demand signal, not as an instantaneous constraint.
  

   use, intrinsic :: iso_fortran_env, only: real64

   implicit none

   private   

   public :: Parameters
   public :: PlantCarbonState
   public :: AllocationOutput
   public :: allocation_residual
   public :: leaf_requirement
   public :: ControlsParam
   public :: allocate_gradual_with_storage


  !--------------------------------------------------------------------------
  ! Numerical constants
  !--------------------------------------------------------------------------

   real(real64), parameter :: pi = 3.1416_real64

   ! Tolerance used only for diagnostic carbon-accounting checks.
   real(real64), parameter :: carbon_accounting_tolerance = 1.0e-10_real64

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
  ! Controls for the gradual daily allocation routine with labile storage
  !--------------------------------------------------------------------------

   type :: ControlsParam

      ! Length of one model time step expressed in years.
      ! For a daily time step, use 1/365.
      real(real64) :: dt_years = 1.0_real64 / 365.0_real64

      ! Time scale used to relax allometric imbalances.
      ! A value of 365 means that a full allometric deficit is not corrected
      ! instantly; instead, roughly 1/365 of the deficit becomes demand per day.
      real(real64) :: allometric_adjustment_days = 365.0_real64

      ! Maximum fraction of current living structural carbon that can become new
      ! structural biomass in one time step. This avoids unrealistic daily jumps.
      ! This Upper bound limits daily structural growth and it is expressed as a fraction of the
      ! current living structural carbon pool:
      !
      !     living_carbon = leaf_mass + root_mass + sapwood_mass
      !
      ! This parameter only limits the maximum amount of new structural biomass 
      ! that can be produced in one timestep, even if storage carbon and structural 
      ! demand are both high. For example, max_allocation_fraction = 0.005 means that, in one daily
      ! timestep, structural growth cannot exceed 0.5% of the current living
      ! structural carbon. This prevents unrealistic jumps in leaf, root, or
      ! sapwood biomass when large storage pools or large allometric deficits
      ! are present.
      real(real64) :: max_allocation_fraction = 0.005_real64

      !!!! TO BE READJUSTED (can express fast/slow growth strategies)
      !---------------------
      ! Background demand time scale for leaves, in years.
      ! This is not an allometric target. It is a small baseline structural demand
      ! used when the plant is close to its allometric constraints.
      real(real64) :: leaf_background_timescale_years = 3.0_real64

      ! Background demand time scale for fine roots, in years.
      real(real64) :: root_background_timescale_years = 3.0_real64

      ! Background demand time scale for sapwood, in years.
      ! This should usually be longer than leaf and fine-root time scales.
      real(real64) :: sapwood_background_timescale_years = 15.0_real64

   end type ControlsParam



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


      ! Diagnostics specific to gradual allocation with labile carbon storage.
      ! These fields remain zero when the original rigid allocate() routine is used.
      real(real64) :: npp_daily = 0.0_real64
      real(real64) :: carbon_storage_before = 0.0_real64
      real(real64) :: carbon_storage_after = 0.0_real64
      real(real64) :: carbon_to_allocate = 0.0_real64
      real(real64) :: total_demand_daily = 0.0_real64
      real(real64) :: leaf_demand_daily = 0.0_real64
      real(real64) :: root_demand_daily = 0.0_real64
      real(real64) :: sapwood_demand_daily = 0.0_real64

      ! Carbon-accounting diagnostics for gradual allocation.
      ! storage_after_npp_unclamped stores the raw storage value after adding NPP.
      ! If this value is negative, storage is clamped to zero and the missing
      ! carbon is reported as unmet_storage_deficit.
      real(real64) :: storage_after_npp_unclamped = 0.0_real64
      real(real64) :: unmet_storage_deficit = 0.0_real64

      ! The structural balance checks whether all carbon assigned to growth was
      ! actually distributed among leaf, fine-root, and sapwood increments.
      real(real64) :: structural_increment_sum = 0.0_real64
      real(real64) :: structural_balance_error = 0.0_real64

      ! The storage balance checks whether the labile storage pool was updated
      ! consistently after adding NPP and subtracting structural allocation.
      real(real64) :: storage_balance_error = 0.0_real64

      ! The whole-plant balance checks the combined structural + storage carbon.
      ! If unmet_storage_deficit is zero, the whole plant should change by NPP.
      ! If unmet_storage_deficit is positive, the allocation routine did not have
      ! enough storage to pay the full negative NPP, and this deficit is reported.
      real(real64) :: whole_plant_balance_error = 0.0_real64

      ! Maximum structural carbon allocation allowed by the daily limiter.
      real(real64) :: max_daily_allocation = 0.0_real64

      ! True when the numerical carbon-accounting equations close within the
      ! tolerance defined by carbon_accounting_tolerance.
      logical :: carbon_accounting_ok = .false.

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
   !> Gradual daily allocation with labile carbon storage.
   !==========================================================================
      subroutine allocate_gradual_with_storage(state, params, controls, npp_rate, &
                                               carbon_storage, result)

         type(PlantCarbonState),          intent(in)    :: state
         type(Parameters),                intent(in)    :: params
         type(ControlsParam), intent(in)    :: controls
         real(real64),                    intent(in)    :: npp_rate
         real(real64),                    intent(inout) :: carbon_storage
         type(AllocationOutput),          intent(inout) :: result

         real(real64) :: storage_after_npp
         real(real64) :: living_carbon
         real(real64) :: max_daily_allocation

         real(real64) :: leaf_required
         real(real64) :: delta_leaf_min_pipe
         real(real64) :: delta_leaf_min_root_nonnegative
         real(real64) :: delta_leaf_min
         real(real64) :: leaf_target

         real(real64) :: root_required
         real(real64) :: delta_root_min
         real(real64) :: root_target

         real(real64) :: sapwood_required_for_leaf_target
         real(real64) :: sapwood_target

         real(real64) :: leaf_deficit
         real(real64) :: root_deficit
         real(real64) :: sapwood_deficit

         real(real64) :: leaf_background_demand
         real(real64) :: root_background_demand
         real(real64) :: sapwood_background_demand

         real(real64) :: frac_leaf
         real(real64) :: frac_root
         real(real64) :: frac_sapwood

         real(real64) :: leaf_area_new
         real(real64) :: sapwood_area_from_mass
         real(real64) :: stem_carbon_total_new
         real(real64) :: height_power_total_stem
         real(real64) :: height_power_exponent
         real(real64) :: pi_over_four

         real(real64) :: structural_carbon_before
         real(real64) :: structural_carbon_after
         real(real64) :: whole_carbon_before
         real(real64) :: whole_carbon_after

         ! Reset the output object to a known state.
         result = AllocationOutput()

         result%carbon_storage_before = carbon_storage

         ! Convert the annualized NPP rate into carbon input over this time step.
         ! If npp_rate is in kgC per area per year and dt_years is 1/365,
         ! npp_daily has units of kgC per area per day.
         result%npp_daily = npp_rate * controls%dt_years

         ! Update the labile carbon storage. Negative NPP consumes storage.
         ! The unclamped value is kept for carbon-accounting diagnostics.
         result%storage_after_npp_unclamped = carbon_storage + result%npp_daily
         storage_after_npp = result%storage_after_npp_unclamped

         ! Storage cannot become negative. If NPP is strongly negative, the
         ! remaining deficit is not paid by this routine. We report that amount
         ! explicitly instead of hiding it inside the storage clamp.
         if (storage_after_npp < 0.0_real64) then
            result%unmet_storage_deficit = -storage_after_npp
            storage_after_npp = 0.0_real64
         else
            result%unmet_storage_deficit = 0.0_real64
         end if

         !--------------------------------------------------------------------
         ! Allometric correction direction following the original allocation
         ! logic.
         !--------------------------------------------------------------------
         ! Important: there is only one primary leaf requirement here.
         ! It comes from the pipe model and is exactly the same quantity used
         ! in the rigid normal-allocation solver:
         !
         !     leaf_required = leaf_requirement(state, params)
         !
         ! This is the leaf mass needed to support the current sapwood mass at
         ! the current height.
         leaf_required = leaf_requirement(state, params)

         ! Minimum leaf increment required by the pipe model.
         ! If the current leaf mass is already above the pipe-model requirement,
         ! this term is negative and will not create demand.
         delta_leaf_min_pipe = leaf_required - state%leaf_mass

         ! Minimum leaf increment required to avoid a negative root increment
         ! when the leaf-root relationship is imposed as:
         !
         !     root_new = leaf_new / leaf_to_root_ratio
         !
         ! In the rigid solver, root increment is computed as:
         !
         !     delta_root = (leaf_old + delta_leaf) / leaf_to_root_ratio
         !                  - root_old
         !
         ! To keep delta_root >= 0, delta_leaf must satisfy:
         !
         !     delta_leaf >= root_old * leaf_to_root_ratio - leaf_old
         !
         ! This is not a second independent leaf requirement. It is a numerical
         ! and biological constraint that prevents negative root allocation in
         ! the normal-growth pathway.
         delta_leaf_min_root_nonnegative = &
            state%root_mass * params%leaf_to_root_ratio - state%leaf_mass

         ! Effective minimum leaf increment implied by the same lower-bound
         ! logic used in the original normal allocation routine.
         ! This is a deficit used as a direction of gradual adjustment, not an
         ! increment that must be applied in a single daily time step.
         delta_leaf_min = max(0.0_real64, &
                              delta_leaf_min_pipe, &
                              delta_leaf_min_root_nonnegative)

         ! Leaf target associated with the minimum admissible leaf increment.
         leaf_target = state%leaf_mass + delta_leaf_min

         ! Root target is derived after the leaf target, following the original
         ! logic:
         !
         !     root_required = leaf_target / leaf_to_root_ratio
         !
         ! This means the pipe-model leaf requirement comes first, and the root
         ! requirement follows from functional balance.
         if (params%leaf_to_root_ratio > 0.0_real64) then
            root_required = leaf_target / params%leaf_to_root_ratio
         else
            root_required = state%root_mass
         end if

         delta_root_min = max(0.0_real64, root_required - state%root_mass)
         root_target = state%root_mass + delta_root_min

         ! Sapwood is different from leaf and root in the original rigid solver:
         ! it is not prescribed by an independent target. In normal allocation,
         ! sapwood increment is the remaining carbon after leaf and root:
         !
         !     delta_sapwood = c_available - delta_leaf - delta_root
         !
         ! For the gradual daily routine, however, we still need a sapwood demand
         ! direction. We therefore compute the sapwood mass that would support
         ! the current leaf target under the pipe model at the current height.
         ! This is a gradual target only; it is not forced instantly.
         sapwood_required_for_leaf_target = leaf_target * params%sla / params%latosa * &
                                            params%wood_density * state%height

         sapwood_target = max(state%sapwood_mass, sapwood_required_for_leaf_target)

         ! Positive deficits define the allometric correction direction.
         leaf_deficit = max(0.0_real64, leaf_target - state%leaf_mass)
         root_deficit = max(0.0_real64, root_target - state%root_mass)
         sapwood_deficit = max(0.0_real64, sapwood_target - state%sapwood_mass)

         ! Convert full allometric deficits into gradual demands. This is the key
         ! step that prevents instant correction of the whole plant structure.
         result%leaf_demand_daily = leaf_deficit / controls%allometric_adjustment_days
         result%root_demand_daily = root_deficit / controls%allometric_adjustment_days
         result%sapwood_demand_daily = sapwood_deficit / controls%allometric_adjustment_days

         !--------------------------------------------------------------------
         ! Small background structural demand.
         !--------------------------------------------------------------------
         ! These terms allow growth to continue when the plant is already close
         ! to the allometric correction targets. They are deliberately separated
         ! from the allometric correction terms above.
         !
         ! The background terms are optional and should be tested carefully. If
         ! they are set too high, they can dominate the allometric correction
         ! signal. If they are set to zero, allocation occurs only when there is
         ! an explicit allometric deficit.
         if (controls%leaf_background_timescale_years > 0.0_real64) then
            leaf_background_demand = leaf_target * controls%dt_years / &
                                     controls%leaf_background_timescale_years
         else
            leaf_background_demand = 0.0_real64
         end if

         if (controls%root_background_timescale_years > 0.0_real64) then
            root_background_demand = root_target * controls%dt_years / &
                                     controls%root_background_timescale_years
         else
            root_background_demand = 0.0_real64
         end if

         if (controls%sapwood_background_timescale_years > 0.0_real64) then
            sapwood_background_demand = sapwood_target * controls%dt_years / &
                                        controls%sapwood_background_timescale_years
         else
            sapwood_background_demand = 0.0_real64
         end if

         result%leaf_demand_daily = result%leaf_demand_daily + leaf_background_demand
         result%root_demand_daily = result%root_demand_daily + root_background_demand
         result%sapwood_demand_daily = result%sapwood_demand_daily + sapwood_background_demand

         result%total_demand_daily = result%leaf_demand_daily + &
                                      result%root_demand_daily + &
                                      result%sapwood_demand_daily

         ! The daily structural allocation is constrained by available storage,
         ! structural demand, and a maximum allowed growth fraction.
         living_carbon = state%leaf_mass + state%root_mass + state%sapwood_mass
         max_daily_allocation = controls%max_allocation_fraction * living_carbon
         result%max_daily_allocation = max_daily_allocation

         if (result%total_demand_daily > 0.0_real64 .and. &
             storage_after_npp > 0.0_real64 .and. &
             max_daily_allocation > 0.0_real64) then

            result%carbon_to_allocate = min(storage_after_npp, &
                                            result%total_demand_daily, &
                                            max_daily_allocation)

            frac_leaf = result%leaf_demand_daily / result%total_demand_daily
            frac_root = result%root_demand_daily / result%total_demand_daily
            frac_sapwood = result%sapwood_demand_daily / result%total_demand_daily

            result%delta_leaf = frac_leaf * result%carbon_to_allocate
            result%delta_root = frac_root * result%carbon_to_allocate
            result%delta_sapwood = frac_sapwood * result%carbon_to_allocate

         else

            result%carbon_to_allocate = 0.0_real64
            result%delta_leaf = 0.0_real64
            result%delta_root = 0.0_real64
            result%delta_sapwood = 0.0_real64

         end if

         ! Update structural carbon pools.
         result%leaf_mass_new = state%leaf_mass + result%delta_leaf
         result%root_mass_new = state%root_mass + result%delta_root
         result%sapwood_mass_new = state%sapwood_mass + result%delta_sapwood
         result%heartwood_mass_new = state%heartwood_mass

         ! Update storage after structural allocation.
         carbon_storage = storage_after_npp - result%carbon_to_allocate
         result%carbon_storage_after = carbon_storage

         ! For gradual allocation, height is computed from total stem carbon and
         ! height-diameter allometry. This avoids forcing the pipe model to be
         ! exactly satisfied in one daily step.
         stem_carbon_total_new = result%sapwood_mass_new + result%heartwood_mass_new
         height_power_exponent = 1.0_real64 + 2.0_real64 / params%allom3
         pi_over_four = pi / 4.0_real64

         if (stem_carbon_total_new > 0.0_real64 .and. params%wood_density > 0.0_real64) then
            height_power_total_stem = params%allom2**(2.0_real64 / params%allom3) * &
                                      (stem_carbon_total_new / params%wood_density) / &
                                      pi_over_four
            result%height_new = height_power_total_stem**(1.0_real64 / height_power_exponent)
         else
            result%height_new = state%height
         end if

         if (result%height_new > 0.0_real64) then
            result%stem_diameter_new = (result%height_new / params%allom2)**(1.0_real64 / params%allom3)
         else
            result%stem_diameter_new = 0.0_real64
         end if

         result%sapwood_area_new = result%leaf_mass_new * params%sla / params%latosa

         ! Diagnostics.
         ! Structural carbon accounting: carbon sent to growth must equal the
         ! sum of the structural increments.
         result%structural_increment_sum = result%delta_leaf + &
            result%delta_root + result%delta_sapwood

         result%structural_balance_error = result%carbon_to_allocate - &
            result%structural_increment_sum

         ! Keep the older diagnostic name for backward compatibility.
         result%carbon_balance_error = result%structural_balance_error

         ! Storage carbon accounting: final storage must equal initial storage
         ! plus NPP input, plus any reported unmet deficit caused by the zero
         ! lower bound, minus the carbon allocated to structure.
         result%storage_balance_error = result%carbon_storage_after - &
            (result%carbon_storage_before + result%npp_daily + &
             result%unmet_storage_deficit - result%carbon_to_allocate)

         ! Whole-plant carbon accounting: structural + storage carbon should
         ! change by NPP, except when negative NPP exceeds available storage.
         structural_carbon_before = state%leaf_mass + state%root_mass + &
            state%sapwood_mass + state%heartwood_mass

         structural_carbon_after = result%leaf_mass_new + result%root_mass_new + &
            result%sapwood_mass_new + result%heartwood_mass_new

         whole_carbon_before = structural_carbon_before + result%carbon_storage_before
         whole_carbon_after = structural_carbon_after + result%carbon_storage_after

         result%whole_plant_balance_error = (whole_carbon_after - whole_carbon_before) - &
            (result%npp_daily + result%unmet_storage_deficit)

         result%carbon_accounting_ok = &
            abs(result%structural_balance_error) <= carbon_accounting_tolerance .and. &
            abs(result%storage_balance_error) <= carbon_accounting_tolerance .and. &
            abs(result%whole_plant_balance_error) <= carbon_accounting_tolerance

         result%leaf_root_residual = result%leaf_mass_new - &
            params%leaf_to_root_ratio * result%root_mass_new

         leaf_area_new = result%leaf_mass_new * params%sla
         sapwood_area_from_mass = 0.0_real64
         if (result%height_new > 0.0_real64) then
            sapwood_area_from_mass = result%sapwood_mass_new / &
                                     (params%wood_density * result%height_new)
         end if

         result%pipe_model_residual = leaf_area_new - params%latosa * sapwood_area_from_mass
         result%allocation_residual_final = 0.0_real64

         result%message = "Gradual storage allocation used; allometry follows the original lower-bound logic."

      end subroutine allocate_gradual_with_storage

end module carbon_allocation_offline_kernel
