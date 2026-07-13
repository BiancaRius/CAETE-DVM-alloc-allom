!===============================================================================
! Offline driver for the carbon-allocation kernel.
!
! Purpose
! -------
! This program runs controlled, daily carbon-allocation experiments outside the
! full vegetation model. It is intended to test the numerical behavior and
! carbon accounting of the allocation routine across combinations of NPP,
! initial labile storage, allometric adjustment speed, construction capacity,
! background structural demand, and plant traits.
!
! For each scenario, the driver:
!   1. initializes plant traits, allocation controls, and carbon pools;
!   2. calls allocate_gradual_with_storage once per simulated day;
!   3. checks pool non-negativity, finite values, allocation limits, and carbon
!      balance errors;
!   4. accumulates growth, starvation, turnover, and storage diagnostics; and
!   5. writes summary, checkpoint, and selected time-series CSV files.
!
! Important scope
! ---------------
! NPP is prescribed as a constant annual rate within each scenario. This driver
! does not simulate climate, photosynthesis, or respiration; those processes
! would determine the NPP supplied to the allocation kernel in the full model.
! All carbon variables use the same mass units as carbon_allocation_offline_kernel.
!===============================================================================
! ATTENTION: this is the complete version of the file, with all the cheks of erros and
! carbon balance. The other version is a simplified one (run_carbon_allocation_offline.f90), 
!with less checks and more focused on the allocation itself.

program run_carbon_allocation_offline

   use, intrinsic :: iso_fortran_env, only: real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use carbon_allocation_offline_kernel
 
   implicit none
 
   ! Simulation length. The current allocation controls assume one call per day.
   integer, parameter :: n_days = 365 * 10

   ! Acceptance tolerance for non-negativity, limiter, and carbon-balance tests.
   real(real64), parameter :: tol = 1.0e-10_real64

   ! Threshold used to classify very small fluxes as zero and to match selected
   ! floating-point scenario values when deciding which traces to write.
   real(real64), parameter :: tiny_positive = 1.0e-12_real64
   

   !________________________________________________________________________
   ! Testing different scenario factors. Each factor is represented by an array of values

   ! Number of levels in each scenario factor. These dimensions must match the
   ! corresponding value arrays declared below (for example, if you want to test 
   ! 10 values of npp, n_npp should be 10). 

   ! Number of NPP levels included in the scenario matrix.
   ! This factor tests how carbon supply affects storage accumulation,
   ! structural growth, starvation risk, and allometric adjustment.
   ! It can be also used to test a theoretical fluctuation in NPP through time
   integer, parameter :: n_npp = 1

   ! Number of initial storage levels included in the scenario matrix.
   ! This factor tests how the amount of labile carbon available at the
   ! beginning of the simulation affects growth, buffering of negative NPP,
   ! and the ability to satisfy structural demand.
   integer, parameter :: n_storage = 1

   ! Number of allometric adjustment timescales included in the scenario matrix.
   ! This factor tests how quickly the plant attempts to correct imbalances
   ! among leaf, fine-root, and sapwood biomass.
   integer, parameter :: n_adjustment = 1

   ! Number of maximum allocation fractions included in the scenario matrix.
   ! This factor tests the effect of the daily construction-capacity limit,
   ! that is, the maximum fraction of living structural carbon that can be
   ! converted into new biomass during one timestep.
   integer, parameter :: n_max_fraction = 1

   ! Number of background-growth configurations included in the scenario matrix.
   ! This factor tests whether continued structural demand is allowed when the
   ! plant is already close to its allometric targets.
   integer, parameter :: n_background = 1

   ! Number of trait combinations included in the scenario matrix.
   ! This factor tests how alternative plant strategies, such as differences
   ! in SLA or wood density, modify allometric requirements and allocation.
   integer, parameter :: n_trait_case = 1

   ! Number of storage-turnover values intended for the scenario matrix.
   ! This factor is meant to test how quickly labile storage is lost through
   ! turnover. In the current implementation, however, these values are not yet
   ! passed to the allocation module, so changing this dimension alone does not
   ! modify the storage-turnover rate actually used by the model.
   integer, parameter :: n_turnover_storage = 1

   ! Container for scenario metadata, end states, cumulative fluxes, numerical
   ! diagnostics, and event counts. Quantities labelled "cumulative" are sums
   ! across all completed daily timesteps.
   type :: ScenarioSummary
 
      ! Scenario identifiers.
      integer :: scenario_id = 0
      logical :: passed = .true.
      integer :: fail_day = 0
      character(len=96) :: failure_reason = "OK"
 
      ! Scenario factors.
      integer :: trait_case = 0
      integer :: background_mode = 0
      real(real64) :: npp_rate = 0.0_real64
      real(real64) :: initial_storage = 0.0_real64
      real(real64) :: allometric_adjustment_days = 0.0_real64
      real(real64) :: max_allocation_fraction = 0.0_real64
 
      ! Initial state. "Living carbon" includes leaf, fine root, and sapwood.
      ! "Structural carbon" also includes heartwood, and "total carbon" adds
      ! labile storage to structural carbon.
      real(real64) :: initial_leaf = 0.0_real64
      real(real64) :: initial_root = 0.0_real64
      real(real64) :: initial_sapwood = 0.0_real64
      real(real64) :: initial_heartwood = 0.0_real64
      real(real64) :: initial_height = 0.0_real64
      real(real64) :: initial_living_carbon = 0.0_real64
      real(real64) :: initial_structural_carbon = 0.0_real64
      real(real64) :: initial_total_carbon = 0.0_real64
 
      ! Final state and net changes relative to the initialized plant. The same
      ! living, structural, and total-carbon definitions are used as above.
      real(real64) :: final_leaf = 0.0_real64
      real(real64) :: final_root = 0.0_real64
      real(real64) :: final_sapwood = 0.0_real64
      real(real64) :: final_heartwood = 0.0_real64
      real(real64) :: final_height = 0.0_real64
      real(real64) :: final_storage = 0.0_real64
      real(real64) :: final_living_carbon = 0.0_real64
      real(real64) :: final_structural_carbon = 0.0_real64
      real(real64) :: final_total_carbon = 0.0_real64
      real(real64) :: net_living_carbon_change = 0.0_real64
      real(real64) :: net_total_carbon_change = 0.0_real64
 
      ! Integrated carbon input and construction. Positive NPP is tracked
      ! separately because structural allocation can also be financed by initial
      ! storage; therefore, allocation / positive NPP is only a diagnostic ratio.
      real(real64) :: cumulative_npp = 0.0_real64
      real(real64) :: cumulative_positive_npp = 0.0_real64
      real(real64) :: cumulative_structural_allocation = 0.0_real64
 
      ! Integrated starvation diagnostics. An unmet storage deficit is the
      ! carbon shortfall created when negative NPP exceeds available storage.
      ! Leaf and root starvation losses leave the plant carbon system, whereas
      ! starved sapwood is transferred to heartwood in the current kernel.
      real(real64) :: cumulative_unmet_storage_deficit = 0.0_real64
      real(real64) :: cumulative_leaf_starvation_loss = 0.0_real64
      real(real64) :: cumulative_root_starvation_loss = 0.0_real64
      real(real64) :: cumulative_sapwood_starvation_loss = 0.0_real64
      real(real64) :: cumulative_sapwood_to_heartwood_starvation = 0.0_real64
      real(real64) :: cumulative_starvation_carbon_loss = 0.0_real64
      real(real64) :: cumulative_unpaid_carbon_deficit = 0.0_real64
 
      ! Integrated turnover diagnostics. Sapwood turnover is transferred to
      ! heartwood and is therefore not counted as an immediate whole-plant carbon
      ! loss; turnover losses from other pools leave the tracked plant system.
      real(real64) :: cumulative_leaf_turnover_loss = 0.0_real64
      real(real64) :: cumulative_root_turnover_loss = 0.0_real64
      real(real64) :: cumulative_sapwood_turnover_loss = 0.0_real64
      real(real64) :: cumulative_storage_turnover_loss = 0.0_real64
      real(real64) :: cumulative_heartwood_turnover_loss = 0.0_real64
      real(real64) :: cumulative_sapwood_to_heartwood_turnover = 0.0_real64
      real(real64) :: cumulative_total_sapwood_to_heartwood = 0.0_real64
      real(real64) :: cumulative_turnover_carbon_loss = 0.0_real64
 
      ! Final allometric and storage diagnostics. Residuals are biological
      ! disequilibrium metrics, not numerical errors, and need not equal zero at
      ! each daily timestep in a gradual-allocation scheme.
      real(real64) :: final_leaf_root_residual = 0.0_real64
      real(real64) :: final_pipe_residual = 0.0_real64
      real(real64) :: final_storage_fraction_total = 0.0_real64
      real(real64) :: final_storage_fraction_living = 0.0_real64
      real(real64) :: max_storage_fraction_total = 0.0_real64
      real(real64) :: max_storage_fraction_living = 0.0_real64
      real(real64) :: structural_fraction_of_positive_npp = 0.0_real64
 
      ! Numerical diagnostics and observed extrema used to detect instability,
      ! bookkeeping errors, negative pools, or unintended allocation spikes.
      real(real64) :: max_abs_structural_balance_error = 0.0_real64
      real(real64) :: max_abs_storage_balance_error = 0.0_real64
      real(real64) :: max_abs_whole_plant_balance_error = 0.0_real64
      real(real64) :: min_living_pool = huge(1.0_real64)
      real(real64) :: min_living_carbon = huge(1.0_real64)
      real(real64) :: max_daily_allocation_observed = 0.0_real64
      real(real64) :: max_daily_demand_observed = 0.0_real64
 
      ! Event counters quantify how often key processes or storage thresholds
      ! occurred during the scenario, rather than only reporting final values.
      integer :: days_with_allocation = 0
      integer :: days_with_unmet_storage_deficit = 0
      integer :: days_with_turnover = 0
      integer :: days_with_storage_above_10pct_living = 0
      integer :: days_with_storage_above_50pct_living = 0
      integer :: days_with_storage_above_100pct_living = 0
 
   end type ScenarioSummary
 
   ! Scenario-factor arrays. Scalar assignment below broadcasts the same value
   ! to every element; array constructors can be used when multiple levels are
   ! enabled. NPP rates are annualized inputs and are scaled by dt_years inside
   ! the allocation kernel.
   real(real64), dimension(n_npp) :: npp_values
   real(real64), dimension(n_storage) :: storage_values
   real(real64), dimension(n_adjustment) :: adjustment_values
   real(real64), dimension(n_max_fraction) :: max_fraction_values
   real(real64), dimension(n_turnover_storage) :: turnover_sto_values

   integer :: i_npp
   integer :: i_storage
   integer :: i_adjustment
   integer :: i_max_fraction
   integer :: i_background
   integer :: i_trait_case
   !integer :: i_turnover_storage !uncomment if you want to test different storage turnover values
   integer :: scenario_id
   integer :: n_pass
   integer :: n_fail
   integer :: summary_unit
   integer :: checkpoint_unit
   integer :: daily_unit
 
   type(Parameters) :: params
   type(ControlsParam) :: controls
   type(PlantCarbonState) :: initial_state
   type(ScenarioSummary) :: summary
 
   ! Active scenario values. Examples of broader sensitivity grids are retained
   ! below each factor for convenient reactivation; update the corresponding n_*
   ! dimension whenever an array constructor contains more than one value.
   npp_values = 3.5_real64
   ! Example: npp_values = [ -0.5_real64, 0.0_real64, 0.5_real64, &
   !                          3.5_real64, 8.0_real64 ]

   storage_values = 0.5_real64
   ! Example: storage_values = [ 0.0_real64, 0.5_real64, 5.0_real64 ]

   adjustment_values = 365.0_real64
   ! Example: adjustment_values = [ 30.0_real64, 365.0_real64, 730.0_real64 ]

   max_fraction_values = 0.005_real64
   ! Example: max_fraction_values = [ 0.001_real64, 0.005_real64, &
   !                                  0.02_real64 ]

   ! IMPORTANT: this array is currently a placeholder only. Its selected value is
   ! not passed to ControlsParam or to allocate_gradual_with_storage because the
   ! storage-turnover rate is defined inside the kernel module. Increasing
   ! n_turnover_storage in the current code would therefore repeat identical
   ! scenarios rather than perform a true storage-turnover sensitivity analysis.
   ! turnover_sto_values = 0.05_real64
   ! Intended future grid: [ 0.0, 0.05, 0.10, 0.20, 0.50 ] yr^-1.

   scenario_id = 0
   n_pass = 0
   n_fail = 0

   ! Create new output files for this run. status="replace" intentionally
   ! overwrites files with the same names in the current working directory.
   open(newunit=summary_unit, &
        file="storage_allocation_sensitivity_summary_turnover.csv", &
        status="replace", action="write")
   call write_summary_header(summary_unit)
 
   open(newunit=checkpoint_unit, &
        file="storage_allocation_sensitivity_checkpoints_turnover.csv", &
        status="replace", action="write")
   call write_checkpoint_header(checkpoint_unit)
 
   open(newunit=daily_unit, &
        file="storage_allocation_daily_selected_turnover.csv", &
        status="replace", action="write")
   call write_daily_header(daily_unit)
 
   ! Evaluate every combination of the enabled scenario factors. The order of
   ! the loops only controls scenario numbering and output order.
   do i_trait_case = 1, n_trait_case
      do i_background = 1, n_background
         do i_max_fraction = 1, n_max_fraction
            do i_adjustment = 1, n_adjustment
               do i_storage = 1, n_storage
                  do i_npp = 1, n_npp
                     
                     ! use this loop to test different storage turnover values in the future
                     !do i_turnover_storage = 1, n_turnover_storage

                        scenario_id = scenario_id + 1
   
                        ! Build all inputs for one independent scenario.
                        call initialize_parameters(params, i_trait_case)
                        call initialize_controls(controls, &
                           adjustment_values(i_adjustment), &
                           max_fraction_values(i_max_fraction), &
                           i_background)
                        call initialize_state(params, initial_state)
   
                        ! Run the complete daily simulation and return an
                        ! aggregated validation and carbon-accounting summary.
                        call run_scenario(scenario_id, params, controls, &
                           initial_state, npp_values(i_npp), &
                           storage_values(i_storage), i_trait_case, &
                           i_background, checkpoint_unit, &
                           daily_unit, summary)

                        call write_summary_row(summary_unit, summary)

                        if (summary%passed) then
                           n_pass = n_pass + 1
                        else
                           n_fail = n_fail + 1
                        end if
                     !end do
                  end do
               end do
            end do
         end do
      end do
   end do
 
   close(summary_unit)
   close(checkpoint_unit)
   close(daily_unit)
 
   write(*,'(a,i0)') "Sensitivity scenarios completed: ", scenario_id
   write(*,'(a,i0)') "PASS: ", n_pass
   write(*,'(a,i0)') "FAIL: ", n_fail
   write(*,'(a)') "Output written to storage_allocation_sensitivity_summary_turnover.csv"
   write(*,'(a)') "Checkpoint output written to storage_allocation_sensitivity_checkpoints_turnover.csv"
   write(*,'(a)') "Selected daily output written to storage_allocation_daily_selected_turnover.csv"
 
   if (n_fail > 0) then
      error stop "At least one sensitivity scenario failed. Inspect the CSV file."
   end if
 
 contains
 
   !---------------------------------------------------------------------------
   ! Initialize the trait and allometric parameter set for one strategy.
   !
   ! The base case is assigned first, after which trait_case changes one focal
   ! trait while keeping all remaining parameters fixed. This design isolates
   ! the effect of the selected trait in sensitivity tests.
   !---------------------------------------------------------------------------
   subroutine initialize_parameters(params, trait_case)
 
     type(Parameters), intent(out) :: params
     integer, intent(in) :: trait_case
 
     ! Assign the common baseline first. SLA controls leaf area per unit leaf
     ! carbon; LATOSA defines the target leaf-area:sapwood-area relationship;
     ! wood_density converts stem carbon to volume; leaf_to_root_ratio defines
     ! the target mass balance; and allom2/allom3 define height-diameter allometry.
     params%sla = 12.0_real64
     params%latosa = 8000.0_real64
     params%wood_density = 250.0_real64
     params%leaf_to_root_ratio = 1.0_real64
     params%allom2 = 40.0_real64
     params%allom3 = 0.5_real64
 
     select case (trait_case)
     case (1)
        ! Base trait combination.
     case (2)
        ! Low-SLA strategy.
        params%sla = 6.0_real64
     case (3)
        ! High-SLA strategy.
        params%sla = 24.0_real64
     case (4)
        ! Denser wood strategy.
        params%wood_density = 500.0_real64
     case default
        error stop "Unknown trait case."
     end select
 
   end subroutine initialize_parameters
 
 
   !---------------------------------------------------------------------------
   ! Initialize timestep and allocation-control parameters.
   !
   ! adjustment_days controls the daily fraction of the current allometric
   ! deficit requested for correction. max_fraction limits daily structural
   ! construction relative to initial living carbon. background_mode switches
   ! the optional baseline structural sink on or off.
   !---------------------------------------------------------------------------
   subroutine initialize_controls(controls, adjustment_days, max_fraction, &
                                  background_mode)
 
     type(ControlsParam), intent(out) :: controls
     real(real64), intent(in) :: adjustment_days
     real(real64), intent(in) :: max_fraction
     integer, intent(in) :: background_mode
 
     ! One kernel call represents one day. NPP, turnover, and background demand
     ! are scaled by this fraction of a year.
     controls%dt_years = 1.0_real64 / 365.0_real64
     ! At each daily call, correction demand equals the current allometric
     ! deficit divided by adjustment_days. This is a relaxation timescale, not a
     ! guarantee that the complete deficit is removed within exactly that period.
     controls%allometric_adjustment_days = adjustment_days

     ! Maximum fraction of initial living structural carbon that can be converted
     ! into new leaf, fine-root, and sapwood biomass during one daily call.
     controls%max_allocation_fraction = max_fraction
 
     select case (background_mode)
     case (1)
        ! Enable a baseline structural sink so that an allometrically balanced
        ! plant can continue growing when carbon is available. Shorter timescales
        ! produce stronger background demand; sapwood is assigned a slower rate.
        controls%leaf_background_timescale_years = 3.0_real64
        controls%root_background_timescale_years = 3.0_real64
        controls%sapwood_background_timescale_years = 15.0_real64
     case (2)
        ! Disable baseline sink demand. Structural construction then occurs only
        ! when positive allometric correction deficits are present.
        controls%leaf_background_timescale_years = 0.0_real64
        controls%root_background_timescale_years = 0.0_real64
        controls%sapwood_background_timescale_years = 0.0_real64
     case default
        error stop "Unknown background mode."
     end select
 
   end subroutine initialize_controls
 
 
   !---------------------------------------------------------------------------
   ! Construct the common initial plant state.
   !
   ! Leaf and fine-root masses are initialized in exact leaf-root balance. The
   ! sapwood mass is then solved numerically so that the initial state is also
   ! approximately consistent with the pipe-model relationship at the height
   ! implied by total stem carbon.
   !---------------------------------------------------------------------------
   subroutine initialize_state(params, state)
 
     type(Parameters), intent(in) :: params
     type(PlantCarbonState), intent(out) :: state
 
     real(real64) :: leaf_mass
     real(real64) :: root_mass
     real(real64) :: sapwood_mass
     real(real64) :: heartwood_mass

     ! Start from a state that is leaf-root balanced by construction. Sapwood is
     ! solved below because height depends on total stem carbon, which makes the
     ! pipe-balance equation implicit in sapwood mass.
     leaf_mass = 1.0_real64
     root_mass = leaf_mass / params%leaf_to_root_ratio
     heartwood_mass = 20.0_real64
     sapwood_mass = solve_sapwood_for_pipe_balance(params, leaf_mass, &
                                                   heartwood_mass)
 
     state = build_state(params, leaf_mass, root_mass, sapwood_mass, &
                         heartwood_mass)
 
   end subroutine initialize_state
 
 
   !---------------------------------------------------------------------------
   ! Assemble a PlantCarbonState from explicit pool masses and diagnose height
   ! from total stem carbon (sapwood + heartwood).
   !---------------------------------------------------------------------------
   function build_state(params, leaf_mass, root_mass, sapwood_mass, &
                        heartwood_mass) result(state)
 
     type(Parameters), intent(in) :: params
     real(real64), intent(in) :: leaf_mass
     real(real64), intent(in) :: root_mass
     real(real64), intent(in) :: sapwood_mass
     real(real64), intent(in) :: heartwood_mass
     type(PlantCarbonState) :: state
 
     state%leaf_mass = leaf_mass
     state%root_mass = root_mass
     state%sapwood_mass = sapwood_mass
     state%heartwood_mass = heartwood_mass
     state%height = height_from_total_stem_carbon(params, &
        sapwood_mass + heartwood_mass)
 
   end function build_state
 
 
   !---------------------------------------------------------------------------
   ! Run one complete scenario for n_days daily timesteps.
   !
   ! carbon_storage is updated in place by the allocation kernel. At each step,
   ! the routine validates the result before committing the new structural state,
   ! then accumulates fluxes and writes selected diagnostics.
   !---------------------------------------------------------------------------
   subroutine run_scenario(scenario_id, params, controls, initial_state, &
                           npp_rate, initial_storage, trait_case, &
                           background_mode, checkpoint_unit, daily_unit, summary)
 
     integer, intent(in) :: scenario_id
     type(Parameters), intent(in) :: params
     type(ControlsParam), intent(in) :: controls
     type(PlantCarbonState), intent(in) :: initial_state
     real(real64), intent(in) :: npp_rate
     real(real64), intent(in) :: initial_storage
     integer, intent(in) :: trait_case
     integer, intent(in) :: background_mode
     integer, intent(in) :: checkpoint_unit
     integer, intent(in) :: daily_unit
     type(ScenarioSummary), intent(out) :: summary
 
     type(PlantCarbonState) :: state
     type(AllocationOutput) :: result
     real(real64) :: carbon_storage
     real(real64) :: storage_fraction_total_now
     real(real64) :: storage_fraction_living_now
     integer :: day
 
     call initialize_summary(summary, scenario_id, params, controls, &
        initial_state, npp_rate, initial_storage, trait_case, &
        background_mode)
 
     state = initial_state
     carbon_storage = initial_storage
 
     do day = 1, n_days

        ! Run one daily allocation step. carbon_storage is an inout argument and
        ! is therefore already updated when the subroutine returns. The structural
        ! pools remain in result until the step has passed all validation checks.
        call allocate_gradual_with_storage(state, params, controls, npp_rate, &
                                           carbon_storage, result)
 
        ! Record extrema even for a step that subsequently fails, so the output
        ! preserves the numerical signal associated with the failure.
        call update_summary_diagnostics(summary, result)

        ! Stop this scenario at the first invalid timestep. Later scenarios still
        ! run, and the final program status reports the total number of failures.
        if (.not. step_is_valid(state, result, carbon_storage, &
                                summary%failure_reason)) then
           summary%passed = .false.
           summary%fail_day = day
           exit
        end if
 
        ! Only validated fluxes contribute to cumulative diagnostics.
        call accumulate_daily_fluxes(summary, result)

        ! Commit the validated post-allocation, post-starvation, and post-turnover
        ! structural state for use as the next day's input.
        state%leaf_mass = result%leaf_mass_new
        state%root_mass = result%root_mass_new
        state%sapwood_mass = result%sapwood_mass_new
        state%heartwood_mass = result%heartwood_mass_new
        state%height = result%height_new
 
        ! Diagnose storage relative to both total tracked carbon and living
        ! structural carbon after the daily state update.
        storage_fraction_total_now = storage_fraction_total(state, carbon_storage)
        storage_fraction_living_now = storage_fraction_living(state, carbon_storage)
 
        summary%max_storage_fraction_total = max(&
           summary%max_storage_fraction_total, storage_fraction_total_now)
        summary%max_storage_fraction_living = max(&
           summary%max_storage_fraction_living, storage_fraction_living_now)
 
        summary%min_living_pool = min(summary%min_living_pool, &
           min(state%leaf_mass, min(state%root_mass, state%sapwood_mass)))
        summary%min_living_carbon = min(summary%min_living_carbon, &
           living_carbon_of_state(state))
 
        if (result%carbon_to_allocate > tiny_positive) then
           summary%days_with_allocation = summary%days_with_allocation + 1
        end if
 
        if (result%unmet_storage_deficit > tiny_positive) then
           summary%days_with_unmet_storage_deficit = &
              summary%days_with_unmet_storage_deficit + 1
        end if
 
        if (result%turnover_carbon_loss > tiny_positive .or. &
            result%sapwood_to_heartwood_turnover > tiny_positive) then
           summary%days_with_turnover = summary%days_with_turnover + 1
        end if
 
        if (storage_fraction_living_now > 0.10_real64) then
           summary%days_with_storage_above_10pct_living = &
              summary%days_with_storage_above_10pct_living + 1
        end if
 
        if (storage_fraction_living_now > 0.50_real64) then
           summary%days_with_storage_above_50pct_living = &
              summary%days_with_storage_above_50pct_living + 1
        end if
 
        if (storage_fraction_living_now > 1.00_real64) then
           summary%days_with_storage_above_100pct_living = &
              summary%days_with_storage_above_100pct_living + 1
        end if
 
        ! Detailed output is intentionally restricted to selected scenarios to
        ! avoid generating an unnecessarily large daily file.
        if (should_write_daily_trace(summary, day)) then
           call write_daily_row(daily_unit, summary, params, state, result, &
                                carbon_storage, day)
        end if
 
        ! Checkpoints provide comparable state snapshots at fixed simulation ages.
        if (is_checkpoint_day(day)) then
           call write_checkpoint_row(checkpoint_unit, summary, params, state, &
                                     carbon_storage, day)
        end if
 
     end do
 
     call finalize_summary(summary, params, state, carbon_storage)
 
   end subroutine run_scenario
 
 
   !---------------------------------------------------------------------------
   ! Initialize all metadata and diagnostic baselines for a scenario.
   !
   ! Assignment from ScenarioSummary() resets every component to its declared
   ! default, preventing values from a previous scenario from being retained.
   !---------------------------------------------------------------------------
   subroutine initialize_summary(summary, scenario_id, params, controls, &
                                 initial_state, npp_rate, initial_storage, &
                                 trait_case, background_mode)
 
     type(ScenarioSummary), intent(out) :: summary
     integer, intent(in) :: scenario_id
     type(Parameters), intent(in) :: params
     type(ControlsParam), intent(in) :: controls
     type(PlantCarbonState), intent(in) :: initial_state
     real(real64), intent(in) :: npp_rate
     real(real64), intent(in) :: initial_storage
     integer, intent(in) :: trait_case
     integer, intent(in) :: background_mode
 
     ! Reset all scalar components to the default initializers declared in the
     ! derived type before assigning scenario-specific values.
     summary = ScenarioSummary()
 
     summary%scenario_id = scenario_id
     summary%passed = .true.
     summary%fail_day = 0
     summary%failure_reason = "OK"
 
     summary%trait_case = trait_case
     summary%background_mode = background_mode
     summary%npp_rate = npp_rate
     summary%initial_storage = initial_storage
     summary%allometric_adjustment_days = controls%allometric_adjustment_days
     summary%max_allocation_fraction = controls%max_allocation_fraction
 
     summary%initial_leaf = initial_state%leaf_mass
     summary%initial_root = initial_state%root_mass
     summary%initial_sapwood = initial_state%sapwood_mass
     summary%initial_heartwood = initial_state%heartwood_mass
     summary%initial_height = initial_state%height
     summary%initial_living_carbon = living_carbon_of_state(initial_state)
     summary%initial_structural_carbon = structural_carbon_of_state(initial_state)
     summary%initial_total_carbon = summary%initial_structural_carbon + &
                                    initial_storage
 
     summary%final_leaf = initial_state%leaf_mass
     summary%final_root = initial_state%root_mass
     summary%final_sapwood = initial_state%sapwood_mass
     summary%final_heartwood = initial_state%heartwood_mass
     summary%final_height = initial_state%height
     summary%final_storage = initial_storage
     summary%final_living_carbon = summary%initial_living_carbon
     summary%final_structural_carbon = summary%initial_structural_carbon
     summary%final_total_carbon = summary%initial_total_carbon
 
     summary%final_leaf_root_residual = leaf_root_residual_for_state(params, &
                                                                     initial_state)
     summary%final_pipe_residual = pipe_residual_for_state(params, initial_state)
     summary%final_storage_fraction_total = storage_fraction_total(initial_state, &
                                                                   initial_storage)
     summary%final_storage_fraction_living = storage_fraction_living(initial_state, &
                                                                     initial_storage)
     summary%max_storage_fraction_total = summary%final_storage_fraction_total
     summary%max_storage_fraction_living = summary%final_storage_fraction_living
     summary%min_living_pool = min(initial_state%leaf_mass, &
                                   min(initial_state%root_mass, &
                                       initial_state%sapwood_mass))
     summary%min_living_carbon = summary%initial_living_carbon
 
     ! Fail immediately if the parameter initialization produced a non-finite
     ! trait value. Additional state/result checks are performed during simulation.
     if (.not. ieee_is_finite(params%sla)) error stop "Invalid SLA."
 
   end subroutine initialize_summary
 
 
   !---------------------------------------------------------------------------
   ! Update maximum numerical errors and maximum observed daily demand/allocation
   ! using the current timestep result.
   !---------------------------------------------------------------------------
   subroutine update_summary_diagnostics(summary, result)
 
     type(ScenarioSummary), intent(inout) :: summary
     type(AllocationOutput), intent(in) :: result
 
     summary%max_abs_structural_balance_error = max(&
        summary%max_abs_structural_balance_error, &
        abs(result%structural_balance_error))
     summary%max_abs_storage_balance_error = max(&
        summary%max_abs_storage_balance_error, &
        abs(result%storage_balance_error))
     summary%max_abs_whole_plant_balance_error = max(&
        summary%max_abs_whole_plant_balance_error, &
        abs(result%whole_plant_balance_error))
     summary%max_daily_allocation_observed = max(&
        summary%max_daily_allocation_observed, result%carbon_to_allocate)
     summary%max_daily_demand_observed = max(&
        summary%max_daily_demand_observed, result%total_demand_daily)
 
   end subroutine update_summary_diagnostics
 
 
   !---------------------------------------------------------------------------
   ! Add the current timestep fluxes to scenario-wide cumulative totals.
   !
   ! These are time-integrated amounts, not final pool sizes. Positive NPP is
   ! accumulated separately from net NPP to support the diagnostic allocation
   ! ratio used in the summary output.
   !---------------------------------------------------------------------------
   subroutine accumulate_daily_fluxes(summary, result)
 
     type(ScenarioSummary), intent(inout) :: summary
     type(AllocationOutput), intent(in) :: result
 
     summary%cumulative_npp = summary%cumulative_npp + result%npp_daily
 
     if (result%npp_daily > 0.0_real64) then
        summary%cumulative_positive_npp = summary%cumulative_positive_npp + &
                                          result%npp_daily
     end if
 
     summary%cumulative_structural_allocation = &
        summary%cumulative_structural_allocation + result%carbon_to_allocate
 
     summary%cumulative_unmet_storage_deficit = &
        summary%cumulative_unmet_storage_deficit + result%unmet_storage_deficit
     summary%cumulative_leaf_starvation_loss = &
        summary%cumulative_leaf_starvation_loss + result%leaf_starvation_loss
     summary%cumulative_root_starvation_loss = &
        summary%cumulative_root_starvation_loss + result%root_starvation_loss
     summary%cumulative_sapwood_starvation_loss = &
        summary%cumulative_sapwood_starvation_loss + &
        result%sapwood_starvation_loss
     summary%cumulative_sapwood_to_heartwood_starvation = &
        summary%cumulative_sapwood_to_heartwood_starvation + &
        result%sapwood_to_heartwood
     summary%cumulative_starvation_carbon_loss = &
        summary%cumulative_starvation_carbon_loss + result%starvation_carbon_loss
     summary%cumulative_unpaid_carbon_deficit = &
        summary%cumulative_unpaid_carbon_deficit + result%unpaid_carbon_deficit
 
     summary%cumulative_leaf_turnover_loss = &
        summary%cumulative_leaf_turnover_loss + result%leaf_turnover_loss
     summary%cumulative_root_turnover_loss = &
        summary%cumulative_root_turnover_loss + result%root_turnover_loss
     summary%cumulative_sapwood_turnover_loss = &
        summary%cumulative_sapwood_turnover_loss + result%sapwood_turnover_loss
     summary%cumulative_storage_turnover_loss = &
        summary%cumulative_storage_turnover_loss + result%storage_turnover_loss
     summary%cumulative_heartwood_turnover_loss = &
        summary%cumulative_heartwood_turnover_loss + result%heartwood_turnover_loss
     summary%cumulative_sapwood_to_heartwood_turnover = &
        summary%cumulative_sapwood_to_heartwood_turnover + &
        result%sapwood_to_heartwood_turnover
     summary%cumulative_total_sapwood_to_heartwood = &
        summary%cumulative_total_sapwood_to_heartwood + &
        result%total_sapwood_to_heartwood
     summary%cumulative_turnover_carbon_loss = &
        summary%cumulative_turnover_carbon_loss + result%turnover_carbon_loss
 
   end subroutine accumulate_daily_fluxes
 
 
   !---------------------------------------------------------------------------
   ! Populate final states, net changes, final residuals, and derived ratios after
   ! the simulation ends or after a validation failure stops the daily loop.
   !---------------------------------------------------------------------------
   subroutine finalize_summary(summary, params, state, carbon_storage)
 
     type(ScenarioSummary), intent(inout) :: summary
     type(Parameters), intent(in) :: params
     type(PlantCarbonState), intent(in) :: state
     real(real64), intent(in) :: carbon_storage
 
     summary%final_leaf = state%leaf_mass
     summary%final_root = state%root_mass
     summary%final_sapwood = state%sapwood_mass
     summary%final_heartwood = state%heartwood_mass
     summary%final_height = state%height
     summary%final_storage = carbon_storage
     summary%final_living_carbon = living_carbon_of_state(state)
     summary%final_structural_carbon = structural_carbon_of_state(state)
     summary%final_total_carbon = summary%final_structural_carbon + &
                                  carbon_storage
     summary%net_living_carbon_change = summary%final_living_carbon - &
                                        summary%initial_living_carbon
     summary%net_total_carbon_change = summary%final_total_carbon - &
                                       summary%initial_total_carbon
 
     summary%final_leaf_root_residual = leaf_root_residual_for_state(params, state)
     summary%final_pipe_residual = pipe_residual_for_state(params, state)
     summary%final_storage_fraction_total = storage_fraction_total(state, &
                                                                   carbon_storage)
     summary%final_storage_fraction_living = storage_fraction_living(state, &
                                                                     carbon_storage)
 
     ! This ratio may exceed one when initial storage also finances construction;
     ! it must not be interpreted as a strict carbon-use efficiency.
     if (summary%cumulative_positive_npp > 0.0_real64) then
        summary%structural_fraction_of_positive_npp = &
           summary%cumulative_structural_allocation / &
           summary%cumulative_positive_npp
     else
        summary%structural_fraction_of_positive_npp = 0.0_real64
     end if
 
   end subroutine finalize_summary
 
 
   !---------------------------------------------------------------------------
   ! Validate one allocation step.
   !
   ! The checks distinguish numerical validity from biological disequilibrium:
   ! structural and storage balance errors must remain within tol, whereas the
   ! leaf-root and pipe-model residuals are allowed to be non-zero because the
   ! kernel performs gradual rather than instantaneous allometric adjustment.
   !---------------------------------------------------------------------------
   function step_is_valid(state, result, carbon_storage, failure_reason) &
        result(is_valid)
 
     type(PlantCarbonState), intent(in) :: state
     type(AllocationOutput), intent(in) :: result
     real(real64), intent(in) :: carbon_storage
     character(len=*), intent(out) :: failure_reason
     logical :: is_valid
 
     is_valid = .true.
     failure_reason = "OK"
 
     if (.not. all_state_values_are_finite(state)) then
        is_valid = .false.
        failure_reason = "Non-finite input state"
        return
     end if
 
     if (.not. result_values_are_finite(result)) then
        is_valid = .false.
        failure_reason = "Non-finite allocation result"
        return
     end if
 
     if (.not. ieee_is_finite(carbon_storage)) then
        is_valid = .false.
        failure_reason = "Non-finite storage"
        return
     end if
 
     if (result%leaf_mass_new < -tol .or. &
         result%root_mass_new < -tol .or. &
         result%sapwood_mass_new < -tol .or. &
         result%heartwood_mass_new < -tol) then
        is_valid = .false.
        failure_reason = "Negative structural pool"
        return
     end if
 
     if (carbon_storage < -tol) then
        is_valid = .false.
        failure_reason = "Negative storage"
        return
     end if
 
     if (result%delta_leaf < -tol .or. result%delta_root < -tol .or. &
         result%delta_sapwood < -tol) then
        is_valid = .false.
        failure_reason = "Negative structural increment"
        return
     end if
 
     if (any_negative_starvation_diagnostic(result)) then
        is_valid = .false.
        failure_reason = "Negative starvation diagnostic"
        return
     end if
 
     if (any_negative_turnover_diagnostic(result)) then
        is_valid = .false.
        failure_reason = "Negative turnover diagnostic"
        return
     end if
 
     if (abs(result%structural_balance_error) > tol) then
        is_valid = .false.
        failure_reason = "Structural balance error"
        return
     end if
 
     if (abs(result%storage_balance_error) > tol) then
        is_valid = .false.
        failure_reason = "Storage balance error"
        return
     end if
 
     if (abs(result%whole_plant_balance_error) > tol) then
        is_valid = .false.
        failure_reason = "Whole-plant balance error"
        return
     end if
 
     if (result%carbon_to_allocate - result%max_daily_allocation > tol) then
        is_valid = .false.
        failure_reason = "Daily allocation limiter violated"
        return
     end if
 
     if (result%carbon_to_allocate - result%total_demand_daily > tol) then
        is_valid = .false.
        failure_reason = "Demand limiter violated"
        return
     end if
 
   end function step_is_valid
 
 
   ! Return true if any starvation flux or deficit diagnostic is negative beyond
   ! numerical tolerance. All quantities in this group are defined as magnitudes.
   function any_negative_starvation_diagnostic(result) result(is_negative)
 
     type(AllocationOutput), intent(in) :: result
     logical :: is_negative
 
     is_negative = result%leaf_starvation_loss < -tol .or. &
                   result%root_starvation_loss < -tol .or. &
                   result%sapwood_starvation_loss < -tol .or. &
                   result%sapwood_to_heartwood < -tol .or. &
                   result%starvation_carbon_loss < -tol .or. &
                   result%unpaid_carbon_deficit < -tol
 
   end function any_negative_starvation_diagnostic
 
 
   ! Return true if any turnover or sapwood-transfer magnitude is negative beyond
   ! numerical tolerance.
   function any_negative_turnover_diagnostic(result) result(is_negative)
 
     type(AllocationOutput), intent(in) :: result
     logical :: is_negative
 
     is_negative = result%leaf_turnover_loss < -tol .or. &
                   result%root_turnover_loss < -tol .or. &
                   result%sapwood_turnover_loss < -tol .or. &
                   result%storage_turnover_loss < -tol .or. &
                   result%heartwood_turnover_loss < -tol .or. &
                   result%sapwood_to_heartwood_turnover < -tol .or. &
                   result%total_sapwood_to_heartwood < -tol .or. &
                   result%turnover_carbon_loss < -tol
 
   end function any_negative_turnover_diagnostic
 
 
   ! Verify that every structural state variable is finite before a timestep is
   ! accepted.
   function all_state_values_are_finite(state) result(is_finite)
 
     type(PlantCarbonState), intent(in) :: state
     logical :: is_finite
 
     is_finite = ieee_is_finite(state%leaf_mass) .and. &
                 ieee_is_finite(state%root_mass) .and. &
                 ieee_is_finite(state%sapwood_mass) .and. &
                 ieee_is_finite(state%heartwood_mass) .and. &
                 ieee_is_finite(state%height)
 
   end function all_state_values_are_finite
 
 
   ! Verify that all allocation outputs used by this driver are finite. This
   ! catches NaN or Inf values before they propagate into subsequent timesteps.
   function result_values_are_finite(result) result(is_finite)
 
     type(AllocationOutput), intent(in) :: result
     logical :: is_finite
 
     is_finite = ieee_is_finite(result%delta_leaf) .and. &
                 ieee_is_finite(result%delta_root) .and. &
                 ieee_is_finite(result%delta_sapwood) .and. &
                 ieee_is_finite(result%leaf_mass_new) .and. &
                 ieee_is_finite(result%root_mass_new) .and. &
                 ieee_is_finite(result%sapwood_mass_new) .and. &
                 ieee_is_finite(result%heartwood_mass_new) .and. &
                 ieee_is_finite(result%height_new) .and. &
                 ieee_is_finite(result%carbon_storage_after) .and. &
                 ieee_is_finite(result%carbon_to_allocate) .and. &
                 ieee_is_finite(result%total_demand_daily) .and. &
                 ieee_is_finite(result%leaf_starvation_loss) .and. &
                 ieee_is_finite(result%root_starvation_loss) .and. &
                 ieee_is_finite(result%sapwood_starvation_loss) .and. &
                 ieee_is_finite(result%sapwood_to_heartwood) .and. &
                 ieee_is_finite(result%starvation_carbon_loss) .and. &
                 ieee_is_finite(result%unpaid_carbon_deficit) .and. &
                 ieee_is_finite(result%leaf_turnover_loss) .and. &
                 ieee_is_finite(result%root_turnover_loss) .and. &
                 ieee_is_finite(result%sapwood_turnover_loss) .and. &
                 ieee_is_finite(result%storage_turnover_loss) .and. &
                 ieee_is_finite(result%heartwood_turnover_loss) .and. &
                 ieee_is_finite(result%sapwood_to_heartwood_turnover) .and. &
                 ieee_is_finite(result%total_sapwood_to_heartwood) .and. &
                 ieee_is_finite(result%turnover_carbon_loss) .and. &
                 ieee_is_finite(result%structural_balance_error) .and. &
                 ieee_is_finite(result%storage_balance_error) .and. &
                 ieee_is_finite(result%whole_plant_balance_error)
 
   end function result_values_are_finite
 
 
   !---------------------------------------------------------------------------
   ! Invert the stem allometry to recover height from total stem carbon.
   !
   ! This helper reproduces the allometric relation used by the kernel so that
   ! initial states and diagnostics are internally consistent.
   !---------------------------------------------------------------------------
   function height_from_total_stem_carbon(params, stem_carbon_total) &
        result(height)
 
     type(Parameters), intent(in) :: params
     real(real64), intent(in) :: stem_carbon_total
     real(real64) :: height
     real(real64) :: height_power_exponent
     real(real64) :: height_power
     real(real64) :: pi_over_four
 
     if (stem_carbon_total <= 0.0_real64 .or. &
         params%wood_density <= 0.0_real64) then
        height = 0.0_real64
        return
     end if
 
     height_power_exponent = 1.0_real64 + 2.0_real64 / params%allom3
     pi_over_four = acos(-1.0_real64) / 4.0_real64
 
     height_power = params%allom2**(2.0_real64 / params%allom3) * &
                    (stem_carbon_total / params%wood_density) / &
                    pi_over_four
 
     height = height_power**(1.0_real64 / height_power_exponent)
 
   end function height_from_total_stem_carbon
 
 
   !---------------------------------------------------------------------------
   ! Diagnose pipe-model disequilibrium for a state:
   !   residual = leaf area - LATOSA * sapwood area inferred from sapwood mass.
   ! Positive values indicate relatively more leaf area; negative values indicate
   ! relatively more sapwood area. The residual is not a carbon-balance error.
   !---------------------------------------------------------------------------
   function pipe_residual_for_state(params, state) result(pipe_residual)
 
     type(Parameters), intent(in) :: params
     type(PlantCarbonState), intent(in) :: state
     real(real64) :: pipe_residual
     real(real64) :: leaf_area
     real(real64) :: sapwood_area_from_mass
 
     leaf_area = state%leaf_mass * params%sla
 
     if (state%height > 0.0_real64 .and. &
         params%wood_density > 0.0_real64) then
        sapwood_area_from_mass = state%sapwood_mass / &
           (params%wood_density * state%height)
     else
        sapwood_area_from_mass = 0.0_real64
     end if
 
     pipe_residual = leaf_area - params%latosa * sapwood_area_from_mass
 
   end function pipe_residual_for_state
 
 
   ! Diagnose leaf-root disequilibrium. Positive values indicate relatively more
   ! leaf carbon than prescribed by the target ratio; negative values indicate
   ! relatively more fine-root carbon.
   function leaf_root_residual_for_state(params, state) result(leaf_root_residual)
 
     type(Parameters), intent(in) :: params
     type(PlantCarbonState), intent(in) :: state
     real(real64) :: leaf_root_residual
 
     leaf_root_residual = state%leaf_mass - &
                          params%leaf_to_root_ratio * state%root_mass
 
   end function leaf_root_residual_for_state
 
 
   ! Return living structural carbon used by the daily construction limiter:
   ! leaf + fine root + sapwood. Heartwood and labile storage are excluded.
   function living_carbon_of_state(state) result(living_carbon)
 
     type(PlantCarbonState), intent(in) :: state
     real(real64) :: living_carbon
 
     living_carbon = state%leaf_mass + state%root_mass + state%sapwood_mass
 
   end function living_carbon_of_state
 
 
   ! Return total structural carbon: leaf + fine root + sapwood + heartwood.
   ! Labile storage is excluded and must be added separately for total plant carbon.
   function structural_carbon_of_state(state) result(structural_carbon)
 
     type(PlantCarbonState), intent(in) :: state
     real(real64) :: structural_carbon
 
     structural_carbon = state%leaf_mass + state%root_mass + &
                         state%sapwood_mass + state%heartwood_mass
 
   end function structural_carbon_of_state
 
 
   ! Return storage as a fraction of total tracked plant carbon:
   ! storage / (structural carbon + storage). This metric is bounded between zero
   ! and one when all pools are non-negative.
   function storage_fraction_total(state, carbon_storage) result(fraction)
 
     type(PlantCarbonState), intent(in) :: state
     real(real64), intent(in) :: carbon_storage
     real(real64) :: fraction
     real(real64) :: total_carbon
 
     total_carbon = structural_carbon_of_state(state) + carbon_storage
 
     if (total_carbon > 0.0_real64) then
        fraction = carbon_storage / total_carbon
     else
        fraction = 0.0_real64
     end if
 
   end function storage_fraction_total
 
 
   ! Return storage relative to living structural carbon:
   ! storage / (leaf + fine root + sapwood). Unlike storage_fraction_total, this
   ! ratio can exceed one when storage is larger than all living structural pools.
   function storage_fraction_living(state, carbon_storage) result(fraction)
 
     type(PlantCarbonState), intent(in) :: state
     real(real64), intent(in) :: carbon_storage
     real(real64) :: fraction
     real(real64) :: living_carbon
 
     living_carbon = living_carbon_of_state(state)
 
     if (living_carbon > 0.0_real64) then
        fraction = carbon_storage / living_carbon
     else
        fraction = 0.0_real64
     end if
 
   end function storage_fraction_living
 
 
   !---------------------------------------------------------------------------
   ! Solve for the sapwood mass that makes the initial pipe-model residual zero.
   !
   ! The solver first expands an upper bound until the residual changes sign, then
   ! applies bisection. Height changes with total stem carbon, so sapwood mass
   ! cannot be obtained from a single direct linear expression here.
   !---------------------------------------------------------------------------
   function solve_sapwood_for_pipe_balance(params, leaf_mass, heartwood_mass) &
        result(sapwood_mass)
 
     type(Parameters), intent(in) :: params
     real(real64), intent(in) :: leaf_mass
     real(real64), intent(in) :: heartwood_mass
     real(real64) :: sapwood_mass
     real(real64) :: left
     real(real64) :: right
     real(real64) :: mid
     real(real64) :: f_mid
     integer :: iter
 
     ! Begin with a strictly positive lower bound and expand the upper bound until
     ! the root is bracketed by residuals with opposite signs.
     left = 1.0e-12_real64
     right = 1.0_real64

     do while (pipe_balance_function(params, leaf_mass, heartwood_mass, &
                                     right) > 0.0_real64)
        right = right * 2.0_real64
        if (right > 1.0e12_real64) then
           error stop "Could not bracket sapwood mass for pipe balance."
        end if
     end do
 
     ! Bisect the bracket until the pipe residual is negligible or the maximum
     ! number of iterations is reached.
     do iter = 1, 300
        mid = 0.5_real64 * (left + right)
        f_mid = pipe_balance_function(params, leaf_mass, heartwood_mass, mid)
 
        if (abs(f_mid) <= 1.0e-12_real64) exit
 
        if (f_mid > 0.0_real64) then
           left = mid
        else
           right = mid
        end if
     end do
 
     sapwood_mass = 0.5_real64 * (left + right)
 
   end function solve_sapwood_for_pipe_balance
 
 
   ! Evaluate the pipe-model residual for a candidate sapwood mass. This scalar
   ! function is used by the bisection solver above.
   function pipe_balance_function(params, leaf_mass, heartwood_mass, &
                                  sapwood_mass) result(residual)
 
     type(Parameters), intent(in) :: params
     real(real64), intent(in) :: leaf_mass
     real(real64), intent(in) :: heartwood_mass
     real(real64), intent(in) :: sapwood_mass
     real(real64) :: residual
     real(real64) :: height
     real(real64) :: leaf_area
     real(real64) :: sapwood_area_from_mass
 
     height = height_from_total_stem_carbon(params, &
        sapwood_mass + heartwood_mass)
     leaf_area = leaf_mass * params%sla
 
     if (height > 0.0_real64) then
        sapwood_area_from_mass = sapwood_mass / &
           (params%wood_density * height)
     else
        sapwood_area_from_mass = 0.0_real64
     end if
 
     residual = leaf_area - params%latosa * sapwood_area_from_mass
 
   end function pipe_balance_function
 
 
   ! Return true on the fixed reporting checkpoints at the end of years 1, 5,
   ! and 10. These dates assume 365-day years and should be updated if n_days or
   ! the timestep convention changes.
   function is_checkpoint_day(day) result(is_checkpoint)
 
     integer, intent(in) :: day
     logical :: is_checkpoint
 
     is_checkpoint = (day == 365) .or. &
                     (day == 365 * 5) .or. &
                     (day == 365 * 10)
 
   end function is_checkpoint_day
 
 
   ! Convert a recognized checkpoint day to its nominal simulation year. Return
   ! -1 for a day that is not one of the predefined checkpoints.
   function checkpoint_year_from_day(day) result(year)
 
     integer, intent(in) :: day
     integer :: year
 
     select case (day)
     case (365)
        year = 1
     case (365 * 5)
        year = 5
     case (365 * 10)
        year = 10
     case default
        year = -1
     end select
 
   end function checkpoint_year_from_day
 
 
   !---------------------------------------------------------------------------
   ! Select a small subset of scenarios for detailed time-series output.
   !
   ! Rows are written every 30 days and on the final day only when all selection
   ! criteria below are satisfied. Note that the current active initial storage
   ! value (0.5) is not selected, so the daily CSV will contain only its header
   ! unless selected_storage or the active scenario values are changed.
   !---------------------------------------------------------------------------
   function should_write_daily_trace(summary, day) result(write_trace)
 
     type(ScenarioSummary), intent(in) :: summary
     integer, intent(in) :: day
     logical :: write_trace
     logical :: selected_npp
     logical :: selected_storage
 
     ! Explicitly select representative NPP and storage levels. Values not listed
     ! here are omitted from the detailed daily output even if they are simulated.
     selected_npp = abs(summary%npp_rate + 0.5_real64) < tiny_positive .or. &
                    abs(summary%npp_rate - 0.5_real64) < tiny_positive .or. &
                    abs(summary%npp_rate - 3.5_real64) < tiny_positive .or. &
                    abs(summary%npp_rate - 8.0_real64) < tiny_positive
 
     selected_storage = abs(summary%initial_storage - 0.0_real64) < &
                        tiny_positive .or. &
                        abs(summary%initial_storage - 5.0_real64) < &
                        tiny_positive
 
 
     write_trace = summary%trait_case == 1 .and. selected_npp .and. &
                   selected_storage .and. &
                   abs(summary%allometric_adjustment_days - 365.0_real64) &
                   < tiny_positive .and. &
                   abs(summary%max_allocation_fraction - 0.005_real64) &
                   < tiny_positive .and. &
                   (mod(day, 30) == 0 .or. day == n_days)
 
   end function should_write_daily_trace
 
 
   ! Write column names for the one-row-per-scenario summary file. Column order
   ! must remain synchronized with write_summary_row.
   subroutine write_summary_header(summary_unit)
 
     integer, intent(in) :: summary_unit
 
     write(summary_unit,'(a)') &
       "scenario_id,status,fail_day,failure_reason," // &
       "trait_case,background_mode,npp_rate,initial_storage," // &
       "allometric_adjustment_days,max_allocation_fraction," // &
       "initial_leaf,initial_root,initial_sapwood,initial_heartwood," // &
       "initial_height,initial_living_carbon,initial_structural_carbon," // &
       "initial_total_carbon,final_leaf,final_root,final_sapwood," // &
       "final_heartwood,final_height,final_storage,final_living_carbon," // &
       "final_structural_carbon,final_total_carbon," // &
       "net_living_carbon_change,net_total_carbon_change," // &
       "cumulative_npp,cumulative_positive_npp," // &
       "cumulative_structural_allocation," // &
       "cumulative_unmet_storage_deficit,cumulative_leaf_starvation_loss," // &
       "cumulative_root_starvation_loss,cumulative_sapwood_starvation_loss," // &
       "cumulative_sapwood_to_heartwood_starvation," // &
       "cumulative_starvation_carbon_loss," // &
       "cumulative_unpaid_carbon_deficit," // &
       "cumulative_leaf_turnover_loss,cumulative_root_turnover_loss," // &
       "cumulative_sapwood_turnover_loss,cumulative_storage_turnover_loss," // &
       "cumulative_heartwood_turnover_loss," // &
       "cumulative_sapwood_to_heartwood_turnover," // &
       "cumulative_total_sapwood_to_heartwood," // &
       "cumulative_turnover_carbon_loss," // &
       "final_leaf_root_residual,final_pipe_residual," // &
       "final_storage_fraction_total,final_storage_fraction_living," // &
       "max_storage_fraction_total,max_storage_fraction_living," // &
       "structural_fraction_of_positive_npp," // &
       "max_abs_structural_balance_error,max_abs_storage_balance_error," // &
       "max_abs_whole_plant_balance_error,min_living_pool," // &
       "min_living_carbon,max_daily_allocation_observed," // &
       "max_daily_demand_observed,days_with_allocation," // &
       "days_with_unmet_storage_deficit,days_with_turnover," // &
       "days_with_storage_above_10pct_living," // &
       "days_with_storage_above_50pct_living," // &
       "days_with_storage_above_100pct_living"
 
   end subroutine write_summary_header
 
 
   ! Write one complete scenario summary in the same order as the header above.
   subroutine write_summary_row(summary_unit, summary)
 
     integer, intent(in) :: summary_unit
     type(ScenarioSummary), intent(in) :: summary
     character(len=8) :: status
 
     if (summary%passed) then
        status = "PASS"
     else
        status = "FAIL"
     end if
 
     ! The unlimited-repeat format writes mixed numeric and character fields
     ! without requiring a manually maintained format descriptor for every column.
     write(summary_unit,'(*(g0))') &
       summary%scenario_id, ",", trim(status), ",", summary%fail_day, ",", &
       trim(summary%failure_reason), ",", summary%trait_case, ",", &
       summary%background_mode, ",", &
       summary%npp_rate, ",", summary%initial_storage, ",", &
       summary%allometric_adjustment_days, ",", &
       summary%max_allocation_fraction, ",", summary%initial_leaf, ",", &
       summary%initial_root, ",", summary%initial_sapwood, ",", &
       summary%initial_heartwood, ",", summary%initial_height, ",", &
       summary%initial_living_carbon, ",", &
       summary%initial_structural_carbon, ",", &
       summary%initial_total_carbon, ",", summary%final_leaf, ",", &
       summary%final_root, ",", summary%final_sapwood, ",", &
       summary%final_heartwood, ",", summary%final_height, ",", &
       summary%final_storage, ",", summary%final_living_carbon, ",", &
       summary%final_structural_carbon, ",", &
       summary%final_total_carbon, ",", &
       summary%net_living_carbon_change, ",", &
       summary%net_total_carbon_change, ",", summary%cumulative_npp, ",", &
       summary%cumulative_positive_npp, ",", &
       summary%cumulative_structural_allocation, ",", &
       summary%cumulative_unmet_storage_deficit, ",", &
       summary%cumulative_leaf_starvation_loss, ",", &
       summary%cumulative_root_starvation_loss, ",", &
       summary%cumulative_sapwood_starvation_loss, ",", &
       summary%cumulative_sapwood_to_heartwood_starvation, ",", &
       summary%cumulative_starvation_carbon_loss, ",", &
       summary%cumulative_unpaid_carbon_deficit, ",", &
       summary%cumulative_leaf_turnover_loss, ",", &
       summary%cumulative_root_turnover_loss, ",", &
       summary%cumulative_sapwood_turnover_loss, ",", &
       summary%cumulative_storage_turnover_loss, ",", &
       summary%cumulative_heartwood_turnover_loss, ",", &
       summary%cumulative_sapwood_to_heartwood_turnover, ",", &
       summary%cumulative_total_sapwood_to_heartwood, ",", &
       summary%cumulative_turnover_carbon_loss, ",", &
       summary%final_leaf_root_residual, ",", summary%final_pipe_residual, ",", &
       summary%final_storage_fraction_total, ",", &
       summary%final_storage_fraction_living, ",", &
       summary%max_storage_fraction_total, ",", &
       summary%max_storage_fraction_living, ",", &
       summary%structural_fraction_of_positive_npp, ",", &
       summary%max_abs_structural_balance_error, ",", &
       summary%max_abs_storage_balance_error, ",", &
       summary%max_abs_whole_plant_balance_error, ",", &
       summary%min_living_pool, ",", summary%min_living_carbon, ",", &
       summary%max_daily_allocation_observed, ",", &
       summary%max_daily_demand_observed, ",", &
       summary%days_with_allocation, ",", &
       summary%days_with_unmet_storage_deficit, ",", &
       summary%days_with_turnover, ",", &
       summary%days_with_storage_above_10pct_living, ",", &
       summary%days_with_storage_above_50pct_living, ",", &
       summary%days_with_storage_above_100pct_living
 
   end subroutine write_summary_row
 
 
   ! Write column names for the checkpoint file. Column order must remain
   ! synchronized with write_checkpoint_row.
   subroutine write_checkpoint_header(checkpoint_unit)
 
     integer, intent(in) :: checkpoint_unit
 
     write(checkpoint_unit,'(a)') &
       "scenario_id,checkpoint_year,checkpoint_day,status,fail_day," // &
       "failure_reason,trait_case,background_mode,npp_rate," // &
       "initial_storage,allometric_adjustment_days,max_allocation_fraction," // &
       "leaf,root,sapwood,heartwood,height,storage,living_carbon," // &
       "structural_carbon,total_carbon,cumulative_npp," // &
       "cumulative_positive_npp,cumulative_structural_allocation," // &
       "cumulative_unmet_storage_deficit,cumulative_leaf_starvation_loss," // &
       "cumulative_root_starvation_loss,cumulative_sapwood_starvation_loss," // &
       "cumulative_sapwood_to_heartwood_starvation," // &
       "cumulative_starvation_carbon_loss," // &
       "cumulative_unpaid_carbon_deficit," // &
       "cumulative_leaf_turnover_loss,cumulative_root_turnover_loss," // &
       "cumulative_sapwood_turnover_loss,cumulative_storage_turnover_loss," // &
       "cumulative_heartwood_turnover_loss," // &
       "cumulative_sapwood_to_heartwood_turnover," // &
       "cumulative_total_sapwood_to_heartwood," // &
       "cumulative_turnover_carbon_loss,leaf_root_residual," // &
       "pipe_residual,storage_fraction_total,storage_fraction_living," // &
       "structural_fraction_of_positive_npp," // &
       "max_abs_structural_balance_error,max_abs_storage_balance_error," // &
       "max_abs_whole_plant_balance_error,max_storage_fraction_total," // &
       "max_storage_fraction_living,min_living_pool,min_living_carbon," // &
       "days_with_allocation,days_with_unmet_storage_deficit," // &
       "days_with_turnover,days_with_storage_above_10pct_living," // &
       "days_with_storage_above_50pct_living," // &
       "days_with_storage_above_100pct_living"
 
   end subroutine write_checkpoint_header
 
 
   ! Write state, cumulative fluxes, and diagnostics at a predefined checkpoint.
   subroutine write_checkpoint_row(checkpoint_unit, summary, params, state, &
                                   carbon_storage, day)
 
     integer, intent(in) :: checkpoint_unit
     type(ScenarioSummary), intent(in) :: summary
     type(Parameters), intent(in) :: params
     type(PlantCarbonState), intent(in) :: state
     real(real64), intent(in) :: carbon_storage
     integer, intent(in) :: day
 
     character(len=8) :: status
     integer :: checkpoint_year
     real(real64) :: current_leaf_root_residual
     real(real64) :: current_pipe_residual
     real(real64) :: current_storage_fraction_total
     real(real64) :: current_storage_fraction_living
     real(real64) :: current_living_carbon
     real(real64) :: current_structural_carbon
     real(real64) :: current_total_carbon
     real(real64) :: structural_fraction_at_checkpoint
 
     if (summary%passed) then
        status = "PASS"
     else
        status = "FAIL"
     end if
 
     ! Recalculate state-dependent diagnostics at the checkpoint rather than
     ! relying only on the final values stored in ScenarioSummary.
     checkpoint_year = checkpoint_year_from_day(day)
     current_leaf_root_residual = leaf_root_residual_for_state(params, state)
     current_pipe_residual = pipe_residual_for_state(params, state)
     current_storage_fraction_total = storage_fraction_total(state, carbon_storage)
     current_storage_fraction_living = storage_fraction_living(state, carbon_storage)
     current_living_carbon = living_carbon_of_state(state)
     current_structural_carbon = structural_carbon_of_state(state)
     current_total_carbon = current_structural_carbon + carbon_storage
 
     ! As in the final summary, this is a diagnostic index and may exceed one
     ! when construction is partly financed by initial storage.
     if (summary%cumulative_positive_npp > 0.0_real64) then
        structural_fraction_at_checkpoint = &
           summary%cumulative_structural_allocation / &
           summary%cumulative_positive_npp
     else
        structural_fraction_at_checkpoint = 0.0_real64
     end if
 
     write(checkpoint_unit,'(*(g0))') &
       summary%scenario_id, ",", checkpoint_year, ",", day, ",", &
       trim(status), ",", summary%fail_day, ",", &
       trim(summary%failure_reason), ",", summary%trait_case, ",", &
       summary%background_mode, ",", &
       summary%npp_rate, ",", summary%initial_storage, ",", &
       summary%allometric_adjustment_days, ",", &
       summary%max_allocation_fraction, ",", state%leaf_mass, ",", &
       state%root_mass, ",", state%sapwood_mass, ",", &
       state%heartwood_mass, ",", state%height, ",", carbon_storage, ",", &
       current_living_carbon, ",", current_structural_carbon, ",", &
       current_total_carbon, ",", summary%cumulative_npp, ",", &
       summary%cumulative_positive_npp, ",", &
       summary%cumulative_structural_allocation, ",", &
       summary%cumulative_unmet_storage_deficit, ",", &
       summary%cumulative_leaf_starvation_loss, ",", &
       summary%cumulative_root_starvation_loss, ",", &
       summary%cumulative_sapwood_starvation_loss, ",", &
       summary%cumulative_sapwood_to_heartwood_starvation, ",", &
       summary%cumulative_starvation_carbon_loss, ",", &
       summary%cumulative_unpaid_carbon_deficit, ",", &
       summary%cumulative_leaf_turnover_loss, ",", &
       summary%cumulative_root_turnover_loss, ",", &
       summary%cumulative_sapwood_turnover_loss, ",", &
       summary%cumulative_storage_turnover_loss, ",", &
       summary%cumulative_heartwood_turnover_loss, ",", &
       summary%cumulative_sapwood_to_heartwood_turnover, ",", &
       summary%cumulative_total_sapwood_to_heartwood, ",", &
       summary%cumulative_turnover_carbon_loss, ",", &
       current_leaf_root_residual, ",", current_pipe_residual, ",", &
       current_storage_fraction_total, ",", &
       current_storage_fraction_living, ",", &
       structural_fraction_at_checkpoint, ",", &
       summary%max_abs_structural_balance_error, ",", &
       summary%max_abs_storage_balance_error, ",", &
       summary%max_abs_whole_plant_balance_error, ",", &
       summary%max_storage_fraction_total, ",", &
       summary%max_storage_fraction_living, ",", &
       summary%min_living_pool, ",", summary%min_living_carbon, ",", &
       summary%days_with_allocation, ",", &
       summary%days_with_unmet_storage_deficit, ",", &
       summary%days_with_turnover, ",", &
       summary%days_with_storage_above_10pct_living, ",", &
       summary%days_with_storage_above_50pct_living, ",", &
       summary%days_with_storage_above_100pct_living
 
   end subroutine write_checkpoint_row
 
 
   ! Write column names for the selected detailed time-series file. Column order
   ! must remain synchronized with write_daily_row.
   subroutine write_daily_header(daily_unit)
 
     integer, intent(in) :: daily_unit
 
     write(daily_unit,'(a)') &
       "scenario_id,day,year,trait_case,background_mode," // &
       "npp_rate,initial_storage,allometric_adjustment_days," // &
       "max_allocation_fraction,leaf,root,sapwood,heartwood," // &
       "height,storage,living_carbon,structural_carbon,total_carbon," // &
       "npp_daily,storage_after_npp_unclamped,unmet_storage_deficit," // &
       "carbon_to_allocate,total_demand_daily,leaf_demand_daily," // &
       "root_demand_daily,sapwood_demand_daily,delta_leaf,delta_root," // &
       "delta_sapwood,leaf_starvation_loss,root_starvation_loss," // &
       "sapwood_starvation_loss,sapwood_to_heartwood_starvation," // &
       "starvation_carbon_loss,unpaid_carbon_deficit," // &
       "leaf_turnover_loss,root_turnover_loss,sapwood_turnover_loss," // &
       "storage_turnover_loss,heartwood_turnover_loss," // &
       "sapwood_to_heartwood_turnover,total_sapwood_to_heartwood," // &
       "turnover_carbon_loss,storage_fraction_total," // &
       "storage_fraction_living,leaf_root_residual,pipe_residual," // &
       "structural_balance_error,storage_balance_error," // &
       "whole_plant_balance_error"
 
   end subroutine write_daily_header
 
 
   ! Write one selected post-step plant state together with the fluxes and
   ! diagnostics generated during that same timestep.
   subroutine write_daily_row(daily_unit, summary, params, state, result, &
                              carbon_storage, day)
 
     integer, intent(in) :: daily_unit
     type(ScenarioSummary), intent(in) :: summary
     type(Parameters), intent(in) :: params
     type(PlantCarbonState), intent(in) :: state
     type(AllocationOutput), intent(in) :: result
     real(real64), intent(in) :: carbon_storage
     integer, intent(in) :: day
 
     ! state and carbon_storage are post-step values, while result contains the
     ! fluxes and balance diagnostics that produced that post-step state.
     write(daily_unit,'(*(g0))') &
       summary%scenario_id, ",", day, ",", &
       real(day, real64) / 365.0_real64, ",", &
       summary%trait_case, ",",  &
       summary%background_mode, ",", summary%npp_rate, ",", &
       summary%initial_storage, ",", &
       summary%allometric_adjustment_days, ",", &
       summary%max_allocation_fraction, ",", state%leaf_mass, ",", &
       state%root_mass, ",", state%sapwood_mass, ",", &
       state%heartwood_mass, ",", state%height, ",", carbon_storage, ",", &
       living_carbon_of_state(state), ",", &
       structural_carbon_of_state(state), ",", &
       structural_carbon_of_state(state) + carbon_storage, ",", &
       result%npp_daily, ",", result%storage_after_npp_unclamped, ",", &
       result%unmet_storage_deficit, ",", result%carbon_to_allocate, ",", &
       result%total_demand_daily, ",", result%leaf_demand_daily, ",", &
       result%root_demand_daily, ",", result%sapwood_demand_daily, ",", &
       result%delta_leaf, ",", result%delta_root, ",", &
       result%delta_sapwood, ",", result%leaf_starvation_loss, ",", &
       result%root_starvation_loss, ",", &
       result%sapwood_starvation_loss, ",", result%sapwood_to_heartwood, ",", &
       result%starvation_carbon_loss, ",", &
       result%unpaid_carbon_deficit, ",", result%leaf_turnover_loss, ",", &
       result%root_turnover_loss, ",", result%sapwood_turnover_loss, ",", &
       result%storage_turnover_loss, ",", &
       result%heartwood_turnover_loss, ",", &
       result%sapwood_to_heartwood_turnover, ",", &
       result%total_sapwood_to_heartwood, ",", &
       result%turnover_carbon_loss, ",", &
       storage_fraction_total(state, carbon_storage), ",", &
       storage_fraction_living(state, carbon_storage), ",", &
       leaf_root_residual_for_state(params, state), ",", &
       pipe_residual_for_state(params, state), ",", &
       result%structural_balance_error, ",", &
       result%storage_balance_error, ",", &
       result%whole_plant_balance_error
 
   end subroutine write_daily_row
 
 end program run_carbon_allocation_offline
 