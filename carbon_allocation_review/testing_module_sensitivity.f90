program test_storage_allocation_sensitivity

  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use carbon_allocation_offline_kernel

  implicit none

  integer, parameter :: n_days = 365 * 10
  real(real64), parameter :: tol = 1.0e-10_real64
  real(real64), parameter :: tiny_positive = 1.0e-12_real64

  integer, parameter :: n_npp = 5
  integer, parameter :: n_storage = 3
  integer, parameter :: n_adjustment = 3
  integer, parameter :: n_max_fraction = 3
  integer, parameter :: n_background = 2
  integer, parameter :: n_trait_case = 4
  integer, parameter :: n_state_case = 5

  type :: ScenarioSummary

     ! Scenario identifiers.
     integer :: scenario_id = 0
     logical :: passed = .true.
     integer :: fail_day = 0
     character(len=96) :: failure_reason = "OK"

     ! Scenario factors.
     integer :: trait_case = 0
     integer :: state_case = 0
     integer :: background_mode = 0
     real(real64) :: npp_rate = 0.0_real64
     real(real64) :: initial_storage = 0.0_real64
     real(real64) :: allometric_adjustment_days = 0.0_real64
     real(real64) :: max_allocation_fraction = 0.0_real64

     ! Initial state.
     real(real64) :: initial_leaf = 0.0_real64
     real(real64) :: initial_root = 0.0_real64
     real(real64) :: initial_sapwood = 0.0_real64
     real(real64) :: initial_heartwood = 0.0_real64
     real(real64) :: initial_height = 0.0_real64

     ! Final state.
     real(real64) :: final_leaf = 0.0_real64
     real(real64) :: final_root = 0.0_real64
     real(real64) :: final_sapwood = 0.0_real64
     real(real64) :: final_heartwood = 0.0_real64
     real(real64) :: final_height = 0.0_real64
     real(real64) :: final_storage = 0.0_real64

     ! Integrated carbon fluxes.
     real(real64) :: cumulative_npp = 0.0_real64
     real(real64) :: cumulative_structural_allocation = 0.0_real64
     real(real64) :: cumulative_unmet_storage_deficit = 0.0_real64

     ! Diagnostic metrics.
     real(real64) :: final_leaf_root_residual = 0.0_real64
     real(real64) :: final_pipe_residual = 0.0_real64
     real(real64) :: final_storage_fraction = 0.0_real64
     real(real64) :: structural_fraction_of_positive_npp = 0.0_real64
     real(real64) :: max_abs_structural_balance_error = 0.0_real64
     real(real64) :: max_abs_storage_balance_error = 0.0_real64
     real(real64) :: max_abs_whole_plant_balance_error = 0.0_real64
     real(real64) :: max_storage_fraction = 0.0_real64
     real(real64) :: min_living_pool = huge(1.0_real64)
     real(real64) :: max_daily_allocation_observed = 0.0_real64
     real(real64) :: max_daily_demand_observed = 0.0_real64
     integer :: days_with_allocation = 0
     integer :: days_with_unmet_storage_deficit = 0

  end type ScenarioSummary

  real(real64), dimension(n_npp) :: npp_values
  real(real64), dimension(n_storage) :: storage_values
  real(real64), dimension(n_adjustment) :: adjustment_values
  real(real64), dimension(n_max_fraction) :: max_fraction_values

  integer :: i_npp
  integer :: i_storage
  integer :: i_adjustment
  integer :: i_max_fraction
  integer :: i_background
  integer :: i_trait_case
  integer :: i_state_case
  integer :: scenario_id
  integer :: n_pass
  integer :: n_fail
  integer :: csv_unit

  type(Parameters) :: params
  type(StorageAllocationControls) :: controls
  type(PlantCarbonState) :: initial_state
  type(ScenarioSummary) :: summary

  npp_values = [ -0.5_real64, 0.0_real64, 0.5_real64, 3.5_real64, 8.0_real64 ]
  storage_values = [ 0.0_real64, 0.5_real64, 5.0_real64 ]
  adjustment_values = [ 30.0_real64, 365.0_real64, 730.0_real64 ]
  max_fraction_values = [ 0.001_real64, 0.005_real64, 0.02_real64 ]

  scenario_id = 0
  n_pass = 0
  n_fail = 0

  open(newunit=csv_unit, file="storage_allocation_sensitivity_summary.csv", &
       status="replace", action="write")

  call write_csv_header(csv_unit)

  do i_trait_case = 1, n_trait_case
     do i_state_case = 1, n_state_case
        do i_background = 1, n_background
           do i_max_fraction = 1, n_max_fraction
              do i_adjustment = 1, n_adjustment
                 do i_storage = 1, n_storage
                    do i_npp = 1, n_npp

                       scenario_id = scenario_id + 1

                       call initialize_parameters(params, i_trait_case)
                       call initialize_controls(controls, adjustment_values(i_adjustment), &
                                                max_fraction_values(i_max_fraction), &
                                                i_background)
                       call initialize_state(params, i_state_case, initial_state)

                       call run_scenario(scenario_id, params, controls, initial_state, &
                                         npp_values(i_npp), storage_values(i_storage), &
                                         i_trait_case, i_state_case, i_background, &
                                         summary)

                       call write_csv_row(csv_unit, summary)

                       if (summary%passed) then
                          n_pass = n_pass + 1
                       else
                          n_fail = n_fail + 1
                       end if

                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  close(csv_unit)

  write(*,'(a,i0)') "Sensitivity scenarios completed: ", scenario_id
  write(*,'(a,i0)') "PASS: ", n_pass
  write(*,'(a,i0)') "FAIL: ", n_fail
  write(*,'(a)') "Output written to storage_allocation_sensitivity_summary.csv"

  if (n_fail > 0) then
     error stop "At least one sensitivity scenario failed. Inspect the CSV file."
  end if

contains

  subroutine initialize_parameters(params, trait_case)

    type(Parameters), intent(out) :: params
    integer, intent(in) :: trait_case

    ! Base parameter set. The sensitivity cases below perturb one trait at a time
    ! while keeping the same allometric coefficients.
    params%sla                = 12.0_real64
    params%latosa             = 8000.0_real64
    params%wood_density       = 250.0_real64
    params%leaf_to_root_ratio = 1.0_real64
    params%allom2             = 40.0_real64
    params%allom3             = 0.5_real64

    select case (trait_case)
    case (1)
       ! Base trait combination.
    case (2)
       ! Low-SLA strategy: lower leaf area per unit leaf carbon.
       params%sla = 6.0_real64
    case (3)
       ! High-SLA strategy: higher leaf area per unit leaf carbon.
       params%sla = 24.0_real64
    case (4)
       ! Denser wood strategy.
       params%wood_density = 500.0_real64
    case default
       error stop "Unknown trait case."
    end select

  end subroutine initialize_parameters


  subroutine initialize_controls(controls, adjustment_days, max_fraction, background_mode)

    type(StorageAllocationControls), intent(out) :: controls
    real(real64), intent(in) :: adjustment_days
    real(real64), intent(in) :: max_fraction
    integer, intent(in) :: background_mode

    controls%dt_years = 1.0_real64 / 365.0_real64
    controls%allometric_adjustment_days = adjustment_days
    controls%max_allocation_fraction = max_fraction

    if (background_mode == 1) then
       ! Background balanced-growth demand enabled.
       controls%leaf_background_timescale_years = 3.0_real64
       controls%root_background_timescale_years = 3.0_real64
       controls%sapwood_background_timescale_years = 15.0_real64
    else if (background_mode == 2) then
       ! Background balanced-growth demand disabled. Allocation then depends only
       ! on explicit allometric correction deficits.
       controls%leaf_background_timescale_years = 0.0_real64
       controls%root_background_timescale_years = 0.0_real64
       controls%sapwood_background_timescale_years = 0.0_real64
    else
       error stop "Unknown background mode."
    end if

  end subroutine initialize_controls


  subroutine initialize_state(params, state_case, state)

    type(Parameters), intent(in) :: params
    integer, intent(in) :: state_case
    type(PlantCarbonState), intent(out) :: state

    real(real64) :: leaf_mass
    real(real64) :: root_mass
    real(real64) :: sapwood_mass
    real(real64) :: heartwood_mass

    select case (state_case)
    case (1)
       ! Original imbalanced test state: leaf mass is higher than root mass,
       ! but leaf area is too low relative to sapwood area.
       leaf_mass = 1.0_real64
       root_mass = 0.8_real64
       sapwood_mass = 10.0_real64
       heartwood_mass = 20.0_real64

    case (2)
       ! Approximately equilibrated state for both leaf-root balance and the
       ! pipe-model relationship.
       leaf_mass = 1.0_real64
       root_mass = leaf_mass / params%leaf_to_root_ratio
       heartwood_mass = 20.0_real64
       sapwood_mass = solve_sapwood_for_pipe_balance(params, leaf_mass, heartwood_mass)

    case (3)
       ! Leaf-rich state.
       leaf_mass = 3.0_real64
       root_mass = 0.8_real64
       sapwood_mass = 10.0_real64
       heartwood_mass = 20.0_real64

    case (4)
       ! Root-rich state.
       leaf_mass = 0.8_real64
       root_mass = 3.0_real64
       sapwood_mass = 10.0_real64
       heartwood_mass = 20.0_real64

    case (5)
       ! Sapwood-rich state: this stresses the pipe-model correction term.
       leaf_mass = 1.0_real64
       root_mass = 1.0_real64
       sapwood_mass = 30.0_real64
       heartwood_mass = 20.0_real64

    case default
       error stop "Unknown state case."
    end select

    state = build_state(params, leaf_mass, root_mass, sapwood_mass, heartwood_mass)

  end subroutine initialize_state


  function build_state(params, leaf_mass, root_mass, sapwood_mass, heartwood_mass) result(state)

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
    state%height = height_from_total_stem_carbon(params, sapwood_mass + heartwood_mass)

  end function build_state


  subroutine run_scenario(scenario_id, params, controls, initial_state, npp_rate, &
                          initial_storage, trait_case, state_case, background_mode, &
                          summary)

    integer, intent(in) :: scenario_id
    type(Parameters), intent(in) :: params
    type(StorageAllocationControls), intent(in) :: controls
    type(PlantCarbonState), intent(in) :: initial_state
    real(real64), intent(in) :: npp_rate
    real(real64), intent(in) :: initial_storage
    integer, intent(in) :: trait_case
    integer, intent(in) :: state_case
    integer, intent(in) :: background_mode
    type(ScenarioSummary), intent(out) :: summary

    type(PlantCarbonState) :: state
    type(AllocationOutput) :: result
    real(real64) :: carbon_storage
    real(real64) :: positive_npp_integral
    real(real64) :: current_storage_fraction
    integer :: day

    call initialize_summary(summary, scenario_id, params, controls, initial_state, &
                            npp_rate, initial_storage, trait_case, state_case, &
                            background_mode)

    state = initial_state
    carbon_storage = initial_storage
    positive_npp_integral = 0.0_real64

    do day = 1, n_days

       call allocate_gradual_with_storage(state, params, controls, npp_rate, &
                                          carbon_storage, result)

       call update_summary_diagnostics(summary, state, result, carbon_storage, &
                                       npp_rate, positive_npp_integral)

       if (.not. step_is_valid(state, result, carbon_storage, summary%failure_reason)) then
          summary%passed = .false.
          summary%fail_day = day
          exit
       end if

       state%leaf_mass = result%leaf_mass_new
       state%root_mass = result%root_mass_new
       state%sapwood_mass = result%sapwood_mass_new
       state%heartwood_mass = result%heartwood_mass_new
       state%height = result%height_new

       summary%cumulative_npp = summary%cumulative_npp + result%npp_daily
       summary%cumulative_structural_allocation = summary%cumulative_structural_allocation + &
                                                  result%carbon_to_allocate
       summary%cumulative_unmet_storage_deficit = summary%cumulative_unmet_storage_deficit + &
                                                  result%unmet_storage_deficit

       if (result%carbon_to_allocate > tiny_positive) then
          summary%days_with_allocation = summary%days_with_allocation + 1
       end if

       if (result%unmet_storage_deficit > tiny_positive) then
          summary%days_with_unmet_storage_deficit = &
             summary%days_with_unmet_storage_deficit + 1
       end if

       if (npp_rate > 0.0_real64) then
          positive_npp_integral = positive_npp_integral + result%npp_daily
       end if

       current_storage_fraction = storage_fraction(state, carbon_storage)
       summary%max_storage_fraction = max(summary%max_storage_fraction, current_storage_fraction)

    end do

    summary%final_leaf = state%leaf_mass
    summary%final_root = state%root_mass
    summary%final_sapwood = state%sapwood_mass
    summary%final_heartwood = state%heartwood_mass
    summary%final_height = state%height
    summary%final_storage = carbon_storage
    summary%final_leaf_root_residual = leaf_root_residual_for_state(params, state)
    summary%final_pipe_residual = pipe_residual_for_state(params, state)
    summary%final_storage_fraction = storage_fraction(state, carbon_storage)

    if (positive_npp_integral > 0.0_real64) then
       summary%structural_fraction_of_positive_npp = &
          summary%cumulative_structural_allocation / positive_npp_integral
    else
       summary%structural_fraction_of_positive_npp = 0.0_real64
    end if

  end subroutine run_scenario


  subroutine initialize_summary(summary, scenario_id, params, controls, initial_state, &
                                npp_rate, initial_storage, trait_case, state_case, &
                                background_mode)

    type(ScenarioSummary), intent(out) :: summary
    integer, intent(in) :: scenario_id
    type(Parameters), intent(in) :: params
    type(StorageAllocationControls), intent(in) :: controls
    type(PlantCarbonState), intent(in) :: initial_state
    real(real64), intent(in) :: npp_rate
    real(real64), intent(in) :: initial_storage
    integer, intent(in) :: trait_case
    integer, intent(in) :: state_case
    integer, intent(in) :: background_mode

    summary%scenario_id = scenario_id
    summary%passed = .true.
    summary%fail_day = 0
    summary%failure_reason = "OK"

    summary%trait_case = trait_case
    summary%state_case = state_case
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

    summary%final_leaf = initial_state%leaf_mass
    summary%final_root = initial_state%root_mass
    summary%final_sapwood = initial_state%sapwood_mass
    summary%final_heartwood = initial_state%heartwood_mass
    summary%final_height = initial_state%height
    summary%final_storage = initial_storage

    summary%cumulative_npp = 0.0_real64
    summary%cumulative_structural_allocation = 0.0_real64
    summary%cumulative_unmet_storage_deficit = 0.0_real64

    summary%final_leaf_root_residual = leaf_root_residual_for_state(params, initial_state)
    summary%final_pipe_residual = pipe_residual_for_state(params, initial_state)
    summary%final_storage_fraction = storage_fraction(initial_state, initial_storage)
    summary%structural_fraction_of_positive_npp = 0.0_real64
    summary%max_abs_structural_balance_error = 0.0_real64
    summary%max_abs_storage_balance_error = 0.0_real64
    summary%max_abs_whole_plant_balance_error = 0.0_real64
    summary%max_storage_fraction = storage_fraction(initial_state, initial_storage)
    summary%min_living_pool = min(initial_state%leaf_mass, &
                                  min(initial_state%root_mass, initial_state%sapwood_mass))
    summary%max_daily_allocation_observed = 0.0_real64
    summary%max_daily_demand_observed = 0.0_real64
    summary%days_with_allocation = 0
    summary%days_with_unmet_storage_deficit = 0

    ! Prevent compiler warnings when the parameter object is not directly used
    ! in this routine beyond the diagnostic helpers above.
    if (.not. ieee_is_finite(params%sla)) error stop "Invalid SLA."

  end subroutine initialize_summary


  subroutine update_summary_diagnostics(summary, state, result, carbon_storage, &
                                        npp_rate, positive_npp_integral)

    type(ScenarioSummary), intent(inout) :: summary
    type(PlantCarbonState), intent(in) :: state
    type(AllocationOutput), intent(in) :: result
    real(real64), intent(in) :: carbon_storage
    real(real64), intent(in) :: npp_rate
    real(real64), intent(in) :: positive_npp_integral

    ! Track maximum numerical accounting errors over the full simulation.
    summary%max_abs_structural_balance_error = max(summary%max_abs_structural_balance_error, &
                                                   abs(result%structural_balance_error))
    summary%max_abs_storage_balance_error = max(summary%max_abs_storage_balance_error, &
                                                abs(result%storage_balance_error))
    summary%max_abs_whole_plant_balance_error = max(summary%max_abs_whole_plant_balance_error, &
                                                    abs(result%whole_plant_balance_error))

    ! Track physically relevant magnitudes.
    summary%min_living_pool = min(summary%min_living_pool, result%leaf_mass_new)
    summary%min_living_pool = min(summary%min_living_pool, result%root_mass_new)
    summary%min_living_pool = min(summary%min_living_pool, result%sapwood_mass_new)
    summary%max_daily_allocation_observed = max(summary%max_daily_allocation_observed, &
                                                result%carbon_to_allocate)
    summary%max_daily_demand_observed = max(summary%max_daily_demand_observed, &
                                            result%total_demand_daily)

    ! These variables are intentionally referenced so that strict compilers do
    ! not warn about unused dummy arguments. They are useful context if this
    ! diagnostic routine is extended later.
    if (.not. ieee_is_finite(carbon_storage)) summary%failure_reason = "Non-finite storage"
    if (.not. ieee_is_finite(npp_rate)) summary%failure_reason = "Non-finite NPP"
    if (.not. ieee_is_finite(positive_npp_integral)) summary%failure_reason = "Non-finite positive NPP integral"
    if (.not. ieee_is_finite(state%leaf_mass)) summary%failure_reason = "Non-finite state"

  end subroutine update_summary_diagnostics


  function step_is_valid(state, result, carbon_storage, failure_reason) result(is_valid)

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

    if (result%leaf_mass_new < -tol .or. result%root_mass_new < -tol .or. &
        result%sapwood_mass_new < -tol .or. result%heartwood_mass_new < -tol) then
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


  function all_state_values_are_finite(state) result(is_finite)

    type(PlantCarbonState), intent(in) :: state
    logical :: is_finite

    is_finite = ieee_is_finite(state%leaf_mass) .and. &
                ieee_is_finite(state%root_mass) .and. &
                ieee_is_finite(state%sapwood_mass) .and. &
                ieee_is_finite(state%heartwood_mass) .and. &
                ieee_is_finite(state%height)

  end function all_state_values_are_finite


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
                ieee_is_finite(result%structural_balance_error) .and. &
                ieee_is_finite(result%storage_balance_error) .and. &
                ieee_is_finite(result%whole_plant_balance_error)

  end function result_values_are_finite


  function height_from_total_stem_carbon(params, stem_carbon_total) result(height)

    type(Parameters), intent(in) :: params
    real(real64), intent(in) :: stem_carbon_total
    real(real64) :: height
    real(real64) :: height_power_exponent
    real(real64) :: height_power
    real(real64) :: pi_over_four

    if (stem_carbon_total <= 0.0_real64 .or. params%wood_density <= 0.0_real64) then
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


  function pipe_residual_for_state(params, state) result(pipe_residual)

    type(Parameters), intent(in) :: params
    type(PlantCarbonState), intent(in) :: state
    real(real64) :: pipe_residual
    real(real64) :: leaf_area
    real(real64) :: sapwood_area_from_mass

    leaf_area = state%leaf_mass * params%sla

    if (state%height > 0.0_real64 .and. params%wood_density > 0.0_real64) then
       sapwood_area_from_mass = state%sapwood_mass / &
                                (params%wood_density * state%height)
    else
       sapwood_area_from_mass = 0.0_real64
    end if

    pipe_residual = leaf_area - params%latosa * sapwood_area_from_mass

  end function pipe_residual_for_state


  function leaf_root_residual_for_state(params, state) result(leaf_root_residual)

    type(Parameters), intent(in) :: params
    type(PlantCarbonState), intent(in) :: state
    real(real64) :: leaf_root_residual

    leaf_root_residual = state%leaf_mass - &
                         params%leaf_to_root_ratio * state%root_mass

  end function leaf_root_residual_for_state


  function storage_fraction(state, carbon_storage) result(fraction)

    type(PlantCarbonState), intent(in) :: state
    real(real64), intent(in) :: carbon_storage
    real(real64) :: fraction
    real(real64) :: total_carbon

    total_carbon = state%leaf_mass + state%root_mass + state%sapwood_mass + &
                   state%heartwood_mass + carbon_storage

    if (total_carbon > 0.0_real64) then
       fraction = carbon_storage / total_carbon
    else
       fraction = 0.0_real64
    end if

  end function storage_fraction


  function solve_sapwood_for_pipe_balance(params, leaf_mass, heartwood_mass) result(sapwood_mass)

    type(Parameters), intent(in) :: params
    real(real64), intent(in) :: leaf_mass
    real(real64), intent(in) :: heartwood_mass
    real(real64) :: sapwood_mass
    real(real64) :: left
    real(real64) :: right
    real(real64) :: mid
    real(real64) :: f_mid
    integer :: iter

    ! Solve for the sapwood mass that makes the pipe-model residual equal to
    ! zero for a given leaf mass and heartwood mass.
    left = 1.0e-12_real64
    right = 1.0_real64

    do while (pipe_balance_function(params, leaf_mass, heartwood_mass, right) > 0.0_real64)
       right = right * 2.0_real64
       if (right > 1.0e12_real64) then
          error stop "Could not bracket sapwood mass for pipe balance."
       end if
    end do

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


  function pipe_balance_function(params, leaf_mass, heartwood_mass, sapwood_mass) result(residual)

    type(Parameters), intent(in) :: params
    real(real64), intent(in) :: leaf_mass
    real(real64), intent(in) :: heartwood_mass
    real(real64), intent(in) :: sapwood_mass
    real(real64) :: residual
    real(real64) :: height
    real(real64) :: leaf_area
    real(real64) :: sapwood_area_from_mass

    height = height_from_total_stem_carbon(params, sapwood_mass + heartwood_mass)
    leaf_area = leaf_mass * params%sla

    if (height > 0.0_real64) then
       sapwood_area_from_mass = sapwood_mass / (params%wood_density * height)
    else
       sapwood_area_from_mass = 0.0_real64
    end if

    residual = leaf_area - params%latosa * sapwood_area_from_mass

  end function pipe_balance_function


  subroutine write_csv_header(csv_unit)

    integer, intent(in) :: csv_unit

    write(csv_unit,'(a)') &
      "scenario_id,status,fail_day,failure_reason,trait_case,state_case,background_mode," // &
      "npp_rate,initial_storage,allometric_adjustment_days,max_allocation_fraction," // &
      "initial_leaf,initial_root,initial_sapwood,initial_heartwood,initial_height," // &
      "final_leaf,final_root,final_sapwood,final_heartwood,final_height,final_storage," // &
      "cumulative_npp,cumulative_structural_allocation,cumulative_unmet_storage_deficit," // &
      "final_leaf_root_residual,final_pipe_residual,final_storage_fraction," // &
      "structural_fraction_of_positive_npp,max_abs_structural_balance_error," // &
      "max_abs_storage_balance_error,max_abs_whole_plant_balance_error," // &
      "max_storage_fraction,min_living_pool,max_daily_allocation_observed," // &
      "max_daily_demand_observed,days_with_allocation,days_with_unmet_storage_deficit"

  end subroutine write_csv_header


  subroutine write_csv_row(csv_unit, summary)

    integer, intent(in) :: csv_unit
    type(ScenarioSummary), intent(in) :: summary
    character(len=8) :: status

    if (summary%passed) then
       status = "PASS"
    else
       status = "FAIL"
    end if

    write(csv_unit,'(*(g0))') &
      summary%scenario_id, ",", trim(status), ",", summary%fail_day, ",", trim(summary%failure_reason), ",", &
      summary%trait_case, ",", summary%state_case, ",", summary%background_mode, ",", &
      summary%npp_rate, ",", summary%initial_storage, ",", summary%allometric_adjustment_days, ",", &
      summary%max_allocation_fraction, ",", &
      summary%initial_leaf, ",", summary%initial_root, ",", summary%initial_sapwood, ",", &
      summary%initial_heartwood, ",", summary%initial_height, ",", &
      summary%final_leaf, ",", summary%final_root, ",", summary%final_sapwood, ",", &
      summary%final_heartwood, ",", summary%final_height, ",", summary%final_storage, ",", &
      summary%cumulative_npp, ",", summary%cumulative_structural_allocation, ",", &
      summary%cumulative_unmet_storage_deficit, ",", &
      summary%final_leaf_root_residual, ",", summary%final_pipe_residual, ",", &
      summary%final_storage_fraction, ",", summary%structural_fraction_of_positive_npp, ",", &
      summary%max_abs_structural_balance_error, ",", summary%max_abs_storage_balance_error, ",", &
      summary%max_abs_whole_plant_balance_error, ",", summary%max_storage_fraction, ",", &
      summary%min_living_pool, ",", summary%max_daily_allocation_observed, ",", &
      summary%max_daily_demand_observed, ",", summary%days_with_allocation, ",", &
      summary%days_with_unmet_storage_deficit

  end subroutine write_csv_row

end program test_storage_allocation_sensitivity
