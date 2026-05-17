program test_storage_allocation_unit_cases

  use, intrinsic :: iso_fortran_env, only: real64
  use carbon_allocation_offline_kernel

  implicit none

  real(real64), parameter :: tol = 1.0e-10_real64
  real(real64), parameter :: loose_tol = 1.0e-8_real64

  type(Parameters) :: params
  type(StorageAllocationControls) :: controls

  call initialize_default_parameters(params)
  call initialize_default_controls(controls)

  write(*,'(a)') "Running unit tests for gradual storage allocation..."

  call test_zero_npp_zero_storage(params, controls)
  call test_positive_npp_imbalanced_plant(params, controls)
  call test_equilibrated_plant_with_background_growth(params, controls)
  call test_equilibrated_plant_without_background_growth(params, controls)
  call test_storage_can_support_growth_when_npp_is_zero(params, controls)
  call test_negative_npp_with_enough_storage(params, controls)
  call test_negative_npp_exceeding_storage(params, controls)
  call test_maximum_daily_allocation_limiter(params, controls)

  write(*,'(a)') "All gradual storage allocation unit tests passed."

contains

  subroutine initialize_default_parameters(params)

    type(Parameters), intent(out) :: params

    ! These values reproduce the simple numerical scenario used in the
    ! previous storage-allocation test driver.
    params%sla                = 12.0_real64
    params%latosa             = 8000.0_real64
    params%wood_density       = 250.0_real64
    params%leaf_to_root_ratio = 1.0_real64
    params%allom2             = 40.0_real64
    params%allom3             = 0.5_real64

  end subroutine initialize_default_parameters


  subroutine initialize_default_controls(controls)

    type(StorageAllocationControls), intent(out) :: controls

    ! Daily time step and gradual allometric relaxation.
    controls%dt_years = 1.0_real64 / 365.0_real64
    controls%allometric_adjustment_days = 365.0_real64
    controls%max_allocation_fraction = 0.005_real64

    ! Background balanced-growth demand. These can be set to zero in tests
    ! where we want allocation to depend only on allometric deficits.
    controls%leaf_background_timescale_years = 3.0_real64
    controls%root_background_timescale_years = 3.0_real64
    controls%sapwood_background_timescale_years = 15.0_real64

  end subroutine initialize_default_controls


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
    ! zero for a given leaf mass and heartwood mass. This helper is used only to
    ! construct a clean initial condition for unit testing.
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


  subroutine check_carbon_accounting(test_name, result)

    character(len=*), intent(in) :: test_name
    type(AllocationOutput), intent(in) :: result

    ! These checks should be true for every call to allocate_gradual_with_storage,
    ! regardless of the biological scenario being tested.
    call assert_close(test_name // ": structural carbon accounting", &
                      result%structural_balance_error, 0.0_real64, tol)

    call assert_close(test_name // ": storage carbon accounting", &
                      result%storage_balance_error, 0.0_real64, tol)

    call assert_close(test_name // ": whole-plant carbon accounting", &
                      result%whole_plant_balance_error, 0.0_real64, tol)

    call assert_true(test_name // ": no negative storage", &
                     result%carbon_storage_after >= -tol)

    call assert_true(test_name // ": no negative leaf increment", &
                     result%delta_leaf >= -tol)

    call assert_true(test_name // ": no negative root increment", &
                     result%delta_root >= -tol)

    call assert_true(test_name // ": no negative sapwood increment", &
                     result%delta_sapwood >= -tol)

  end subroutine check_carbon_accounting


  subroutine assert_true(test_name, condition)

    character(len=*), intent(in) :: test_name
    logical, intent(in) :: condition

    if (.not. condition) then
       write(*,'(a)') "FAILED: " // trim(test_name)
       error stop "Unit test failed."
    end if

  end subroutine assert_true


  subroutine assert_close(test_name, value, expected, tolerance)

    character(len=*), intent(in) :: test_name
    real(real64), intent(in) :: value
    real(real64), intent(in) :: expected
    real(real64), intent(in) :: tolerance

    if (abs(value - expected) > tolerance) then
       write(*,'(a)') "FAILED: " // trim(test_name)
       write(*,'(a,es18.8)') "  value    = ", value
       write(*,'(a,es18.8)') "  expected = ", expected
       write(*,'(a,es18.8)') "  tolerance= ", tolerance
       error stop "Unit test failed."
    end if

  end subroutine assert_close


  subroutine print_pass(test_name)

    character(len=*), intent(in) :: test_name

    write(*,'(a)') "PASS: " // trim(test_name)

  end subroutine print_pass


  subroutine test_zero_npp_zero_storage(params, controls)

    type(Parameters), intent(in) :: params
    type(StorageAllocationControls), intent(in) :: controls
    type(PlantCarbonState) :: state
    type(AllocationOutput) :: result
    real(real64) :: carbon_storage
    real(real64) :: npp_rate

    ! Test A: with zero NPP and zero storage, the routine must not create
    ! structural growth, even if the plant is allometrically imbalanced.
    state = build_state(params, 1.0_real64, 0.8_real64, 10.0_real64, 20.0_real64)
    carbon_storage = 0.0_real64
    npp_rate = 0.0_real64

    call allocate_gradual_with_storage(state, params, controls, npp_rate, carbon_storage, result)
    call check_carbon_accounting("A zero NPP and zero storage", result)

    call assert_close("A alloc", result%carbon_to_allocate, 0.0_real64, tol)
    call assert_close("A delta leaf", result%delta_leaf, 0.0_real64, tol)
    call assert_close("A delta root", result%delta_root, 0.0_real64, tol)
    call assert_close("A delta sapwood", result%delta_sapwood, 0.0_real64, tol)
    call assert_close("A storage", carbon_storage, 0.0_real64, tol)
    call assert_close("A leaf unchanged", result%leaf_mass_new, state%leaf_mass, tol)
    call assert_close("A root unchanged", result%root_mass_new, state%root_mass, tol)
    call assert_close("A sapwood unchanged", result%sapwood_mass_new, state%sapwood_mass, tol)

    call print_pass("A zero NPP and zero storage")

  end subroutine test_zero_npp_zero_storage


  subroutine test_positive_npp_imbalanced_plant(params, controls)

    type(Parameters), intent(in) :: params
    type(StorageAllocationControls), intent(in) :: controls
    type(PlantCarbonState) :: state
    type(PlantCarbonState) :: initial_state
    type(AllocationOutput) :: result
    integer :: day
    real(real64) :: carbon_storage
    real(real64) :: npp_rate
    real(real64) :: initial_lr_residual
    real(real64) :: initial_pipe_residual

    ! Test B: with positive NPP and an imbalanced plant, structural pools should
    ! increase and allometric residuals should move toward zero over one year.
    initial_state = build_state(params, 1.0_real64, 0.8_real64, 10.0_real64, 20.0_real64)
    state = initial_state
    carbon_storage = 0.0_real64
    npp_rate = 3.5_real64

    initial_lr_residual = leaf_root_residual_for_state(params, initial_state)
    initial_pipe_residual = pipe_residual_for_state(params, initial_state)

    do day = 1, 365
       call allocate_gradual_with_storage(state, params, controls, npp_rate, carbon_storage, result)
       call check_carbon_accounting("B positive NPP imbalanced plant", result)
       state%leaf_mass = result%leaf_mass_new
       state%root_mass = result%root_mass_new
       state%sapwood_mass = result%sapwood_mass_new
       state%heartwood_mass = result%heartwood_mass_new
       state%height = result%height_new
    end do

    call assert_true("B leaf increased", state%leaf_mass > initial_state%leaf_mass)
    call assert_true("B root increased", state%root_mass > initial_state%root_mass)
    call assert_true("B sapwood increased", state%sapwood_mass > initial_state%sapwood_mass)
    call assert_true("B storage non-negative", carbon_storage >= -tol)
    call assert_true("B leaf-root residual improved", &
                     abs(result%leaf_root_residual) < abs(initial_lr_residual))
    call assert_true("B pipe residual improved", &
                     abs(result%pipe_model_residual) < abs(initial_pipe_residual))

    call print_pass("B positive NPP with imbalanced plant")

  end subroutine test_positive_npp_imbalanced_plant


  subroutine test_equilibrated_plant_with_background_growth(params, controls)

    type(Parameters), intent(in) :: params
    type(StorageAllocationControls), intent(in) :: controls
    type(PlantCarbonState) :: state
    type(AllocationOutput) :: result
    real(real64) :: carbon_storage
    real(real64) :: npp_rate
    real(real64) :: sapwood_mass

    ! Test C1: an equilibrated plant should still grow when background balanced
    ! growth demand is enabled and carbon is available.
    sapwood_mass = solve_sapwood_for_pipe_balance(params, 1.0_real64, 0.0_real64)
    state = build_state(params, 1.0_real64, 1.0_real64, sapwood_mass, 0.0_real64)
    carbon_storage = 0.0_real64
    npp_rate = 3.5_real64

    call assert_close("C1 initial leaf-root balance", &
                      leaf_root_residual_for_state(params, state), 0.0_real64, loose_tol)
    call assert_close("C1 initial pipe balance", &
                      pipe_residual_for_state(params, state), 0.0_real64, loose_tol)

    call allocate_gradual_with_storage(state, params, controls, npp_rate, carbon_storage, result)
    call check_carbon_accounting("C1 equilibrated plant with background growth", result)

    call assert_true("C1 positive allocation", result%carbon_to_allocate > 0.0_real64)
    call assert_true("C1 positive leaf growth", result%delta_leaf > 0.0_real64)
    call assert_true("C1 positive root growth", result%delta_root > 0.0_real64)
    call assert_true("C1 positive sapwood growth", result%delta_sapwood > 0.0_real64)

    call print_pass("C1 equilibrated plant with background growth")

  end subroutine test_equilibrated_plant_with_background_growth


  subroutine test_equilibrated_plant_without_background_growth(params, controls)

    type(Parameters), intent(in) :: params
    type(StorageAllocationControls), intent(in) :: controls
    type(StorageAllocationControls) :: local_controls
    type(PlantCarbonState) :: state
    type(AllocationOutput) :: result
    real(real64) :: carbon_storage
    real(real64) :: npp_rate
    real(real64) :: sapwood_mass

    ! Test C2: if an equilibrated plant has no allometric deficit and background
    ! demand is disabled, positive NPP should remain in storage rather than being
    ! converted into structural biomass.
    local_controls = controls
    local_controls%leaf_background_timescale_years = 0.0_real64
    local_controls%root_background_timescale_years = 0.0_real64
    local_controls%sapwood_background_timescale_years = 0.0_real64

    sapwood_mass = solve_sapwood_for_pipe_balance(params, 1.0_real64, 0.0_real64)
    state = build_state(params, 1.0_real64, 1.0_real64, sapwood_mass, 0.0_real64)
    carbon_storage = 0.0_real64
    npp_rate = 3.5_real64

    call allocate_gradual_with_storage(state, params, local_controls, npp_rate, carbon_storage, result)
    call check_carbon_accounting("C2 equilibrated plant without background growth", result)

    call assert_close("C2 near-zero allocation", result%carbon_to_allocate, 0.0_real64, loose_tol)
    call assert_close("C2 storage equals daily NPP", carbon_storage, result%npp_daily, loose_tol)

    call print_pass("C2 equilibrated plant without background growth")

  end subroutine test_equilibrated_plant_without_background_growth


  subroutine test_storage_can_support_growth_when_npp_is_zero(params, controls)

    type(Parameters), intent(in) :: params
    type(StorageAllocationControls), intent(in) :: controls
    type(PlantCarbonState) :: state
    type(AllocationOutput) :: result
    real(real64) :: carbon_storage_before
    real(real64) :: carbon_storage
    real(real64) :: npp_rate

    ! Test D: when storage is positive, the routine can use stored carbon for
    ! structural growth even if the current time step has zero NPP.
    state = build_state(params, 1.0_real64, 0.8_real64, 10.0_real64, 20.0_real64)
    carbon_storage_before = 1.0_real64
    carbon_storage = carbon_storage_before
    npp_rate = 0.0_real64

    call allocate_gradual_with_storage(state, params, controls, npp_rate, carbon_storage, result)
    call check_carbon_accounting("D storage supports growth with zero NPP", result)

    call assert_true("D positive allocation", result%carbon_to_allocate > 0.0_real64)
    call assert_true("D storage decreased", carbon_storage < carbon_storage_before)
    call assert_true("D leaf increased", result%leaf_mass_new > state%leaf_mass)
    call assert_true("D root increased", result%root_mass_new > state%root_mass)

    call print_pass("D storage supports growth when NPP is zero")

  end subroutine test_storage_can_support_growth_when_npp_is_zero


  subroutine test_negative_npp_with_enough_storage(params, controls)

    type(Parameters), intent(in) :: params
    type(StorageAllocationControls), intent(in) :: controls
    type(PlantCarbonState) :: state
    type(AllocationOutput) :: result
    real(real64) :: carbon_storage_before
    real(real64) :: carbon_storage
    real(real64) :: npp_rate

    ! Test E1: negative NPP should first consume labile storage. If enough
    ! storage remains and structural demand is positive, allocation can still
    ! occur, but the reported unmet deficit must remain zero.
    state = build_state(params, 1.0_real64, 0.8_real64, 10.0_real64, 20.0_real64)
    carbon_storage_before = 1.0_real64
    carbon_storage = carbon_storage_before
    npp_rate = -3.65_real64

    call allocate_gradual_with_storage(state, params, controls, npp_rate, carbon_storage, result)
    call check_carbon_accounting("E1 negative NPP with enough storage", result)

    call assert_close("E1 no unmet deficit", result%unmet_storage_deficit, 0.0_real64, tol)
    call assert_true("E1 storage decreased", carbon_storage < carbon_storage_before)
    call assert_true("E1 allocation non-negative", result%carbon_to_allocate >= -tol)

    call print_pass("E1 negative NPP with enough storage")

  end subroutine test_negative_npp_with_enough_storage


  subroutine test_negative_npp_exceeding_storage(params, controls)

    type(Parameters), intent(in) :: params
    type(StorageAllocationControls), intent(in) :: controls
    type(PlantCarbonState) :: state
    type(AllocationOutput) :: result
    real(real64) :: carbon_storage
    real(real64) :: npp_rate

    ! Test E2: if negative NPP is larger than available storage, storage is
    ! clamped at zero and the missing carbon is reported as unmet_storage_deficit.
    state = build_state(params, 1.0_real64, 0.8_real64, 10.0_real64, 20.0_real64)
    carbon_storage = 0.001_real64
    npp_rate = -3.65_real64

    call allocate_gradual_with_storage(state, params, controls, npp_rate, carbon_storage, result)
    call check_carbon_accounting("E2 negative NPP exceeding storage", result)

    call assert_true("E2 positive unmet deficit", result%unmet_storage_deficit > 0.0_real64)
    call assert_close("E2 zero final storage", carbon_storage, 0.0_real64, tol)
    call assert_close("E2 zero allocation", result%carbon_to_allocate, 0.0_real64, tol)

    call print_pass("E2 negative NPP exceeding storage")

  end subroutine test_negative_npp_exceeding_storage


  subroutine test_maximum_daily_allocation_limiter(params, controls)

    type(Parameters), intent(in) :: params
    type(StorageAllocationControls), intent(in) :: controls
    type(StorageAllocationControls) :: local_controls
    type(PlantCarbonState) :: state
    type(AllocationOutput) :: result
    real(real64) :: carbon_storage
    real(real64) :: npp_rate
    real(real64) :: expected_maximum_allocation

    ! Test F: when storage and demand are both high, the allocation must be
    ! capped by max_allocation_fraction times living structural carbon.
    local_controls = controls
    local_controls%max_allocation_fraction = 1.0e-6_real64

    state = build_state(params, 1.0_real64, 0.8_real64, 10.0_real64, 20.0_real64)
    carbon_storage = 1000.0_real64
    npp_rate = 1000.0_real64

    expected_maximum_allocation = local_controls%max_allocation_fraction * &
                                  (state%leaf_mass + state%root_mass + state%sapwood_mass)

    call allocate_gradual_with_storage(state, params, local_controls, npp_rate, carbon_storage, result)
    call check_carbon_accounting("F maximum daily allocation limiter", result)

    call assert_close("F allocation equals limiter", &
                      result%carbon_to_allocate, expected_maximum_allocation, tol)

    call assert_close("F stored maximum allocation diagnostic", &
                      result%max_daily_allocation, expected_maximum_allocation, tol)

    call print_pass("F maximum daily allocation limiter")

  end subroutine test_maximum_daily_allocation_limiter

end program test_storage_allocation_unit_cases
