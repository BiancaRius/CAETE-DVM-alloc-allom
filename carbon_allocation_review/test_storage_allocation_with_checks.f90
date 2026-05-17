program test_storage_allocation

  use, intrinsic :: iso_fortran_env, only: real64
  use carbon_allocation_offline_kernel

  implicit none

  type(Parameters) :: params
  type(StorageAllocationControls) :: controls
  type(PlantCarbonState) :: state
  type(AllocationOutput) :: result

  integer :: day
  real(real64) :: npp_rate
  real(real64) :: carbon_storage
  real(real64) :: stem_carbon_total_initial
  real(real64) :: height_power_exponent
  real(real64) :: height_power_initial
  real(real64) :: pi_over_four

  !--------------------------------------------------------------------------
  ! Parameter values for a first numerical test.
  ! Units must be consistent across all carbon pools and allometric parameters.
  !--------------------------------------------------------------------------
  params%sla                = 12.0_real64
  params%latosa             = 8000.0_real64
  params%wood_density       = 250.0_real64
  params%leaf_to_root_ratio = 1.0_real64
  params%allom2             = 40.0_real64
  params%allom3             = 0.5_real64

  !--------------------------------------------------------------------------
  ! Initial carbon state before daily allocation.
  !--------------------------------------------------------------------------
  state%leaf_mass      = 1.0_real64
  state%root_mass      = 0.8_real64
  state%sapwood_mass   = 10.0_real64
  state%heartwood_mass = 20.0_real64

  ! Compute initial height from total stem carbon and height-diameter allometry.
  stem_carbon_total_initial = state%sapwood_mass + state%heartwood_mass
  height_power_exponent = 1.0_real64 + 2.0_real64 / params%allom3
  pi_over_four = acos(-1.0_real64) / 4.0_real64
  height_power_initial = params%allom2**(2.0_real64 / params%allom3) * &
                         (stem_carbon_total_initial / params%wood_density) / &
                         pi_over_four
  state%height = height_power_initial**(1.0_real64 / height_power_exponent)

  !--------------------------------------------------------------------------
  ! Gradual allocation controls.
  !--------------------------------------------------------------------------
  controls%dt_years = 1.0_real64 / 365.0_real64
  controls%allometric_adjustment_days = 365.0_real64
  controls%max_allocation_fraction = 0.005_real64
  controls%leaf_background_timescale_years = 3.0_real64
  controls%root_background_timescale_years = 3.0_real64
  controls%sapwood_background_timescale_years = 15.0_real64

  ! Annualized NPP rate used every day in this simple test.
  npp_rate = 3.5_real64
  carbon_storage = 0.0_real64

  write(*,'(a)') " day      leaf      root   sapwood    height   storage      alloc    LR_resid   pipe_resid   Cbal_err"

  do day = 1, 365

     call allocate_gradual_with_storage(state, params, controls, npp_rate, carbon_storage, result)

     call check_carbon_accounting(day, result, 1.0e-10_real64)

     ! Update the state with the gradual-allocation result.
     state%leaf_mass = result%leaf_mass_new
     state%root_mass = result%root_mass_new
     state%sapwood_mass = result%sapwood_mass_new
     state%heartwood_mass = result%heartwood_mass_new
     state%height = result%height_new

     ! Print monthly diagnostics plus the first day.
     if (day == 1 .or. mod(day, 30) == 0 .or. day == 365) then
        write(*,'(i4,9es11.3)') day, state%leaf_mass, state%root_mass, &
             state%sapwood_mass, state%height, carbon_storage, &
             result%carbon_to_allocate, result%leaf_root_residual, &
             result%pipe_model_residual, result%whole_plant_balance_error
     end if

  end do

contains

  subroutine check_carbon_accounting(day, result, tolerance)

    integer, intent(in) :: day
    type(AllocationOutput), intent(in) :: result
    real(real64), intent(in) :: tolerance

    ! Stop immediately if the structural allocation does not close.
    if (abs(result%structural_balance_error) > tolerance) then
       write(*,'(a,i0)') "Structural carbon accounting failed at day ", day
       write(*,'(a,es14.6)') "  structural_balance_error = ", result%structural_balance_error
       error stop "Structural carbon accounting failed."
    end if

    ! Stop immediately if the storage update does not close.
    if (abs(result%storage_balance_error) > tolerance) then
       write(*,'(a,i0)') "Storage carbon accounting failed at day ", day
       write(*,'(a,es14.6)') "  storage_balance_error = ", result%storage_balance_error
       error stop "Storage carbon accounting failed."
    end if

    ! Stop immediately if the whole-plant accounting does not close.
    if (abs(result%whole_plant_balance_error) > tolerance) then
       write(*,'(a,i0)') "Whole-plant carbon accounting failed at day ", day
       write(*,'(a,es14.6)') "  whole_plant_balance_error = ", result%whole_plant_balance_error
       error stop "Whole-plant carbon accounting failed."
    end if

    ! Warn when negative NPP exceeded available storage. This is not a numerical
    ! error, but it means that another process must close the missing carbon sink.
    if (result%unmet_storage_deficit > tolerance) then
       write(*,'(a,i0,a,es14.6)') "Warning at day ", day, &
          ": unmet storage deficit = ", result%unmet_storage_deficit
    end if

    ! These guards catch non-physical outputs early.
    if (result%carbon_storage_after < -tolerance) error stop "Negative storage detected."
    if (result%delta_leaf < -tolerance) error stop "Negative leaf increment detected."
    if (result%delta_root < -tolerance) error stop "Negative root increment detected."
    if (result%delta_sapwood < -tolerance) error stop "Negative sapwood increment detected."

  end subroutine check_carbon_accounting

end program test_storage_allocation
