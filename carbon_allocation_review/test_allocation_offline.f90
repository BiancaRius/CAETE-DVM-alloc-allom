program test_allocation

  use, intrinsic :: iso_fortran_env, only: real64
  use carbon_allocation_offline_kernel

  implicit none

  type(Parameters) :: params
  type(PlantCarbonState)     :: state
  type(AllocationOutput)     :: result

  real(real64) :: c_available
  real(real64) :: sapwood_area_initial

  !---------------------------------------------------------------------------
  ! Parameter values for a first numerical test.
  ! These values are only for testing the solver, not for final model use.
  ! Carbon pools are assumed to be in kgC per average individual.
  !---------------------------------------------------------------------------

  params%sla                = 12.0_real64
  params%latosa             = 8000.0_real64
  params%wood_density       = 250.0_real64
  params%leaf_to_root_ratio = 1.0_real64
  params%allom2             = 40.0_real64
  params%allom3             = 0.5_real64

  !---------------------------------------------------------------------------
  ! Initial carbon state before allocation.
  !---------------------------------------------------------------------------

  state%leaf_mass      = 1.0_real64
  state%root_mass      = 0.8_real64
  state%sapwood_mass   = 10.0_real64
  state%heartwood_mass = 20.0_real64

  ! Compute initial sapwood area from the pipe model ( to calculate initial height)
  sapwood_area_initial = state%leaf_mass * params%sla / params%latosa

  ! Compute initial plant height from the pipe model.
  state%height = state%sapwood_mass / &
    (params%wood_density * sapwood_area_initial)

  ! Carbon available for allocation over this test period.
  c_available = 5.0_real64

  call allocate(state, params, c_available, result)

  write(*,'(a,l1)')       "normal_allocation = ", result%normal_allocation
  write(*,'(a,l1)')       "converged         = ", result%converged
  write(*,'(a,i0)')       "iterations        = ", result%iterations

  write(*,'(a,es16.8)')   "delta_leaf        = ", result%delta_leaf
  write(*,'(a,es16.8)')   "delta_root        = ", result%delta_root
  write(*,'(a,es16.8)')   "delta_sapwood     = ", result%delta_sapwood

  write(*,'(a,es16.8)')   "leaf_new          = ", result%leaf_mass_new
  write(*,'(a,es16.8)')   "root_new          = ", result%root_mass_new
  write(*,'(a,es16.8)')   "sapwood_new       = ", result%sapwood_mass_new
  write(*,'(a,es16.8)')   "heartwood_new     = ", result%heartwood_mass_new

  write(*,'(a,es16.8)')   "height_old        = ", state%height
  write(*,'(a,es16.8)')   "height_new        = ", result%height_new
  write(*,'(a,es16.8)')   "diameter_new      = ", result%stem_diameter_new

  write(*,'(a,es16.8)')   "carbon_error      = ", result%carbon_balance_error
  write(*,'(a,es16.8)')   "leaf_root_resid   = ", result%leaf_root_residual
  write(*,'(a,es16.8)')   "pipe_resid        = ", result%pipe_model_residual
  write(*,'(a,es16.8)')   "alloc_resid       = ", result%allocation_residual_final

  write(*,'(a,a)')        "message           = ", trim(result%message)

end program test_allocation