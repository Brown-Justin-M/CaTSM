module rkcn2_module
  !/**********************************************************************|
  !| RKCN2 MODULE                                                         |
  !| ------------                                                         |
  !| This module contains the implementation of the Runge--Kutta/Crank--  |
  !| Nicolson scheme to update the state of the system using the prior    |
  !| states and calculated rhs of the system of equations.                |
  !\**********************************************************************/
  use grid_module, only: dp,li,grid_type
  use params_module, only: params_type
  use fields_module, only: state_type,rhs_type,field_type
  implicit none
  private
  public rkcn2_dt_limit,rkcn2_run

contains
  !/**********************************************************************|
  !| FUNCTIONS                                                            |
  !\**********************************************************************/ 
  function rkcn2_dt_limit(grid,state,diff)
    !/********************************************************************|
    !| RKCN2_DT_LIMIT                                                     |
    !| --------------                                                     |
    !| Calculates the limit on dt to ensure numerical stability. This is  |
    !| adapted from the 1D version described in 4.1.1(e) of Peyret (2002) |
    !| Spectral Methods for Incompressible Viscous Flow, Equation 4.97.   |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   state: a state_type object containing the simulation data        |
    !|   diff: the smallest diffusion coeffient for the equations         |
    !|   returns: the maximum length of the next timestep                 |
    !\********************************************************************/
    use equations_module, only: dt_limit
    implicit none

    type(grid_type), intent(in) :: grid
    type(state_type), intent(in) :: state
    real(kind=dp), intent(in) :: diff
    real(kind=dp) :: rkcn2_dt_limit
    
    ! TODO: These numbers are currently just conservative estimates
    rkcn2_dt_limit = dt_limit(grid,state,diff,1.0_dp,0.5_dp)
  end function

  !/**********************************************************************|
  !| SUBROUTINES                                                          |
  !\**********************************************************************/ 
  subroutine rkcn2_run(grid,params,state,prior,rhs)  
    !/********************************************************************|
    !| RKCN2_RUN                                                          |
    !| ---------                                                          |
    !| Run the RKCN2 scheme, update the rhs and the state using prior     |
    !| values of the rhs and state using the following equation           |
    !|   u* = u(n) + dt * (E(n) + I(n)),                                  |
    !|   2 * u(n+1) = u(n) + u* + dt * (E(*) + I(n+1)),                   |
    !| where E(i) is the explicit term of the ith step and I(i) is the    |
    !| implicit term of the ith step.                                     |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   prior: an array of prior state_type objects from previous states |
    !|   rhs: an array of rhs_type objects current and past               |
    !\********************************************************************/
    use comm_module, only: get_id
    use grid_module, only: allocate_s,allocate_p,allocate_s_vector
    use fields_module, only: init_state,copy_state,free_state
    use fields_module, only: init_rhs,free_rhs
    use fields_module, only: step_time,transform_all_to_p
    use transform_module, only: transform_s_to_p
    use equations_module, only: explicit_rhs,explicit_implicit_rhs,implicit_update,update_eos
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(inout) :: state
    type(state_type), intent(in) :: prior(:)
    type(rhs_type), intent(inout) :: rhs(:)

    type(state_type) :: mid_state
    type(rhs_type) :: mid_rhs
    complex(kind=dp), dimension(:,:,:), pointer :: work_scalar
    complex(kind=dp), dimension(:,:,:,:), pointer :: work_vector

    real (kind=dp), pointer :: work_p(:,:,:)

    call init_state(grid,mid_state)
    call copy_state(prior(1),mid_state)
    call init_rhs(grid,mid_rhs)

    mid_state%dt = state%dt
    mid_state%time = prior(1)%time + mid_state%dt

    call explicit_rhs(grid,params,prior(1),mid_state,rhs(1))
    call explicit_implicit_rhs(grid,params,prior(1),mid_state,mid_rhs)

    mid_state%all_s = prior(1)%all_s + mid_state%dt * (rhs(1)%all + mid_rhs%all)

    call transform_all_to_p(grid,mid_state)
    call update_eos(grid,params,mid_state)

    state%all_s = prior(1)%all_s + mid_state%all_s

    call explicit_rhs(grid,params,mid_state,state,mid_rhs)
    call free_state(mid_state)

    state%all_s = state%all_s + state%dt * mid_rhs%all
    call free_rhs(mid_rhs)

    call implicit_update(grid,params,state,2.0_dp)
  end subroutine

end module