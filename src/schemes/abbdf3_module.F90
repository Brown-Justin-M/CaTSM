module abbdf3_module
  !/**********************************************************************|
  !| ABBDF3 MODULE                                                        |
  !| -------------                                                        |
  !| This module contains the implementation of the Adams--Bashforth/     |
  !| Backward Difference Formula to update the state of the system using  |
  !| the prior states and calculated rhs of the system of equations.      |
  !\**********************************************************************/
  use grid_module, only: dp,li,grid_type
  use params_module, only: params_type
  use fields_module, only: state_type,rhs_type,field_type
  implicit none
  private
  public abbdf3_dt_limit,abbdf3_run

  !/**********************************************************************|
  !| DERIVED TYPES                                                        |
  !\**********************************************************************/ 
  type coefficients_type
    !/********************************************************************| 
    !| COEFFICIENTS_TYPE                                                  |
    !| -----------------                                                  |
    !| The coefficients type contains the necessary coefficients for the  |
    !| ABBDF3 solver equation. These depend on the length of the time     |
    !| steps and so need to be calculated at every time.                  |
    !\********************************************************************/ 
    real(kind=dp) :: a0,a1,a2,a3
    real(kind=dp) :: b0,b1,b2
  end type coefficients_type

contains
  !/**********************************************************************|
  !| FUNCTIONS                                                            |
  !\**********************************************************************/ 
  function abbdf3_dt_limit(grid,state,diff)
    !/********************************************************************| 
    !| ABBDF3_DT_LIMIT                                                    |
    !| ---------------                                                    |
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
    real(kind=dp) :: abbdf3_dt_limit
    
    abbdf3_dt_limit = dt_limit(grid,state,diff,5.61_dp,0.69_dp)
  end function

  !/**********************************************************************|
  !| SUBROUTINES                                                          |
  !\**********************************************************************/ 
  subroutine abbdf3_run(grid,params,state,prior,rhs)
    !/********************************************************************| 
    !| ABBDF3_RUN                                                         |
    !| ----------                                                         |
    !| Run the ABBDF3 scheme, update the rhs and the state using prior    |
    !| values of the rhs and state using the following equation           |
    !|   sum(a_j * u(n+1-j)) = sum(b_j * E(n-j)) + dt * I(n+1),           |
    !| where E(i) is the explicit term of the ith step and I(i) is the    |
    !| implicit term of the ith step.                                     |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   prior: an array of prior state_type objects from previous states |
    !|   rhs: an array of rhs_type objects current and past               |
    !\********************************************************************/
    use comm_module, only: get_id
    use fields_module, only: step_time,transform_all_to_p
    use equations_module, only: explicit_rhs,implicit_update
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(inout) :: state
    type(state_type), intent(in) :: prior(:)
    type(rhs_type), intent(inout) :: rhs(:)

    type(coefficients_type) :: coeffs
    
    call abbdf3_coeffs(state,prior,coeffs)

    state%all_s = -coeffs%a3 * prior(3)%all_s
    state%all_s = state%all_s - coeffs%a2 * prior(2)%all_s
    state%all_s = state%all_s - coeffs%a1 * prior(1)%all_s                

    state%all_s = state%all_s + coeffs%b2 * rhs(3)%all
    state%all_s = state%all_s + coeffs%b1 * rhs(2)%all

    call explicit_rhs(grid,params,prior(1),state,rhs(1))

    state%all_s = state%all_s + coeffs%b0 * rhs(1)%all

    call implicit_update(grid,params,state,coeffs%a0)
  end subroutine

  subroutine abbdf3_coeffs(state,prior,coeffs)
    !/********************************************************************| 
    !| ABBDF3_COEFFS                                                      |
    !| -------------                                                      |
    !| Calculates the coefficients of the following equation              |
    !|   sum(a_j * u(n+1-j)) = sum(b_j * E(n-j)) + dt * I(n+1),           |
    !| where E(i) is the explicit term of the ith step and I(i) is the    |
    !| implicit term of the ith step. The coefficients are given in       |
    !| Equation 4.83 of Peyret (2002). (Note that b in this equation is   |
    !| actually b * dt from Equation 4.83)                                |
    !|   state: a state_type object containing the simulation data        |
    !|   prior: an array of prior state_type objects from previous states |
    !|   coeffs: a coefficients_type object that contains the coefficients|
    !\********************************************************************/
    implicit none

    type(state_type), intent(in) :: state
    type(state_type), intent(in) :: prior(:)
    type(coefficients_type), intent(out) :: coeffs

    real(kind=dp) :: r1,r2,opr1,opr1pr2,r1pr2

    r1 = prior(1)%dt / state%dt
    r2 = prior(2)%dt / state%dt

    opr1 = 1.0_dp + r1
    opr1pr2 = 1.0_dp + r1 + r2
    r1pr2 = r1 + r2

    coeffs%a0 = 1.0_dp + 1.0_dp / opr1 + 1.0_dp / opr1pr2
    coeffs%a1 = -(opr1 * opr1pr2) / (r1 * r1pr2)
    coeffs%a2 = opr1pr2 / (r1 * r2 * opr1)
    coeffs%a3 = -opr1 / (r2 * r1pr2 * opr1pr2)

    coeffs%b0 = (opr1 * opr1pr2) / (r1 * r1pr2) * state%dt
    coeffs%b1 = -opr1pr2 / (r1 * r2) * state%dt
    coeffs%b2 = opr1 / (r2 * r1pr2) * state%dt
  end subroutine

end module abbdf3_module
