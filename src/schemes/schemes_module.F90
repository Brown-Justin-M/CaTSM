module schemes_module
  !/**********************************************************************|
  !| SCHEMES MODULE                                                       |
  !| --------------                                                       |
  !| This module contains the interface and selectors to run the various  |
  !| numerical schemes made available by the program. The primary         |
  !| available interface is the update subroutine, which steps the        |
  !| program forward one timestep.                                        |
  !\**********************************************************************/
  use grid_module, only: dp,li,grid_type
  use params_module, only: params_type,scheme_rkcn2,scheme_abbdf3
  use fields_module, only: state_type,rhs_type,field_type
  implicit none
  private
  public update

contains
  !/**********************************************************************|
  !| FUNCTIONS                                                            |
  !\**********************************************************************/
  function get_dt(grid,params,state,scheme)
    !/********************************************************************|
    !| GET_DT                                                             |
    !| ------                                                             |
    !| Calculates the length of the next timestep, using the limit        |
    !| appropriate to the chosen scheme.                                  |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   scheme: the numerical scheme                                     |
    !|   returns: the length of the next timestep                         |
    !\********************************************************************/
    use rkcn2_module, only: rkcn2_dt_limit
    use abbdf3_module, only: abbdf3_dt_limit
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    integer, intent(in) :: scheme
    real(kind=dp) :: get_dt

    real(kind=dp) :: diff

    diff = min(params%diff_t,params%diff_c,params%diff_vel)
    get_dt = state%dt

    ! Calculate the limit appropriate to scheme
    select case (scheme)
    case (scheme_rkcn2)
      get_dt = params%limit_mult_dt * rkcn2_dt_limit(grid,state,diff)
    case (scheme_abbdf3)
      get_dt = params%limit_mult_dt * abbdf3_dt_limit(grid,state,diff)
    end select

    ! The timestep can only increase by a fixed factor each step
    if (get_dt .le. params%max_increase_dt * state%dt) then
      get_dt = min(get_dt,params%max_dt)
    else
      get_dt = min(params%max_increase_dt * state%dt,params%max_dt)
    end if
  end function

  !/**********************************************************************|
  !| SUBROUTINES                                                          |
  !\**********************************************************************/
  subroutine update(grid,params,state,prior,rhs)
    !/********************************************************************|
    !| UPDATE                                                             |
    !| ------                                                             |
    !| Update the simulation to the next simulation state by choosing the |
    !| appropriate scheme from the parameters and executing run_scheme.   |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   prior: the prior state_type objects for the previous times       |
    !|   rhs: an array of rhs_type objects for prior times                |
    !\********************************************************************/
    type(grid_type), intent(inout) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(inout) :: state
    type(state_type), intent(inout) :: prior(:)
    type(rhs_type), intent(inout) :: rhs(:)

    state%timestep = state%timestep + 1

    ! Choose the appropriate scheme.
    if (state%timestep - state%restart .lt. 3) then
      call run_scheme(grid,params,state,prior,rhs,params%start_schemein)
    else 
      call run_scheme(grid,params,state,prior,rhs,params%schemein)
    end if
  end subroutine

  subroutine run_scheme(grid,params,state,prior,rhs,scheme)
    !/********************************************************************|
    !| RUN_SCHEME                                                         |
    !| ----------                                                         |
    !| Update the simulation to the next simulation state.                |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   prior: the prior state_type objects for the previous times       |
    !|   rhs: an array of rhs_type objects for prior times                |
    !|   scheme: the scheme to be run                                     |
    !\********************************************************************/
    use comm_module, only: get_id,abort_comm
    use fields_module, only: step_time,transform_all_to_p
    use equations_module, only: update_eos
    use rkcn2_module, only: rkcn2_run
    use abbdf3_module, only: abbdf3_run
    implicit none

    type(grid_type), intent(inout) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(inout) :: state
    type(state_type), intent(inout) :: prior(:)
    type(rhs_type), intent(inout) :: rhs(:)
    integer, intent(in) :: scheme

    real(kind=dp) :: rtmp

    ! If the timestep hasn't been set yet, set it now
    if (state%dt == 0.0_dp) state%dt = params%start_dt

    ! Cycle the arrays to the next time step and copy 
    call step_time(state,prior,rhs)

    ! Calculate the timestep limit, based on the scheme stability region
    state%dt = get_dt(grid,params,state,scheme)
    state%time = prior(1)%time + state%dt

    ! Update the simulation in spectral space
    select case (scheme)
    case (scheme_rkcn2)
      call rkcn2_run(grid,params,state,prior,rhs)
    case (scheme_abbdf3)
      call abbdf3_run(grid,params,state,prior,rhs)
    end select

    call filter_fields(grid,params,state,prior,rhs)

    ! Check if the simulation has crashed
    rtmp = sum(abs(state%all_s))
    if (rtmp /= rtmp) then
      print*,"FATAL: NaN detected, ending simulation."
      call abort_comm
    end if
    if (rtmp > huge(1.0_dp) .or. rtmp < -huge(1.0_dp)) then
      print*,"FATAL: +/- infinity detected, ending simulation."
      call abort_comm
    end if

    ! Update the boundary conditions if they exist
    call enforce_bounds(grid,params,state)

    ! Update the physical space representation of the state
    call transform_all_to_p(grid,state)
    call update_eos(grid,params,state)
  end subroutine

  subroutine filter_fields(grid,params,state,prior,rhs)
    !/********************************************************************|
    !| FILTER_FIELDS                                                      |
    !| -------------                                                      |
    !| Using an exponential filter of order params%order_filter, filter   |
    !| the current fields and rhs fields of the system every              |
    !| params%steps_filter timesteps. This prevents the development of    |
    !| high-spectral-order features that can arise in spectral solutions  |
    !| of the compressible equations.                                     |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   prior: the prior state_type objects for the previous times       |
    !|   rhs: an array of rhs_type objects for prior times                |
    !\********************************************************************/
    use grid_module, only: pi
    implicit none

    type(grid_type), intent(inout) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(inout) :: state
    type(state_type), intent(inout) :: prior(:)
    type(rhs_type), intent(inout) :: rhs(:)

    real(kind=dp) :: alpha,filter,xfactor,zfactor,hkx,hky,hkz
    integer :: i,j,k,l,m
    logical :: ltmp

    if (params%steps_filter <= 0) then  
      return
    end if

    if (modulo(state%timestep,params%steps_filter) /= 0) then
      return
    end if

    ! Filter the timestep occasionally using an exponential filter
    alpha = -log(1.e-14_dp) / pi**params%order_filter
    xfactor = grid%x_length / grid%x_n / 2.0_dp
    zfactor = grid%z_length / grid%z_n / 2.0_dp

    do l = 1,state%n_vars
      do j = 1,grid%y_s_count
        hky = abs(grid%lky(j))
        do i = 1,grid%x_s_count
          hkx = abs(grid%lkx(i))
          if (hkx == 0.0_dp) cycle
          
          do k = 1,2 * grid%z_modes
            hkz = abs(grid%lkz(k))

            filter = exp(-alpha * (hkx * xfactor) ** params%order_filter)
            filter = filter * exp(-alpha * (hkz * zfactor) ** params%order_filter)
    
            state%all_s(k,i,j,l) = filter * state%all_s(k,i,j,l)
            do m = 1, size(prior)
              prior(m)%all_s(k,i,j,l) = filter * prior(m)%all_s(k,i,j,l)
            end do

            do m = 1, size(rhs)
              rhs(m)%all(k,i,j,l) = filter * rhs(m)%all(k,i,j,l)
            end do
          end do
        end do
      end do
    end do
  end subroutine

  subroutine enforce_bounds(grid,params,state)
    !/********************************************************************|
    !| ENFORCE_BOUNDS                                                     |
    !| --------------                                                     |
    !| Use the specified constraints from the parameter file to enforce   |
    !| the boundaries of the system. This includes options for mirror     |
    !| symmetry in z (params%reflection) and fixed boundaries at the      |
    !| top and bottom of the reflected domain (params%bound_record). The  |
    !| fixed boundaries are recorded from the initial conditions of the   |
    !| simulation. If the system has mirror symmetry, the vertical net    |
    !| flow is set to zero. In addition, the net horizontal flow may also |
    !| be set to zero (params%no_horizontal_flow).                        |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !\********************************************************************/
    use comm_module, only: get_id
    use grid_module, only: pi
    use transform_module, only: transform_z_p_to_s,transform_z_s_to_p
    implicit none

    type(grid_type), intent(inout) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(inout) :: state

    complex(kind=dp) :: cwork(grid%z_n,grid%x_s_count,grid%y_s_count)

    integer :: k1,k2,k3,i,j,k,l,m

    ! If the domain should have mirror symmetry, reflect it in Fourier space
    if (params%reflection) then
      do l = 1, state%n_vars
        if (associated(state%rhow%s,state%all_s(:,:,:,l)) .or. &
            associated(state%w%s,state%all_s(:,:,:,l))) then
          state%all_s(2:grid%z_modes,:,:,l) = -state%all_s(2*grid%z_modes:grid%z_modes+2:-1,:,:,l)
        else
          state%all_s(2:grid%z_modes,:,:,l) = state%all_s(2*grid%z_modes:grid%z_modes+2:-1,:,:,l)
        end if
      end do

      if (get_id() == 0) then
        state%rhow%s(1,1,1) = 0.0_dp
        state%w%s(1,1,1) = 0.0_dp
      end if
    end if

    ! If no horizontal flow is specified, zero the constant Fourier mode of u
    if (params%no_horizontal_flow .and. get_id() == 0) then
      state%rhou%s(:,1,1) = 0.0_dp
      state%u%s(:,1,1) = 0.0_dp
    end if

    ! Record and set the bounds of the entropy field
    k1 = grid%z_n/2-params%bound_n
    k2 = grid%z_n/2+params%bound_n
    k3 = grid%z_n-params%bound_n

    if (params%bound_record) then
      call transform_z_s_to_p(state%en%s,cwork)
      if (get_id() == 0 .and. state%init_en_setup) then
        state%init_en = cwork(:,1,1)
        state%init_en_setup = .false.
      end if
      if (get_id() == 0 ) then
        cwork(k1:k2,1,1) = state%init_en(k1:k2)
        cwork(1:params%bound_n,1,1) = state%init_en(1:params%bound_n)
        cwork(k3:grid%z_n,1,1) = state%init_en(k3:grid%z_n)
      end if
      call transform_z_p_to_s(cwork,state%en%s)
    end if

    ! Record and set the bounds of the chemical field 
    if (params%bound_record) then
      call transform_z_s_to_p(state%c%s,cwork)
      if (get_id() == 0 .and. state%init_c_setup) then
        state%init_c = cwork(:,1,1)
        state%init_c_setup = .false.
      end if
      if (get_id() == 0 ) then
        cwork(k1:k2,1,1) = state%init_c(k1:k2)
        cwork(1:params%bound_n,1,1) = state%init_c(1:params%bound_n)
        cwork(k3:grid%z_n,1,1) = state%init_c(k3:grid%z_n)
      end if
      call transform_z_p_to_s(cwork,state%c%s)
    end if
  end subroutine

end module
