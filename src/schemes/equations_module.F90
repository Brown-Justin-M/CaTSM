module equations_module
  !/**********************************************************************| 
  !| EQUATIONS MODULE                                                     |
  !| ----------------                                                     |
  !| This module contains the subroutines specific to the equations. This |
  !| includes the calculation of the explicit and implicit terms of the   |
  !| equations and the general form for the timestepping limit.           |
  !\**********************************************************************/
    use grid_module, only: dp,li,pi,grid_type
    use params_module, only: params_type
    use fields_module, only: state_type,rhs_type,field_type,vector_type
    implicit none
    private
    public dt_limit,explicit_rhs,explicit_implicit_rhs,implicit_update,update_eos
  
contains
  !/**********************************************************************| 
  !| FUNCTIONS                                                            |
  !\**********************************************************************/
  function dt_limit(grid,state,diff,a,b)
    !/********************************************************************|
    !| DT_LIMIT                                                           |
    !| --------                                                           |
    !| Calculates the limit on dt to ensure numerical stability. This is  |
    !| adapted from the 1D version described in 4.1.1(e) of Peyret (2002) |
    !| Spectral Methods for Incompressible Viscous Flow, Equation 4.97.   |
    !|   grid: a grid_type object that describes the grid                 |
    !|   state: a state_type object containing the simulation data        |
    !|   diff: the smallest diffusion coeffient for the equations         |
    !|   returns: the length of the next timestep                         |
    !\********************************************************************/
    use comm_module, only: collect_max,collect_array_max
    implicit none

    type(grid_type), intent(in) :: grid
    type(state_type), intent(in) :: state
    real(kind=dp), intent(in) :: diff
    real(kind=dp), intent(in) :: a,b
    real(kind=dp) :: dt_limit
    
    real(kind=dp) :: hkx,hky,hkz,ksquared
    real(kind=dp) :: umax(3),advect,rtmp
    integer :: i,j,k

    ! Get the maximum velocity on the grid
    umax(1) = collect_array_max(abs(state%u%p))
    umax(2) = collect_array_max(abs(state%v%p))
    umax(3) = collect_array_max(abs(state%w%p))

    
    ! Calculate the maximum denominator of Equation 4.97
    rtmp = epsilon(0.0_dp)
    do j = 1,grid%y_s_count
      hky = abs(grid%lky(j))
      do i = 1,grid%x_s_count
        hkx = abs(grid%lkx(i))
        do k = 1,grid%z_s_count
          hkz = abs(grid%lkz(k))
          ksquared = hkx**2 + hky**2 + hkz**2
          advect = hkx * umax(1) + hky * umax(2) + hkz * umax(3)
          rtmp = max(rtmp, a * advect - diff * ksquared)
        end do
      end do
    end do
    rtmp = collect_max(rtmp)

    if (rtmp .le. 0.0_dp) then
      dt_limit = huge(0.0_dp)
    else
      dt_limit = a * b / rtmp
    end if
  end function

  !/**********************************************************************|
  !| SUBROUTINES                                                          |
  !\**********************************************************************/ 
  subroutine explicit_rhs(grid,params,prior,state,rhs)
    !/********************************************************************|
    !| EXPLICICT_RHS                                                      |
    !| -------------                                                      |
    !| Calculate the explicit parts of the right hand sides of the        |
    !| equations. Store them in rhs. Note the warning in explicit_rhs_vel |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   prior: the prior state_type objects for the previous time        |
    !|   rhs: an rhs_type object for the current time                     |
    !|   b_rhs: the multiplier on the rhs in the final scheme             |
    !\********************************************************************/
    use grid_module, only: allocate_p
    use transform_module, only: transform_s_to_p,transform_p_to_s
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: prior
    type(state_type), intent(in) :: state
    type(rhs_type), intent(inout) :: rhs

    real (kind=dp), pointer :: work_p(:,:,:)

    call allocate_p(grid,work_p)

    rhs%all = 0.0_dp

    call explicit_rhs_vel(grid,params,prior,state,rhs%rhovel)
    call explicit_rhs_c(grid,params,prior,rhs%c)

    call transform_s_to_p(rhs%c,work_p)
    work_p = -prior%mu%p * work_p / prior%t%p
    call transform_p_to_s(work_p,rhs%en)
    call advection_rhs(grid,params,prior,prior%c,rhs%c)

    call explicit_rhs_rho(grid,params,prior,rhs%rho)
    call explicit_rhs_en(grid,params,prior,rhs%en)

    deallocate(work_p)
  end subroutine

  subroutine explicit_implicit_rhs(grid,params,prior,state,rhs)
    !/********************************************************************|
    !| EXPLICIT_IMPLICIT_RHS                                              |
    !| ---------------------                                              |
    !| Calculate the implicit part of the equations, but treat them as    |
    !| though they were explicit. This is useful for the Crank--Nicolson  |
    !| scheme.                                                            |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   prior: the prior state_type objects for the previous time        |
    !|   rhs: an rhs_type object for the current time                     |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: prior
    type(state_type), intent(in) :: state
    type(rhs_type), intent(out) :: rhs

    rhs%all = 0.0_dp

    if (params%implicit) then
      call div_rhs(grid,params,prior,prior%rhovel,-1.0_dp,rhs%rho)
      call gradient_rhs(grid,params,prior,prior%rho,-params%cs**2,rhs%rhovel)
    end if
  end subroutine
  
  subroutine explicit_rhs_rho(grid,params,state,out)
    !/********************************************************************|
    !| EXPLICICT_RHS_T                                                    |
    !| ---------------                                                    |
    !| Calculate the explicit parts of the right hand side of the         |
    !| temperature equation. Store it in out. Note that the state is the  |
    !| state for which the right hand side should be evaluated            |
    !| (typically the previous state).                                    |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   out: an array used to contain the rhs information                |
    !\********************************************************************/
    use comm_module, only: get_id
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    complex(kind=dp), intent(inout) :: out(:,:,:)

    if (.not. params%implicit) then
      call div_rhs(grid,params,state,state%rhovel,-1.0_dp,out)
    end if
  end subroutine

  subroutine explicit_rhs_en(grid,params,state,out)
    !/********************************************************************|
    !| EXPLICICT_RHS_T                                                    |
    !| ---------------                                                    |
    !| Calculate the explicit parts of the right hand side of the         |
    !| temperature equation. Store it in out. Note that the state is the  |
    !| state for which the right hand side should be evaluated            |
    !| (typically the previous state).                                    |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   out: an array used to contain the rhs information                |
    !\********************************************************************/
    use grid_module, only: allocate_p,allocate_s
    use comm_module, only: get_id
    use transform_module, only: transform_p_to_s,transform_s_to_p
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    complex(kind=dp), intent(inout) :: out(:,:,:)

    complex (kind=dp), pointer :: work_s(:,:,:)
    real (kind=dp), pointer :: work_p(:,:,:)

    call advection_rhs(grid,params,state,state%en,out)

    call diffusion_rhs(grid,params,state,state%cp,state%t,out,params%diff_t,state%rho%p*state%t%p)

    call epsilon_rhs(grid,params,state,state%vel,out,params%diff_vel,state%t)
  end subroutine

  subroutine explicit_rhs_c(grid,params,state,out)
    !/********************************************************************|
    !| EXPLICICT_RHS_C                                                    |
    !| ---------------                                                    |
    !| Calculate the explicit parts of the right hand side of the         |
    !| concentration equation. Store it in out. Note that the state is    |
    !| the state for which the right hand side should be evaluated        |
    !| (typically the previous state).                                    |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   out: an array used to contain the rhs information                |
    !\********************************************************************/
    use comm_module, only: get_id
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    complex(kind=dp), intent(inout) :: out(:,:,:)
 
    call diffusion_rhs(grid,params,state,state%rho,state%c,out,params%diff_c,state%rho%p)
  end subroutine

  subroutine explicit_rhs_vel(grid,params,prior,state,out)
    !/********************************************************************|
    !| EXPLICIT_RHS_VEL                                                   |
    !| ----------------                                                   |
    !| **WARNING** For now, the velocity in "state%vel%s" is expected to  |
    !| already include the scheme estimates from prior times.             |
    !| Calculate the explicit parts of the right hand side of the         |
    !| velocity equation. Store it in out. Note that the state is the     |
    !| state for which the right hand side should be evaluated (typically |
    !| the previous state).                                               |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   out: an array used to contain the rhs information                |
    !\********************************************************************/
    use grid_module, only: allocate_p,allocate_s
    use comm_module, only: get_id
    use transform_module, only: transform_p_to_s,transform_s_to_p 
    use mpi
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: prior
    type(state_type), intent(in) :: state
    complex(kind=dp), intent(inout) ::  out(:,:,:,:)

    real(kind=dp), pointer :: work_p(:,:,:)
    complex(kind=dp), pointer :: work_s(:,:,:)

    integer :: k

    call allocate_p(grid,work_p)
    call allocate_s(grid,work_s)

    ! Evaluate the RHS of the momentum components
    call advection_div_rhs(grid,params,prior,prior%rhou,out(:,:,:,1))
    call advection_div_rhs(grid,params,prior,prior%rhov,out(:,:,:,2))
    call advection_div_rhs(grid,params,prior,prior%rhow,out(:,:,:,3))

    ! Add in buoyancy due to vertical gravity component
    do k = 1,grid%z_p_count
      work_p(:,:,k) = params%buoy*prior%rho%p(:,:,k) * grid%lg(k)
    enddo
    call transform_p_to_s(work_p,work_s)
    out(:,:,:,3) = out(:,:,:,3) + work_s

    call gradient_rhs(grid,params,prior,prior%p,-1.0_dp,out)
    if (params%implicit) then
      call gradient_rhs(grid,params,prior,prior%rho,params%cs**2,out)
    end if

    ! Evaluate the RHS of the diffusion components
    call viscous_stress_rhs(grid,params,prior,prior%vel,out)

    deallocate(work_p)
    deallocate(work_s)
  end subroutine

  subroutine advection_div_rhs(grid,params,state,field,out)
    !/********************************************************************|
    !| ADVECTION_DIV_RHS                                                  |
    !| -----------------                                                  |
    !| Calculate the advection term in the right hand side of an equation |
    !| using the div (field * state%u) form.                              |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   field: the vector field for which to calculate the advection     |
    !|   out: an array used to contain the rhs information                |
    !\********************************************************************/
    use grid_module, only: allocate_p,allocate_s
    use grid_ops_module, only: deriv_x,deriv_y,deriv_z
    use transform_module, only: transform_p_to_s
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    type(field_type), intent(in) :: field
    complex(kind=dp), intent(inout) :: out(:,:,:)

    real(kind=dp), pointer :: work_p(:,:,:)
    complex (kind=dp), pointer :: work_s(:,:,:)

    call allocate_p(grid,work_p) 
    call allocate_s(grid,work_s)

    work_p = field%p * state%u%p
    call transform_p_to_s(work_p,work_s)
    call deriv_x(grid,work_s,work_s)
    out = out - work_s

    work_p = field%p * state%v%p
    call transform_p_to_s(work_p,work_s)
    call deriv_y(grid,work_s,work_s)
    out = out - work_s

    work_p = field%p * state%w%p
    call transform_p_to_s(work_p,work_s)
    call deriv_z(grid,work_s,work_s)
    out = out - work_s
    
    deallocate(work_p,work_s)
  end subroutine

  subroutine advection_rhs(grid,params,state,field,out)
    !/********************************************************************|
    !| ADVECTION_RHS                                                      |
    !| -------------                                                      |
    !| Calculate the advection term in the right hand side of an equation |
    !| using the state%u dot grad(field) form.                            |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   field: the field for which to calculate the advection            |
    !|   out: an array used to contain the rhs information                |
    !\********************************************************************/
    use grid_module, only: allocate_p,allocate_s
    use grid_ops_module, only: deriv_x,deriv_y,deriv_z
    use transform_module, only: transform_s_to_p,transform_p_to_s
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    type(field_type), intent(in) :: field
    complex(kind=dp), intent(inout) :: out(:,:,:)

    real(kind=dp), pointer :: work_p(:,:,:)
    complex (kind=dp), pointer :: work_s(:,:,:)

    call allocate_p(grid,work_p) 
    call allocate_s(grid,work_s)

    call deriv_x(grid,field%s,work_s)
    call transform_s_to_p(work_s,work_p)
    work_p = work_p * state%u%p
    call transform_p_to_s(work_p,work_s)
    out = out - work_s

    call deriv_y(grid,field%s,work_s)
    call transform_s_to_p(work_s,work_p)
    work_p = work_p * state%v%p
    call transform_p_to_s(work_p,work_s)
    out = out - work_s

    call deriv_z(grid,field%s,work_s)
    call transform_s_to_p(work_s,work_p)
    work_p = work_p * state%w%p
    call transform_p_to_s(work_p,work_s)
    out = out - work_s
    
    deallocate(work_p,work_s)
  end subroutine

  subroutine diffusion_rhs(grid,params,state,factor,field,out,diff,denom)
    !/********************************************************************|
    !| DIFFUSION_RHS                                                      |
    !| -------------                                                      |
    !| Calculate the diffusion term in the right hand side of an equation.|
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   field: the field for which to calculate the advection            |
    !|   out: an array used to contain the rhs information                |
    !|   diff: the diffusion coefficient                                  |
    !\********************************************************************/
    use grid_module, only: allocate_p,allocate_s
    use grid_ops_module, only: deriv_x,deriv_y,deriv_z
    use transform_module, only: transform_p_to_s,transform_s_to_p
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    type(field_type), intent(in) :: factor
    type(field_type), intent(in) :: field
    complex(kind=dp), intent(inout) :: out(:,:,:)
    real(kind=dp), intent(in) :: diff
    real(kind=dp), intent(in), optional :: denom(:,:,:)

    real(kind=dp), pointer :: work_p(:,:,:)
    complex (kind=dp), pointer :: work_s(:,:,:)
    complex (kind=dp), pointer :: work_rhs(:,:,:)

    real(kind=dp) :: rtmp,ksquared,hkx,hky,hkz
    integer :: i,j,k

    complex(kind=dp), parameter :: iu = (0._dp,1._dp)

    call allocate_p(grid,work_p)
    call allocate_s(grid,work_s)
    call allocate_s(grid,work_rhs)

    ! Calculate x-derivatives
    call deriv_x(grid,field%s,work_s)

    call transform_s_to_p(work_s,work_p)
    work_p = work_p * diff * factor%p
    call transform_p_to_s(work_p,work_s)

    call deriv_x(grid,work_s,work_rhs)

    ! Calculate y-derivatives
    call deriv_y(grid, field%s, work_s)

    call transform_s_to_p(work_s,work_p)
    work_p = work_p * diff * factor%p
    call transform_p_to_s(work_p,work_s)

    call deriv_y(grid,work_s,work_s)

    work_rhs = work_rhs + work_s

    ! Calculate z-derivatives
    call deriv_z(grid,field%s,work_s)

    call transform_s_to_p(work_s,work_p)
    work_p = work_p * diff * factor%p
    call transform_p_to_s(work_p,work_s)

    call deriv_z(grid,work_s,work_s)

    work_rhs = work_rhs + work_s

    ! Optionally divide by denominator
    if (present(denom)) then
      call transform_s_to_p(work_rhs,work_p)
      work_p = work_p / denom
      call transform_p_to_s(work_p,work_rhs)
    end if

    out = out + work_rhs

    deallocate(work_p,work_s,work_rhs)
  end subroutine

  subroutine viscous_stress_rhs(grid,params,state,field,out)
    !/********************************************************************|
    !| VISCOUS_STRESS_RHS                                                 |
    !| ------------------                                                 |
    !| Calculate the viscous stress term in the right hand side of the    |
    !| momentum equation.                                                 |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   field: the velocity field                                        |
    !|   out: an array used to contain the rhs information                |
    !\********************************************************************/
    use grid_module, only: allocate_p,allocate_s
    use transform_module, only: transform_p_to_s,transform_s_to_p
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    type(vector_type), intent(in) :: field
    complex(kind=dp), intent(inout) :: out(:,:,:,:)

    real(kind=dp), pointer :: work_p(:,:,:)
    complex (kind=dp), pointer :: work_s(:,:,:)
    complex (kind=dp), pointer :: div(:,:,:)

    real(kind=dp) :: rtmp,ksquared,hkx,hky,hkz
    integer :: i,j,k

    complex(kind=dp), parameter :: iu = (0._dp,1._dp)

    call allocate_p(grid,work_p)
    call allocate_s(grid,work_s)
    call allocate_s(grid,div)

    div = 0.0_dp
    call div_rhs(grid,params,state,state%vel,1.0_dp,div)

    work_s = 0.0_dp
    do j = 1,grid%y_s_count
      hky = grid%lky(j)
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        work_s(:,i,j) = work_s(:,i,j) + iu * hkx * div(:,i,j)
        do k = 1,2 * grid%z_modes
          hkz = grid%lkz(k)
          ksquared = hkx**2 + hky**2 + hkz**2
          work_s(k,i,j) = work_s(k,i,j) - ksquared * field%s(k,i,j,1)
        end do
      end do
    end do

    out(:,:,:,1) = out(:,:,:,1) + params%rho_ref * params%diff_vel * work_s

    work_s = 0.0_dp
    do j = 1,grid%y_s_count
      hky = grid%lky(j)
      work_s(:,:,j) = work_s(:,:,j) + iu * hky * div(:,:,j)
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        do k = 1,2 * grid%z_modes
          hkz = grid%lkz(k)
          ksquared = hkx**2 + hky**2 + hkz**2
          work_s(k,i,j) = work_s(k,i,j) - ksquared * field%s(k,i,j,2)
        end do
      end do
    end do

    out(:,:,:,2) = out(:,:,:,2) + params%rho_ref * params%diff_vel * work_s

    work_s = 0.0_dp
    do j = 1,grid%y_s_count
      hky = grid%lky(j)
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        do k = 1,2 * grid%z_modes
          hkz = grid%lkz(k)
          ksquared = hkx**2 + hky**2 + hkz**2
          work_s(k,i,j) = work_s(k,i,j) + iu * hkz * div(k,i,j)
          work_s(k,i,j) = work_s(k,i,j) - ksquared * field%s(k,i,j,3)
        end do
      end do
    end do

    out(:,:,:,3) = out(:,:,:,3) + params%rho_ref * params%diff_vel * work_s

    deallocate(work_p,work_s,div)
  end subroutine

  subroutine epsilon_rhs(grid,params,state,field,out,diff,denom)
    !/********************************************************************|
    !| EPSILON_RHS                                                        |
    !| -----------                                                        |
    !| Calculate the dissipative heating term in the right hand side of   |
    !| the energy equation. If present, divide through by denom.          |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   field: the field for which to calculate the advection            |
    !|   out: an array used to contain the rhs information                |
    !|   denom: the field to divide by at the end                         |
    !\********************************************************************/
    use grid_module, only: allocate_p,allocate_s
    use grid_ops_module, only: deriv_x,deriv_y,deriv_z
    use transform_module, only: transform_p_to_s,transform_s_to_p
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    type(field_type), intent(in), optional :: denom
    type(vector_type), intent(in) :: field
    complex(kind=dp), intent(inout) :: out(:,:,:)
    real(kind=dp), intent(in) :: diff

    real(kind=dp), pointer :: work_p(:,:,:)
    complex (kind=dp), pointer :: work_s(:,:,:)
    complex (kind=dp), pointer :: work_div(:,:,:)

    real(kind=dp) :: rtmp,ksquared,hkx,hky,hkz
    integer :: i,j,k,l

    complex(kind=dp), parameter :: iu = (0._dp,1._dp)

    call allocate_p(grid,work_p)
    call allocate_s(grid,work_s)
    call allocate_s(grid,work_div)

    do l = 1,3
      call deriv_x(grid,field%s(:,:,:,l),work_s)

      call transform_s_to_p(work_s,work_p)
      if (present(denom)) then
        work_p = diff * work_p * work_p / denom%p
      else
        work_p = diff * work_p * work_p
      end if
      call transform_p_to_s(work_p,work_s)

      out = out + work_s

      call deriv_y(grid,field%s(:,:,:,l),work_s)

      call transform_s_to_p(work_s,work_p)
      if (present(denom)) then
        work_p = diff * work_p * work_p / denom%p
      else
        work_p = diff * work_p * work_p
      end if
      call transform_p_to_s(work_p,work_s)

      out = out + work_s

      call deriv_z(grid,field%s(:,:,:,l),work_s)

      call transform_s_to_p(work_s,work_p)
      if (present(denom)) then
        work_p = diff * work_p * work_p / denom%p
      else
        work_p = diff * work_p * work_p
      end if
      call transform_p_to_s(work_p,work_s)

      out = out + work_s
    end do
    
    call deriv_x(grid,field%s(:,:,:,1),work_div)
    call deriv_y(grid,field%s(:,:,:,2),work_s)
    work_div = work_div + work_s
    call deriv_z(grid,field%s(:,:,:,3),work_s)
    work_div = work_div + work_s

    call transform_s_to_p(work_div,work_p)
    if (present(denom)) then
      work_p = diff * work_p * work_p / denom%p
    else
      work_p = diff * work_p * work_p
    end if
    call transform_p_to_s(work_p,work_s)

    out = out + work_s

    deallocate(work_p,work_s,work_div)
  end subroutine

  subroutine gradient_rhs(grid,params,state,field,factor,out)
    !/********************************************************************|
    !| GRADIENT_RHS                                                       |
    !| ------------                                                       |
    !| Calculate the viscous stress term in the right hand side of the    |
    !| momentum equation.                                                 |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   field: the field for which to calculate the advection            |
    !|   factor: a scalar factor multiplied on the term                   |
    !|   out: an array used to contain the rhs information                |
    !\********************************************************************/
    use grid_module, only: allocate_p,allocate_s
    use transform_module, only: transform_p_to_s,transform_s_to_p
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    type(field_type), intent(in) :: field
    real(kind=dp), intent(in) :: factor
    complex(kind=dp), intent(inout) :: out(:,:,:,:)

    real(kind=dp) :: hkx,hky,hkz
    integer :: i,j,k

    complex(kind=dp), parameter :: iu = (0._dp,1._dp)

    do j = 1,grid%y_s_count
      hky = grid%lky(j)
      out(:,:,j,2) = out(:,:,j,2) + factor * iu * hky * field%s(:,:,j)
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        out(:,i,j,1) = out(:,i,j,1) + factor * iu * hkx * field%s(:,i,j)
        do k = 1,2 * grid%z_modes
          hkz = grid%lkz(k)
          out(k,i,j,3) = out(k,i,j,3) + factor * iu * hkz * field%s(k,i,j)
        end do
      end do
    end do
  end subroutine

  subroutine div_rhs(grid,params,state,field,factor,out)
    !/********************************************************************|
    !| DIV_RHS                                                            |
    !| -------                                                            |
    !| Calculate the viscous stress term in the right hand side of the    |
    !| momentum equation.                                                 |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   field: the field for which to calculate the advection            |
    !|   factor: a scalar factor multiplied on the term                   |
    !|   out: an array used to contain the rhs information                |
    !\********************************************************************/
    use grid_module, only: allocate_p,allocate_s
    use transform_module, only: transform_p_to_s,transform_s_to_p
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state
    type(vector_type), intent(in) :: field
    real(kind=dp), intent(in) :: factor
    complex(kind=dp), intent(inout) :: out(:,:,:)

    real(kind=dp) :: hkx,hky,hkz
    integer :: i,j,k

    complex(kind=dp), parameter :: iu = (0._dp,1._dp)

    do j = 1,grid%y_s_count
      hky = grid%lky(j)
      out(:,:,j) = out(:,:,j) + factor * iu * hky * field%s(:,:,j,2)
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        out(:,i,j) = out(:,i,j) + factor * iu * hkx * field%s(:,i,j,1)
        do k = 1,2 * grid%z_modes
          hkz = grid%lkz(k)
          out(k,i,j) = out(k,i,j) + factor * iu * hkz * field%s(k,i,j,3)
        end do
      end do
    end do
  end subroutine

  subroutine implicit_update(grid,params,state,a0)
    !/********************************************************************|
    !| IMPLICIT_UPDATE                                                    |
    !| ---------------                                                    |
    !| Once the right hand side has been calculated and the effects of    |
    !| the prior steps have been included in state, use this to perform   |
    !| the final update to the new time with                              |
    !|     a0 * u(n+1) + dt * I(n+1) = b_rhs * RHS + PRIOR_TERMS          |
    !| This subroutine updates every variable in state.                   |
    !|   grid: a grid_type object that describes the grid                 |
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !|   a0: the multiplier on the new value in the final scheme          |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(inout) :: state
    real(kind=dp), optional :: a0

    real(kind=dp) :: a0in,ksquared,hkx,hky,hkz,rtmp,nu
    complex(kind=dp) :: tmp,ctmp,crl,div
    integer :: i,j,k

    complex(kind=dp), parameter :: iu = (0._dp,1._dp)

    a0in = 1.0_dp
    if (present(a0)) a0in = a0

    ! In the absence of an implicit rhs, scale the result by a0
    state%all_s = state%all_s / a0in

    if (.not. params%implicit) then
      return
    end if

    ! We are using the implicit rhs, so reverse the previous scaling
    ! Note: this wastes a small amount of time for readability
    state%rho%s = state%rho%s * a0in
    state%rhou%s = state%rhou%s * a0in
    state%rhov%s = state%rhov%s * a0in
    state%rhow%s = state%rhow%s * a0in

    ! An artificial viscosity may be used to relax initial conditions 
    ! to a hydrostatic state. This viscosity may be specificed to decay 
    ! over time.
    nu = params%visc_stab
    if (params%steps_stab > 0) then
      nu = nu * exp(-real(state%timestep) / real(params%steps_stab))
    end if

    ! TODO: Implicit sound waves are handled only in 2D but shound be 
    ! generalized to 3D
    do j = 1,grid%y_s_count
      hky = grid%lky(j)
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)

        if (hkx == 0.0_dp) then
          ! Handle horizontally invariant modes
          do k = 1,2 * grid%z_modes
            hkz = grid%lkz(k)
            ksquared = hkx**2 + hky**2 + hkz**2

            rtmp = a0in + state%dt * ksquared * params%visc_stab
            ctmp = iu * params%cs**2 * state%dt

            ! For horizontally invariant modes, the u solution is trivial
            state%rhou%s(k,i,j) = state%rhou%s(k,i,j) / rtmp

            rtmp = rtmp + state%dt**2 * ksquared * params%cs**2

            ! Solve for the updated w field
            state%rhow%s(k,i,j) = (state%rhow%s(k,1,j) &
            - ctmp * hkz * state%rho%s(k,i,j) / a0in) / rtmp
            
            ! Solve for the updated density field, using the updated vel
            state%rho%s(k,i,j) = (state%rho%s(k,i,j) &
            - iu * hkz * state%dt * state%rhow%s(k,i,j)) / a0in
          end do
        else
          ! Handle horizontally varying modes
          do k = 1,2 * grid%z_modes
            hkz = grid%lkz(k)
            ksquared = hkx**2 + hky**2 + hkz**2

            rtmp = a0in + state%dt**2 * ksquared * params%cs**2 &
                  / (a0in + state%dt * params%visc_stab * ksquared)
    
            ! Calculate the updated divergence field
            tmp = hkx * state%rhou%s(k,i,j) + hkz * state%rhow%s(k,i,j)
            div = iu * tmp / (a0in + state%dt * ksquared * params%visc_stab)

            ! Using the updated divergence field, solve for the rho field
            state%rho%s(k,i,j) = (state%rho%s(k,i,j) - state%dt * div) / rtmp

            rtmp = a0in + state%dt * params%visc_stab * ksquared
            ctmp = iu * params%cs**2 * state%dt

            ! Using the updated rho field, solve for momentum
            state%rhou%s(k,i,j) = (state%rhou%s(k,i,j) &
            - ctmp * hkx * state%rho%s(k,i,j)) / rtmp
            state%rhow%s(k,i,j) = (state%rhow%s(k,i,j) &
            - ctmp * hkz * state%rho%s(k,i,j)) / rtmp
          end do
        end if
      end do
    end do
  end subroutine

  subroutine diffuse(grid,state,in,out,diff,a0)
    !/********************************************************************|
    !| DIFFUSE                                                            |
    !| -------                                                            |
    !| Once the right hand side has been calculated and the effects of    |
    !| the prior steps have been included in state, use this to perform   |
    !| the final update to the new time with                              |
    !|     a0 * u(n+1) + dt * I(n+1) = b_rhs * RHS + PRIOR_TERMS          |
    !| This subroutine updates one variable.                              |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   state: a state_type object containing the simulation data        |
    !|   in: the input array containing the prior time information        |
    !|   out: the output array to contain the new variable state          |
    !|   diff: the diffusion coefficient                                  |
    !|   a0: the multiplier on the new value in the final scheme          |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    type(state_type), intent(in) :: state
    complex(kind=dp), intent(in) :: in(:,:,:)
    complex(kind=dp), intent(out) :: out(:,:,:)
    real(kind=dp), intent(in) :: diff
    real(kind=dp), intent(in) :: a0

    real(kind=dp) :: rtmp,ksquared,hkx,hky,hkz
    integer :: i,j,k

    do j = 1,grid%y_s_count
      hky = grid%lky(j)
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        do k = 1,2 * grid%z_modes
          hkz = grid%lkz(k)
          ksquared = hkx**2 + hky**2 + hkz**2
          rtmp = a0 + state%dt * ksquared * diff
          out(k,i,j) = in(k,i,j) / rtmp
        end do
      end do
    end do
  end subroutine

  subroutine update_eos(grid,params,state)
    !/********************************************************************|
    !| UPDATE_EOS                                                         |
    !| ----------                                                         |
    !| Update the entropy, pressure, specific heat capacity, and chemical |
    !| potential of the system in spectral space. Also reevaluate the     |
    !| velocities. Update the physical representations of those fields.   |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !\********************************************************************/
    use params_module, only: eos_ideal,eos_linear
    use transform_module, only: transform_p_to_s
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state

    select case (params%eosin)
    case (eos_linear)
      call update_eos_linear(grid,params,state)
    case (eos_ideal)
      call update_eos_ideal(grid,params,state)
    end select

    state%u%p = state%rhou%p / state%rho%p
    state%v%p = state%rhov%p / state%rho%p
    state%w%p = state%rhow%p / state%rho%p

    call transform_p_to_s(state%u%p,state%u%s)
    call transform_p_to_s(state%v%p,state%v%s)
    call transform_p_to_s(state%w%p,state%w%s)

    call transform_p_to_s(state%p%p,state%p%s)
    call transform_p_to_s(state%t%p,state%t%s)
    call transform_p_to_s(state%cp%p,state%cp%s)
    call transform_p_to_s(state%mu%p,state%mu%s)
  end subroutine

  subroutine update_eos_linear(grid,params,state)
    !/********************************************************************|
    !| UPDATE_EOS_LINEAR                                                  |
    !| -----------------                                                  |
    !| Update the entropy, pressure, specific heat capacity, and chemical |
    !| potential of the system in spectral space using the linear         |
    !| equation of state from Dewar (PC).                                 |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state

    real(kind=dp) :: ratio,en_ref
    integer :: i,j,k

    ratio = params%beta_c / params%beta_t
    en_ref = (params%beta_t * params%t_ref / params%rho_ref + params%p_ref / params%cs**2 / params%rho_ref**2) / params%alr

    state%p%p = 1.0_dp / params%rho_ref - 1.0_dp / state%rho%p + params%alr * state%en%p
    state%p%p = state%p%p * params%rho_ref**2 * params%cs**2

    state%t%p = params%t_ref + params%rho_ref * params%alr / params%beta_t * state%en%p
    state%t%p = state%t%p + params%alr * state%p%p + ratio * state%c%p
    state%t%p = state%t%p + params%p_ref / params%cs**2 / params%rho_ref / params%beta_t

    state%cp%p = state%rho%p * params%beta_t * state%t%p / params%alr / params%rho_ref

    state%mu%p = ratio * (state%en%p + en_ref)
  end subroutine

  subroutine update_eos_ideal(grid,params,state)
    !/********************************************************************|
    !| UPDATE_EOS_IDEAL                                                   |
    !| ----------------                                                   |
    !| Update the entropy, pressure, specific heat capacity, and chemical |
    !| potential of the system in spectral space using the ideal gas      |
    !| equation of state.                                                 |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   state: a state_type object containing the simulation data        |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(state_type), intent(in) :: state

    real(kind=dp) :: gam_m1,cp,rtmp
    integer :: i,j,k

    gam_m1 = params%gamma - 1.0_dp
    cp = params%cv * params%gamma
    rtmp = params%cv * gam_m1 * params%t_ref

    state%p%p = (rtmp * state%rho%p * exp(state%en%p / cp))**params%gamma &
                * params%p_ref**(-gam_m1)

    state%t%p = params%t_ref * exp(state%en%p / cp) & 
                * (state%p%p / params%p_ref)**(gam_m1 / params%gamma)

    state%cp%p = state%rho%p * cp

    state%mu%p = 0.0_dp
  end subroutine

end module equations_module
