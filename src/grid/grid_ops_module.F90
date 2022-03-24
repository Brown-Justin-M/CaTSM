module grid_ops_module
  !/**********************************************************************|
  !| GRID OPS MODULE                                                      |
  !| ---------------                                                      |
  !| This module contains several functions and subroutines designed to   |
  !| operate with the grid_type, including several basic vector calculus  |
  !| operations.                                                          |
  !\**********************************************************************/
  use grid_module, only: dp,li,grid_type
  implicit none
  private
  public div,curl,deriv_x,deriv_y,deriv_z
  public cross,remove_div

contains
  !/**********************************************************************|
  !| FUNCTIONS                                                            |
  !\**********************************************************************/
  function div(grid,time,x)
    !/********************************************************************|
    !| DIV                                                                |
    !| ---                                                                |
    !| Calculates the divergence of the spectral vector x.                |
    !|   grid: a grid_type object that describes the grid                 |
    !|   time: the current time                                           |
    !|   x: a spectral vector array                                       |
    !|   returns: a 3D array of the divergence of x in spectral space     |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    real(kind=dp), intent(in) :: time
    complex(kind=dp), intent(in) :: x(:,:,:,:)
    complex(kind=dp) :: div(grid%z_s_count,grid%x_s_count,grid%y_s_count)

    complex(kind=dp), parameter :: iu = (0.0_dp,1.0_dp) 
    real(kind=dp) :: hkx,hky,hkz
    integer :: i,j,k

    do j = 1,grid%y_s_count
      hky = grid%lky(j)
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        do k = 1,grid%z_s_count
          hkz = grid%lkz(k)
          div(k,i,j) = iu * hkx * x(k,i,j,1)
          div(k,i,j) = div(k,i,j) + iu * hky * x(k,i,j,2)
          div(k,i,j) = div(k,i,j) + iu * hkz * x(k,i,j,3)
        end do
      end do
    end do
  end function

  function get_shear(grid,time)
    !/********************************************************************|
    !| GET_SHEAR                                                          |
    !| ---------                                                          |
    !| Calculate the values du/dz and dv/dz of the background shear at    |
    !| time.                                                              |
    !|   grid: a grid_type object that describes the grid                 |
    !|   time: the current time                                           |
    !|   returns: an array containing (du/dz,dv/dz)                       |
    !\********************************************************************/
    implicit none

    type(grid_type) :: grid
    real(kind=dp), intent(in) :: time
    real(kind=dp) :: get_shear(2)

    real(kind=dp) :: rtmp

    rtmp = time + grid%shear_time_offset

    get_shear(1) = grid%x_shear + grid%x_shear_ramp * rtmp
    get_shear(2) = grid%y_shear + grid%y_shear_ramp * rtmp

    rtmp = grid%shear_freq * (time + grid%shear_time_offset)

    get_shear(1) = get_shear(1) * cos(rtmp)

    rtmp = grid%shear_freq * (time + grid%shear_time_offset) + grid%y_shear_phase

    get_shear(2) = get_shear(2) * cos(rtmp)
  end function

  !/**********************************************************************|
  !| SUBROUTINES                                                          |
  !\**********************************************************************/ 
  subroutine curl(grid,time,in,out)
    !/********************************************************************|
    !| CURL                                                               |
    !| ----                                                               |
    !| Calculates the curl of the spectral vector in.                     |
    !|   grid: a grid_type object that describes the grid                 |
    !|   time: the current time                                           |
    !|   in: a spectral vector array                                      |
    !|   out: a 4D array of the curl of in in spectral space              |
    !\********************************************************************/
    implicit none 

    type(grid_type), intent(in) :: grid
    real(kind=dp), intent(in) :: time
    complex(kind=dp), intent(in) :: in(:,:,:,:)
    complex(kind=dp), intent(out) :: out(:,:,:,:)

    complex(kind=dp), parameter :: iu = (0.0_dp,1.0_dp) 
    real(kind=dp) :: hkx,hky,hkz
    integer :: i,j,k

    do j = 1,grid%y_s_count
      hky = grid%lky(j)
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        do k = 1,grid%z_s_count
          hkz = grid%lkz(k)
          out(k,i,j,1) = iu * (hky * in(k,i,j,3) - hkz * in(k,i,j,2))
          out(k,i,j,2) = iu * (hkz * in(k,i,j,1) - hkx * in(k,i,j,3))
          out(k,i,j,3) = iu * (hkx * in(k,i,j,2) - hky * in(k,i,j,1))
        end do
      end do
    end do
  end subroutine

  subroutine deriv_x(grid,in,out)
    !/********************************************************************|
    !| DERIV_X                                                            |
    !| -------                                                            |
    !| Calculates the x derivative of the spectral scalar in.             |
    !|   grid: a grid_type object that describes the grid                 |
    !|   time: the current time                                           |
    !|   in: a spectral scalar array                                      |
    !|   out: a 3D spectral array of d(in)/dx                             |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    complex(kind=dp), intent(in) :: in(:,:,:)
    complex(kind=dp), intent(out) :: out(:,:,:)

    complex(kind=dp), parameter :: iu = (0.0_dp,1.0_dp) 
    real(kind=dp) :: hkx
    integer :: i,j

    do j = 1,grid%y_s_count
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        out(:,i,j) = iu * hkx * in(:,i,j)
      end do
    end do
  end subroutine

  subroutine deriv_y(grid,in,out)
    !/********************************************************************|
    !| DERIV_Y                                                            |
    !| -------                                                            |
    !| Calculates the y derivative of the spectral scalar in.             |
    !|   grid: a grid_type object that describes the grid                 |
    !|   time: the current time                                           |
    !|   in: a spectral scalar array                                      |
    !|   out: a 3D spectral array of d(in)/dy                             |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    complex(kind=dp), intent(in) :: in(:,:,:)
    complex(kind=dp), intent(out) :: out(:,:,:)

    complex(kind=dp), parameter :: iu = (0.0_dp,1.0_dp) 
    real(kind=dp) :: hky
    integer :: j

    do j = 1,grid%y_s_count
      hky = grid%lky(j)
      out(:,:,j) = iu * hky * in(:,:,j)
    end do
  end subroutine

  subroutine deriv_z(grid,in,out)
    !/********************************************************************|
    !| DERIV_Z                                                            |
    !| -------                                                            |
    !| Calculates the z derivative of the spectral scalar in.             |
    !|   grid: a grid_type object that describes the grid                 |
    !|   time: the current time                                           |
    !|   in: a spectral scalar array                                      |
    !|   out: a 3D spectral array of d(in)/dz                             |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    complex(kind=dp), intent(in) :: in(:,:,:)
    complex(kind=dp), intent(out) :: out(:,:,:)

    complex(kind=dp), parameter :: iu = (0.0_dp,1.0_dp) 
    real(kind=dp) :: hkx,hky,hkz
    integer :: i,j,k

    do j = 1,grid%y_s_count
      hky = grid%lky(j)
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        do k = 1,grid%z_s_count
          hkz = grid%lkz(k)
          out(k,i,j) = iu * hkz * in(k,i,j)
        end do
      end do
    end do
  end subroutine

  subroutine cross(a,b,out,dim)
    !/********************************************************************|
    !| CROSS                                                              |
    !| -----                                                              |
    !| Calculates one component of the cross product of a with b.         |
    !|   a: a spectral scalar array                                       |
    !|   b: a spectral scalar array                                       |
    !|   dim: the desired component of the cross product (1 for x, etc.)  |
    !|   out: a 3D spectral array of the component of the cross product   |
    !\********************************************************************/
    implicit none

    real(kind=dp), pointer, intent(in) :: a(:,:,:,:)
    real(kind=dp), pointer, intent(in) :: b(:,:,:,:)
    real(kind=dp), pointer, intent(out) :: out(:,:,:)
    integer, intent(in) :: dim

    select case (dim)
    case (1)
      out = a(:,:,:,2) * b(:,:,:,3) - a(:,:,:,3) * b(:,:,:,2)
    case (2)
      out = a(:,:,:,3) * b(:,:,:,1) - a(:,:,:,1) * b(:,:,:,3)
    case (3)
      out = a(:,:,:,1) * b(:,:,:,2) - a(:,:,:,2) * b(:,:,:,1)
    end select
  end subroutine

  subroutine remove_div(grid,time,in,x,factor)
    !/********************************************************************|
    !| REMOVE_DIV                                                         |
    !| ----------                                                         |
    !| Removes divergence by solving the following for x                  |
    !|       laplace(p) = div(in)                                         |
    !|     factor * x = factor * x - grad(p) - (in - factor * x)          |
    !| Often in incompressible systems, it is sufficient to solve the     |
    !| pressure such that the rhs of the system is divergence free. This  |
    !| is not the case if the coordinate system is changing with time. In |
    !| that case, the velocity update equation is expected to look like   |
    !|     u = in + factor * x,                                           |
    !| where "factor" is the multiplier on the rhs of the current time,   |
    !| "x" is the divergence-removed rhs, and "in" is the estimate of the |
    !| velocity (accounting for any prior steps in a multi-step method).  |
    !|   grid: a grid_type object that describes the grid                 |
    !|   time: the current time                                           |
    !|   in: a 4D array containing the estimate of the velocity           |
    !|   x: a 4D array for the quantity from which to remove divergence   |
    !|   factor: the factor on x, if omitted, defaults to 1               |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    real(kind=dp), intent(in) :: time
    complex(kind=dp), intent(in) :: in(:,:,:,:)
    complex(kind=dp), intent(inout) :: x(:,:,:,:)
    real(kind=dp), intent(in), optional :: factor

    real(kind=dp) :: factorin

    complex(kind=dp), parameter :: iu = (0.0_dp,1.0_dp) 
    real(kind=dp) :: hkx,hky,hkz,ksquared
    integer :: i,j,k
    complex(kind=dp) :: ctmp

    factorin = 1.0_dp
    if (present(factor)) factorin = factor

    do j = 1,grid%y_s_count
      hky = grid%lky(j)
      do i = 1,grid%x_s_count
        hkx = grid%lkx(i)
        do k = 1,grid%z_s_count
          hkz = grid%lkz(k)
          ksquared = max(hkx**2 + hky**2 + hkz**2,epsilon(1.0_dp))

          ctmp = hkx * in(k,i,j,1)
          ctmp = ctmp + hky * in(k,i,j,2)
          ctmp = ctmp + hkz * in(k,i,j,3)

          ctmp = ctmp / ksquared / factorin

          x(k,i,j,1) = x(k,i,j,1) - hkx * ctmp
          x(k,i,j,2) = x(k,i,j,2) - hky * ctmp
          x(k,i,j,3) = x(k,i,j,3) - hkz * ctmp
        end do
      end do
    end do
  end subroutine

end module grid_ops_module