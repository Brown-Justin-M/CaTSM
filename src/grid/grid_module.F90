module grid_module
  !/**********************************************************************|
  !| GRID MODULE                                                          |
  !| -----------                                                          |
  !| This module contains the grid_type derived type, which contains the  |
  !| major parameters and arrays governing the grid. The module also has  |
  !| several subroutines useful for allocating arrays and measuring any   |
  !| grid deformation caused by shear.                                    |
  !\**********************************************************************/
  use transform_module, only: dp,li
  implicit none
  private
  public dp,li,pi
  public grid_type
  public init_grid,free_grid
  public allocate_s,allocate_p,allocate_s_vector,allocate_p_vector

  !/**********************************************************************|
  !| PARAMETERS                                                           |
  !\**********************************************************************/
  real(kind=dp), parameter :: pi = acos(-1.0_dp)

  !/**********************************************************************|
  !| DERIVED TYPES                                                        |
  !\**********************************************************************/
  type grid_type
    !/********************************************************************|
    !| GRID_TYPE                                                          |
    !| ---------                                                          |
    !| This derived type contains the quantities and arrays that govern   |
    !| the general properties of the grid. It is expected that the grid   |
    !| will first be read from file or that init_grid is called with the  |
    !| the x_modes, y_modes, and z_modes arguments explicitly. To read    |
    !| from file, use                                                     |
    !|     read(unit,fmt="(DT)") grid                                     |
    !| The namelist file should contain a &grid_values namelist with      |
    !| x_modes, y_modes, z_modes, and any of the remaining variables      |
    !| below that aren't listed as "derived." Once init_grid is called,   |
    !| The values in the grid object should not be changed by the user.   |
    !| Notes that each "mode" includes the + and - wavenumber, so it      |
    !| takes 8 complex numbers to track 4 modes, for example.             |
    !\********************************************************************/
    ! The following are read in from file and determine the grid structure
    integer :: x_modes = 0 ! The number of Fourier modes in the x direction
    integer :: y_modes = 0 ! The number of Fourier modes in the y direction
    integer :: z_modes = 0 ! The number of Fourier modes in the z direction
    real(kind=dp) :: x_length = 1.0_dp ! The length of the domain in x
    real(kind=dp) :: y_length = 1.0_dp ! The length of the domain in y
    real(kind=dp) :: z_length = 1.0_dp ! The length of the domain in z
    integer :: y_pencils = 1 ! The number of pencils in physical space in y
    integer :: z_pencils = 1 ! The number of pencils in physical space in z
    integer :: dealias_factor = 3 ! The factor of points used to dealias

    ! The following are read in from file and determine the grid motion
    real(kind=dp) :: x_shear = 0.0_dp ! The vertical shear in u, du/dz
    real(kind=dp) :: y_shear = 0.0_dp ! The vertical shear in v, du/dz
    real(kind=dp) :: shear_freq = 0.0_dp ! The frequency of the shear (cycles)
    real(kind=dp) :: shear_time_offset = 0.0_dp ! The offset of shear time
    real(kind=dp) :: x_shear_ramp = 0.0_dp ! The rate of change of x_shear
    real(kind=dp) :: y_shear_ramp = 0.0_dp ! The rate of change of y_shear
    real(kind=dp) :: y_shear_phase = 0.0_dp ! The starting phase of y_shear

    ! All following are derived from the above and should not be changed
    ! The derived properties of the grid structure
    integer :: x_n ! The number of physical grid points in the x direction
    integer :: y_n ! The number of physical grid points in the y direction
    integer :: z_n ! The number of physical grid points in the z direction
    real(kind=dp) :: dx,dy,dz ! The grid spacing in x, y, and z
    real(kind=dp) :: dv ! The volume of a single cell
    real(kind=dp) :: area ! The total horizontal area of the domain
    real(kind=dp) :: volume ! The total volume of the domain
    real(kind=dp) :: xy_n ! The number of grid points in a horizontal slice
    real(kind=dp) :: xyz_n ! The number of grid points in the simulation

    real(kind=dp), pointer :: kx(:),ky(:),kz(:) ! global x,y,z wavenumbers
    real(kind=dp), pointer :: x(:),y(:),z(:) ! global x,y,z values
    real(kind=dp), pointer :: lkx(:),lky(:),lkz(:) ! local x,y,z wavenumbers
    real(kind=dp), pointer :: lx(:),ly(:),lz(:) ! local x,y,z values
    real(kind=dp), pointer :: g(:) ! global normalized gravitational acceleration
    real(kind=dp), pointer :: lg(:) ! local normalized gravitational acceleration

    ! The following are used internally to track the offset between local
    ! and global quantities during parallel simulations in spectral and 
    ! physical space
    integer :: x_s_offset=0,y_s_offset=0,z_s_offset=0
    integer :: x_p_offset=0,y_p_offset=0,z_p_offset=0
    ! The following are used internally to track the total number of modes
    ! or cells in spectral or physical space
    integer :: x_s_count,y_s_count,z_s_count
    integer :: x_p_count,y_p_count,z_p_count

    contains
      private
      procedure, pass :: read => read_grid
      generic, public :: read(formatted) => read
  end type grid_type

contains
  !/**********************************************************************|
  !| SUBROUTINES                                                          |
  !\**********************************************************************/ 
  subroutine read_grid(gridin,unit,iotype,v_list,iostat,iomsg)
    !/********************************************************************|
    !| READ_GRID                                                          |
    !| ---------                                                          |
    !| Read the grid in from a namelist called "grid." To use this, run   |
    !|     read(unit,fmt="(DT)") grid                                     |
    !|   gridin: the grid object to read into                             |
    !|   unit: the file unit to be read from                              |
    !\********************************************************************/
    implicit none

    class(grid_type), intent(inout) :: gridin
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg

    integer :: y_pencils,z_pencils
    integer :: x_modes,y_modes,z_modes
    real(kind=dp) :: x_length,y_length,z_length
    real(kind=dp) :: x_shear,shear_freq,x_shear_ramp
    real(kind=dp) :: y_shear,y_shear_phase,y_shear_ramp
    real(kind=dp) :: shear_time_offset
    integer :: dealias_factor

    namelist /grid/ y_pencils,z_pencils
    namelist /grid/ x_modes,y_modes,z_modes
    namelist /grid/ x_length,y_length,z_length
    namelist /grid/ x_shear,shear_freq,x_shear_ramp
    namelist /grid/ y_shear,y_shear_ramp,y_shear_phase
    namelist /grid/ shear_time_offset
    namelist /grid/ dealias_factor

    ! Load any previous (default) values from the grid object
    x_modes = gridin%x_modes
    y_modes = gridin%y_modes
    z_modes = gridin%z_modes
    y_pencils = gridin%y_pencils
    z_pencils = gridin%z_pencils
    x_length = gridin%x_length
    y_length = gridin%y_length
    z_length = gridin%z_length
    dealias_factor = gridin%dealias_factor

    x_shear = gridin%x_shear
    y_shear = gridin%y_shear
    shear_freq = gridin%shear_freq
    shear_time_offset = gridin%shear_time_offset
    y_shear_phase = gridin%y_shear_phase
    x_shear_ramp = gridin%x_shear_ramp
    y_shear_ramp = gridin%y_shear_ramp

    read(unit,nml=grid,iostat=iostat) 

    ! Copy the results from the file
    gridin%x_modes = x_modes
    gridin%y_modes = y_modes
    gridin%z_modes = z_modes
    gridin%x_length = x_length
    gridin%y_length = y_length
    gridin%z_length = z_length
    gridin%y_pencils = y_pencils
    gridin%z_pencils = z_pencils
    
    gridin%dealias_factor = dealias_factor

    gridin%x_shear = x_shear
    gridin%y_shear = y_shear
    gridin%shear_freq = shear_freq
    gridin%shear_time_offset = shear_time_offset
    gridin%x_shear_ramp = x_shear_ramp
    gridin%y_shear_phase = y_shear_phase
    gridin%y_shear_ramp = y_shear_ramp
  end subroutine

  subroutine init_grid(grid,x_modes,y_modes,z_modes)
    !/********************************************************************|
    !| INIT_GRID                                                          |
    !| ---------                                                          |
    !| This subroutine allocates the arrays contained within a grid       |
    !| object. This includes the wave vector arrays and x, y, z for       |
    !| calculations in physical space. It also initializes the values of  |
    !| these arrays and other calculated grid quantities. This should be  |
    !| called after read_grid if the user wants to specify any custom     |
    !| grid properties like shear. If the modes arguments are omitted,    |
    !| the previous values in the grid object are kept.                   |
    !|   grid: a grid_type object that describes the grid                 |
    !|   x_modes: the number of modes in the x direction                  |
    !|   y_modes: the number of modes in the y direction                  |
    !|   z_modes: the number of modes in the z direction                  |
    !\********************************************************************/
    use comm_module, only: get_tasks,abort_comm
    use transform_module, only: init_plans
    implicit none

    type(grid_type), intent(inout) :: grid
    real(kind=dp), intent(in), optional :: x_modes,y_modes,z_modes

    if (present(x_modes)) grid%x_modes = x_modes
    if (present(y_modes)) grid%y_modes = y_modes
    if (present(z_modes)) grid%z_modes = z_modes

    if (grid%x_modes == 0) then
      print*,"FATAL: grid%x_modes,y_modes,z_modes should be set before"
      print*,"FATAL: init_grid callor they should be manually specified"
      print*,"FATAL: as arguments to init_grid.(The code does not support"
      print*,"FATAL: 0 modes in x.)"
      call abort_comm
    end if

    grid%x_n = max(grid%dealias_factor * (2 * grid%x_modes) / 2,1) 
    grid%y_n = max(grid%dealias_factor * (2 * grid%y_modes) / 2,1) 
    grid%z_n = max(grid%dealias_factor * (2 * grid%z_modes) / 2,1) 

    grid%dx = grid%x_length / grid%x_n
    grid%dy = grid%y_length / grid%y_n
    grid%dz = grid%z_length / grid%z_n
    grid%dv = grid%dx * grid%dy * grid%dz

    grid%area = grid%x_length * grid%y_length
    grid%volume = grid%x_length * grid%y_length * grid%z_length
    grid%xy_n = grid%x_n * grid%y_n
    grid%xyz_n = grid%x_n * grid%y_n * grid%z_n

      if (grid%y_modes == 0) then 
      grid%y_pencils = 1
      grid%z_pencils = get_tasks()
    end if

    if (grid%y_pencils * grid%z_pencils /= get_tasks()) then
      print*,"FATAL: Incorrect number of processors,",get_tasks(), &
        & ", expected",grid%y_pencils * grid%z_pencils
      call abort_comm
    end if

    call init_plans(grid%x_modes,grid%y_modes,grid%z_modes, &
      & grid%x_n,grid%y_n,grid%z_n,grid%y_pencils,grid%z_pencils)

    ! Read the offsets and counts from the transforms for parallel runs
    grid%x_s_count = grid%x_modes + 1
    grid%y_s_count = max(2 * grid%y_modes,1)
    grid%z_s_count = max(2 * grid%z_modes,1)
    grid%x_p_count = grid%x_n
    grid%y_p_count = grid%y_n
    grid%z_p_count = grid%z_n

    call init_wavenumbers(grid)
    call init_xyz(grid)
  end subroutine

  subroutine init_wavenumbers(grid)
    !/********************************************************************|
    !| INIT_WAVENUMBERS                                                   |
    !| ----------------                                                   |
    !| Initialize the wavenumber arrays within a grid object.             |
    !|   grid: a grid_type object that describes the grid                 |
    !\********************************************************************/
    implicit none
    type(grid_type), intent(inout) :: grid

    real(kind=dp) :: rtmp
    integer :: i

    allocate(grid%kx(0:grid%x_modes))
    allocate(grid%ky(0:max(2*grid%y_modes - 1,0))) 
    allocate(grid%kz(0:max(2*grid%z_modes - 1,0)))

    ! Initialize the kx array
    rtmp = 2.0_dp * pi / grid%x_length
    do i = 0,grid%x_modes
      grid%kx(i) = i * rtmp
    end do

    ! Initialize the ky array
    rtmp = 2.0_dp * pi / grid%y_length
    if (grid%y_modes == 0) then 
      grid%ky = 0.0_dp 
    else 
      do i = 0,grid%y_modes
        grid%ky(i) = i * rtmp
      end do
      do i = grid%y_modes + 1,2 * grid%y_modes - 1
        grid%ky(i) = -(2 * grid%y_modes - i) * rtmp
      end do
    end if

    ! Initialize the kz array
    rtmp = 2.0_dp * pi / grid%z_length
    if (grid%z_modes == 0) then 
      grid%kz = 0.0_dp 
    else
      do i = 0,grid%z_modes
        grid%kz(i) = i * rtmp
      end do
      do i = grid%z_modes + 1,2 * grid%z_modes - 1
        grid%kz(i) = -(2 * grid%z_modes - i) * rtmp
      end do
    end if

    ! Set up the local wavenumbers
    grid%lkx(1:) => grid%kx(grid%x_s_offset:)
    grid%lky(1:) => grid%ky(grid%y_s_offset:)
    grid%lkz(1:) => grid%kz(grid%z_s_offset:)
  end subroutine

  subroutine init_xyz(grid)
    !/********************************************************************|
    !| INIT_XYZ                                                           |
    !| --------                                                           |
    !| Initialize the x, y, and z arrays within a grid object.            |
    !|   grid: a grid_type object that describes the grid                 |
    !\********************************************************************/
    implicit none
    type(grid_type), intent(inout) :: grid

    integer :: i

    allocate(grid%x(grid%x_n))
    allocate(grid%y(grid%y_n))
    allocate(grid%z(grid%z_n))
    allocate(grid%g(grid%z_n))
    
    ! Initialize the x array
    do i = 1,grid%x_n
      grid%x(i) = (i - 1) * grid%dx
    end do

    ! Initialize the y array
    if (grid%y_modes == 0) then 
      grid%y = 0.0_dp 
    else 
      do i = 1,grid%y_n
        grid%y(i) = (i - 1) * grid%dy
      end do
    end if

    ! Initialize the z array
    if (grid%z_n == 0) then 
      grid%z = 0.0_dp 
    else
      do i = 1,grid%z_n
        grid%z(i) = (i - 1) * grid%dz
      end do
    end if
    grid%g = -2.0_dp*atan(sin(2.0_dp*pi*grid%z/grid%z_length)/0.01_dp)/pi

    ! Set up the local position arrays
    grid%lx(1:) => grid%x(grid%x_p_offset + 1:)
    grid%ly(1:) => grid%y(grid%y_p_offset + 1:)
    grid%lz(1:) => grid%z(grid%z_p_offset + 1:)

    grid%lg(1:) => grid%g(grid%z_p_offset + 1:)
  end subroutine
    
  subroutine  free_grid(grid)
    !/********************************************************************|
    !| FREE_GRID                                                          |
    !| ---------                                                          |
    !| This subroutine deallocates the arrays contained within a grid     |
    !| object.                                                            |
    !\********************************************************************/
    use transform_module, only: free_plans
    implicit none

    type(grid_type), intent(inout) :: grid

    call free_plans

    deallocate(grid%kx,grid%ky,grid%kz)
    deallocate(grid%x,grid%y,grid%z)
  end subroutine

  subroutine allocate_s(grid,cwork)
    !/********************************************************************|
    !| ALLOCATE_S                                                         |
    !| ----------                                                         |
    !| Allocate an array with the appropriate dimensions for a scalar in  |
    !| spectral coordinates.                                              |
    !|   grid: a grid_type object that describes the grid                 |
    !|   cwork: the 3D array to allocate                                  |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    complex(kind=dp), pointer, intent(out) :: cwork(:,:,:)

    allocate(cwork(2*grid%z_modes,grid%x_s_count,grid%y_s_count))
    cwork = 0.0_dp
  end subroutine

  subroutine allocate_s_vector(grid,cwork,dim)
    !/********************************************************************|
    !| ALLOCATE_S_VECTOR                                                  |
    !| -----------------                                                  |
    !| Allocate an array with the appropriate dimensions for a vector in  |
    !| spectral coordinates.                                              |
    !|   grid: a grid_type object that describes the grid                 |
    !|   cwork: the 4D array to allocate                                  |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    complex(kind=dp), pointer, intent(out) :: cwork(:,:,:,:)
    integer :: dim

    allocate(cwork(2*grid%z_modes,grid%x_s_count,grid%y_s_count,dim))
    cwork = 0.0_dp
  end subroutine

  subroutine allocate_p(grid,rwork)
    !/********************************************************************|
    !| ALLOCATE_P                                                         |
    !| ----------                                                         |
    !| Allocate an array with the appropriate dimensions for a scalar in  |
    !| physical coordinates.                                              |
    !|   grid: a grid_type object that describes the grid                 |
    !|   cwork: the 3D array to allocate                                  |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    real(kind=dp), pointer, intent(out) :: rwork(:,:,:)

    allocate(rwork(grid%x_n,grid%y_p_count,grid%z_p_count))
    rwork = 0.0_dp
  end subroutine

  subroutine allocate_p_vector(grid,rwork,dim)
    !/********************************************************************|
    !| ALLOCATE_P_VECTOR                                                  |
    !| -----------------                                                  |
    !| Allocate an array with the appropriate dimensions for a vector in  |
    !| physical coordinates.                                              |
    !|   grid: a grid_type object that describes the grid                 |
    !|   cwork: the 4D array to allocate                                  |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    real(kind=dp), pointer, intent(out) :: rwork(:,:,:,:)
    integer :: dim

    allocate(rwork(grid%x_n,grid%y_p_count,grid%z_p_count,dim))
    rwork = 0.0_dp
  end subroutine

end module grid_module