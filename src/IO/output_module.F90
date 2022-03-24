module output_module
  !/**********************************************************************|
  !| OUTPUT MODULE                                                        |
  !| -------------                                                        |
  !| This module contains the interface to write to an output file from a |
  !| state_type derived type. The output module supports spectral and     |
  !| physical output files in 3D in NetCDF format and vertical profiles   |
  !| in NetCDF format. It also supports a diagnostic file in CSV format.  |
  !\**********************************************************************/
  use grid_module, only: dp,li,grid_type
  use params_module, only: params_type
  use fields_module, only: state_type
  implicit none
  private
  public nc_ids_type,controllers_type
  public open_files,close_files,write_output_files
  public check

  !/**********************************************************************|
  !| FLAGS                                                                |
  !\**********************************************************************/
  integer, parameter :: io_spectral = 1
  integer, parameter :: io_vertical = 2
  integer, parameter :: io_single_record = 4

  !/**********************************************************************|
  !| DERIVED TYPES                                                        |
  !\**********************************************************************/
  type nc_ids_type
    !/********************************************************************|
    !| NC_IDS_TYPE                                                        |
    !| -----------                                                        |
    !| This derived type contains the NetCDF IDs for the variables output |
    !| in any NetCDF file. These are not intended to be manipulated by    |
    !| the user but rather controlled by open_file and used by            |
    !| write_file.                                                        |
    !\********************************************************************/
    character(len=100) :: file_name ! The name of the file
    integer :: ncid ! The netcdf ID of the file
    integer :: time_dimid !| The netcdf ID of the time dimension
    integer :: x_dimid,y_dimid,z_dimid ! The netcdf IDs of the dimensions
    integer :: x,y,z ! The netcdf IDs of the x, y, and z variables
    integer :: kx,ky,kz ! The netcdf IDs of the wavenumber variables
    integer :: time,timestep,dt ! The netcdf IDs of the time variables
    integer :: rho,en,c,u,v,w,p,t ! The netcdf IDs of the basic arrays
    integer :: rhoaux,enaux,caux,uaux,vaux,waux,paux,taux ! The netcdf IDs of the aux arrays
    integer :: flags ! The io flags of the output file
    integer :: step ! The current step of the output
    integer :: dims ! The number of dimensions in the output
    integer :: dimids(4) ! The array of netcdf dimension IDs
  end type nc_ids_type

  type controllers_type
    !/********************************************************************|
    !| CONTROLLERS_TYPE                                                   |
    !| ----------------                                                   |
    !| This derived type contains the parameters governing the IO         |
    !| routines. It contains the handles to the relevant data files and   |
    !| the information on the file names and steps in between outputs.    |
    !| This derived type is typically read by open_files in a namelist    |
    !| called io, but the object can be edited manually prior to calling  |
    !| open_files instead.                                                |
    !\********************************************************************/
    type(nc_ids_type) :: restart ! The handle on the restart file
    type(nc_ids_type) :: data ! The handle on the data file
    type(nc_ids_type) :: profile ! The handle on the profile file

    integer :: steps_diag = 1 ! The number of steps between diag files
    integer :: steps_data = 100 ! The number of steps between data files
    integer :: steps_restart = 100 ! The number of steps between restart files
    integer :: steps_profile = 100 ! The number of steps between profile files
    
    character(len=100) :: file_diag = "diag.csv" ! The name of the diag file
    character(len=100) :: file_data = "data.nc" ! The name of the data file
    character(len=100) :: file_profile = "profiles.nc" ! The name of the profile file
    character(len=100) :: file_restart = "restart.nc" ! The name of the restart file

    character(len=100) :: file_input = "" ! The name of the input file

    integer :: unit_diag = 79 ! The unit of the diag file
  end type controllers_type

contains
  !/**********************************************************************|
  !| FUNCTIONS                                                            |
  !\**********************************************************************/
  function mean_flux(grid,x,u,factor)
    !/********************************************************************|
    !| MEAN_FLUX                                                          |
    !| --------                                                           |
    !| Compute the mean turbulent flux of quantity x in the direction of  |
    !| velocity component u. The horizontal mean of both x and u are      |
    !| subtracted to get the local perturbation of each. An optional      |
    !| factor may be specified as well. The mean is only taken over the   |
    !| region spanning from 1/8 to 3/8 of the domain in z.
    !|   grid: a grid_type object that describes the grid                 | 
    !|   x: the physical array for which the turbulent mean is desired    |
    !|   u: the physical array for the component of velocity              |
    !|   factor: an optional physical array containing a scaling quantity |
    !|   returns: the average turbulent flux of x in the domain           |
    !\********************************************************************/
    use comm_module, only: get_comm,collect_sum
    implicit none

    type(grid_type) :: grid
    real(kind=dp) :: x(:,:,:)
    real(kind=dp) :: u(:,:,:)
    real(kind=dp) :: x_dx_mean,x_dz_mean
    real(kind=dp) :: mean_flux
    real(kind=dp), optional :: factor(:,:,:)

    real(kind=dp) :: rtmp,vel,vol
    integer :: i,j,k,ierror

    mean_flux = 0.0_dp
    vol = 0.0_dp

    if (present(factor)) then
      do k = 1,grid%z_p_count
        if (grid%lz(k) < 0.125_dp * grid%z_length) cycle
        if (grid%lz(k) > 0.375_dp * grid%z_length) cycle
        do j = 1,grid%y_p_count
          vel = sum(u(:,j,k)) / grid%x_n
          rtmp = sum(x(:,j,k)) / grid%x_n
          do i = 1,grid%x_n
            mean_flux = mean_flux + factor(i,j,k) * grid%dv * (u(i,j,k) - vel) * (x(i,j,k) - rtmp)
            vol = vol + grid%dv
          end do
        end do
      end do
    else
      do k = 1,grid%z_p_count
        if (grid%lz(k) < 0.125_dp * grid%z_length) cycle
        if (grid%lz(k) > 0.375_dp * grid%z_length) cycle
        do j = 1,grid%y_p_count
          vel = sum(u(:,j,k)) / grid%x_n
          rtmp = sum(x(:,j,k)) / grid%x_n
          do i = 1,grid%x_n
            mean_flux = mean_flux + grid%dv * (u(i,j,k) - vel) * (x(i,j,k) - rtmp)
            vol = vol + grid%dv
          end do
        end do
      end do
    end if

    mean_flux = collect_sum(mean_flux)
    mean_flux = mean_flux / collect_sum(vol)
  end function

  function dissipation(grid,x)
    !/********************************************************************|
    !| DISSIPATION                                                        |
    !| -----------                                                        |
    !| This function returns the mean dissipation a quantity x, that is,  |
    !| it returns (gradient(x))^2, where the square indicates the inner   |
    !| product of a vector with itself. This is useful in understanding   |
    !| the turbulence in the system and can serve as a proxy for the flux |
    !| in some circumstances.                                             |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   x: the physical array for which the turbulent mean is desired    |
    !|   returns: the average dissipation of x in the domain              |
    !\********************************************************************/
    use grid_module, only: allocate_s
    use grid_ops_module, only: deriv_x,deriv_y,deriv_z
    use fields_module, only: sum_squares
    implicit none

    type(grid_type) :: grid
    complex(kind=dp) :: x(:,:,:)

    complex(kind=dp), pointer :: work_s(:,:,:)
    real(kind=dp) :: dissipation

    call allocate_s(grid,work_s)
    
    call deriv_x(grid,x,work_s)
    dissipation = sum_squares(grid,work_s)

    call deriv_y(grid,x,work_s)
    dissipation = dissipation + sum_squares(grid,work_s)

    call deriv_z(grid,x,work_s)
    dissipation = dissipation + sum_squares(grid,work_s)

    dissipation = dissipation / grid%xyz_n

    deallocate(work_s)
  end function

  function ke(grid,state)
    !/********************************************************************|
    !| KE                                                                 |
    !| --                                                                 |
    !| This function calculates the total kinetic energy within the       |
    !| domain of the simulation.                                          |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   state: a state_type object that contains the simulation data     |
    !|   returns: the total kinetic energy contained in the domain        |
    !\********************************************************************/
    use comm_module, only: get_comm,collect_sum
    implicit none

    type(grid_type) :: grid
    type(state_type) :: state
    real(kind=dp) :: ke

    real(kind=dp) :: rtmp
    integer :: i,j,k,ierror

    ke = 0.0_dp

    do k = 1,grid%z_p_count
      do j = 1,grid%y_p_count
        do i = 1,grid%x_n
          ke = ke + state%rhou%p(i,j,k) * state%u%p(i,j,k)
          ke = ke + state%rhov%p(i,j,k) * state%v%p(i,j,k)
          ke = ke + state%rhow%p(i,j,k) * state%w%p(i,j,k)
        end do
      end do
    end do

    ke = grid%dv * collect_sum(ke) / 2.0_dp
  end function

  function pe(grid,params,state)
    !/********************************************************************|
    !| PE                                                                 |
    !| --                                                                 |
    !| This function calculates the total potential energy within the     |
    !| domain of the simulation.                                          |
    !|   grid: a grid_type object that describes the grid                 |
    !|   params: a params_type object that specifies the parameters       |  
    !|   state: a state_type object that contains the simulation data     |
    !|   returns: the total potential energy contained in the domain      |
    !\********************************************************************/
    use comm_module, only: get_comm,collect_sum
    implicit none

    type(grid_type) :: grid
    type(params_type) :: params
    type(state_type) :: state
    real(kind=dp) :: pe

    real(kind=dp) :: rtmp
    integer :: i,j,k,ierror

    pe = 0.0_dp

    do k = 1,grid%z_p_count
      rtmp = sum(grid%g(:k+grid%z_p_offset)) * grid%dz
      do j = 1,grid%y_p_count
        do i = 1,grid%x_n
          pe = pe - (state%rho%p(i,j,k) - params%rho_ref) * rtmp
        end do
      end do
    end do

    pe = grid%dv * params%buoy * collect_sum(pe)
  end function

  function ie(grid,params,state)
    !/********************************************************************|
    !| IE                                                                 |
    !| --                                                                 |
    !| This function calculates the total internal energy within the      |
    !| domain of the simulation.                                          |
    !|   grid: a grid_type object that describes the grid                 |
    !|   params: a params_type object that specifies the parameters       |  
    !|   state: a state_type object that contains the simulation data     |
    !|   returns: the total internal energy contained in the domain       |
    !\********************************************************************/
    use comm_module, only: get_comm,collect_sum,abort_comm
    use params_module, only: eos_ideal,eos_linear
    implicit none

    type(grid_type) :: grid
    type(params_type) :: params
    type(state_type) :: state
    real(kind=dp) :: ie

    real(kind=dp) :: rtmp,rtmp1,rtmp2,rtmp3,rtmp4,en_ref,ratio
    integer :: i,j,k,ierror

    ie = 0.0_dp

    if (params%eosin == eos_ideal) then
      rtmp = params%rho_ref * params%cv * params%gamma * params%t_ref
      do k = 1,grid%z_p_count
        do j = 1,grid%y_p_count
          do i = 1,grid%x_n
            ie = ie + state%cp%p(i,j,k) * state%t%p(i,j,k) - rtmp
          end do
        end do
      end do
      ie = ie / params%gamma
    else if (params%eosin == eos_linear) then
      rtmp1 = 1.0_dp / params%rho_ref - params%beta_t * params%t_ref / params%rho_ref
      rtmp2 = 2.0_dp * params%cs**2 * params%rho_ref**2
      rtmp3 = params%rho_ref * params%alr / 2.0_dp / params%beta_c
      en_ref = (params%beta_t * params%t_ref / params%rho_ref + params%p_ref / params%cs**2 / params%rho_ref**2) / params%alr
      rtmp4 = params%t_ref - rtmp3 * ratio * en_ref 
      do k = 1,grid%z_p_count
        do j = 1,grid%y_p_count
          do i = 1,grid%x_n
            rtmp = rtmp1 * state%p%p(i,j,k) - (state%p%p(i,j,k) + params%p_ref)**2 / rtmp2
            rtmp = rtmp + (state%t%p(i,j,k) - params%t_ref - rtmp3 * state%en%p(i,j,k)) * (state%en%p(i,j,k) + en_ref)
            rtmp = rtmp + params%p_ref**2 / rtmp2 + rtmp4 * state%en%p(i,j,k)
            ie = ie + state%rho%p(i,j,k) * rtmp - state%p%p(i,j,k)
          end do
        end do
      end do
    else
      print*, "FATAL: Invalid eos model on call to IE"
      call abort_comm()
    end if

    ie = grid%dv * collect_sum(ie)
  end function

  !/**********************************************************************|
  !| SUBROUTINES                                                          |
  !\**********************************************************************/
  subroutine open_files(grid,params,controls)
    !/********************************************************************|
    !| OPEN_FILES                                                         |
    !| ----------                                                         |
    !| This subroutine reads in the io parameters and opens all output    |
    !| files.                                                             |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   controls: a controllers_type object that handles io parameters   |
    !\********************************************************************/
    use comm_module, only: get_id,abort_comm
    implicit none
    
    type(grid_type) :: grid
    type(params_type) :: params
    type(controllers_type) :: controls
    
    ! Open the diagnostic file and write the header
    if (get_id() == 0) then 
      open(unit=controls%unit_diag,FILE=controls%file_diag,action="write")
      write(controls%unit_diag,*) "# timestep,time,dt,", &
        & "rms_u,rms_v,rms_w,rms_t,rms_c,", &
        & "flux_t,flux_c,flux_en,diss_t,diss_c,", &
        & "eflux_t,eflux_c,eflux_en,", &
        & "ke,pe,ie,", &
        & "probe_w,probe_t,probe_c,probe_dt,probe_dc"
    end if

    ! Open the restart file
    call open_file(grid,params,controls%restart,controls%file_restart,ior(io_single_record,io_spectral))

    ! Open the data file
    call open_file(grid,params,controls%data,controls%file_data)
    
    ! Open the profile file
    call open_file(grid,params,controls%profile,controls%file_profile,io_vertical)
  end subroutine

  subroutine close_files(controls)
    !/********************************************************************|
    !| CLOSE_FILES                                                        |
    !| -----------                                                        |
    !| This subroutine closes all output files.                           |
    !|   controls: a controllers_type object that handles io parameters   |
    !\********************************************************************/
    use comm_module, only: get_id
    implicit none

    type(controllers_type) :: controls
    integer :: i,ierror
  
    if (get_id() == 0) then
        CLOSE(UNIT=controls%unit_diag)
    end if
  
    call close_file(controls%restart)
    call close_file(controls%data)
    call close_file(controls%profile)
  end subroutine

  subroutine write_output_files(grid,params,controls,state)
    !/********************************************************************|
    !| WRITE_FILES                                                        |
    !| ----------                                                         |
    !| This subroutine evaluates whether a step should be output to the   |
    !| output files and then calls the relevant subroutine to write that  |
    !| step to file.                                                      |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   controls: a controllers_type object that handles io parameters   |
    !|   state: a state_type object that contains the simulation data     |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(controllers_type), intent(in) :: controls
    type(state_type), intent(in) :: state
  
    ! Write to the diafnostics file
    if (MOD(state%timestep,controls%steps_diag) == 0) then
      call write_diagnostics(grid,params,controls%unit_diag,state)
    end if

    ! Write to the profile file
    if (MOD(state%timestep,controls%steps_profile) == 0) then
      call write_file(grid,params,controls%profile,state)
    end if

    ! Write to the restart file
    if (MOD(state%timestep,controls%steps_restart) == 0) then
      call write_file(grid,params,controls%restart,state)
    end if

    ! Write to the data file
    if (MOD(state%timestep,controls%steps_data) == 0) then
      call write_file(grid,params,controls%data,state)
    end if
  end subroutine

  subroutine write_diagnostics(grid,params,unit,state)
    !/********************************************************************|
    !| WRITE_DIAGNOSTICS                                                  |
    !| -----------------                                                  |
    !| This subroutine computes various diagnostics regarding the state   |
    !| of the simulation and then writes those to a CSV file.             |
    !|   grid: a grid_type object that describes the grid                 |
    !|   params: a params_type object that specifies the parameters       | 
    !|   unit: the unit of the output file to write to                    |
    !|   state: a state_type object that contains the simulation data     |
    !\********************************************************************/
    use comm_module, only: get_id,get_comm
    use grid_ops_module, only: div
    use fields_module, only: sum_squares,sum_squares_vector,sum_product_p
    implicit none

    type(grid_type) :: grid
    type(params_type) :: params
    type(state_type) :: state
    integer :: unit

    real(kind=dp) :: flux_t,flux_c,flux_en,diss_t,diss_c
    real(kind=dp) :: eflux_t,eflux_c,eflux_en
    real(kind=dp) :: ke_sum,pe_sum,ie_sum
    real(kind=dp) :: rms_u,rms_v,rms_w,rms_t,rms_c,prod_uw,prod_vw
    real(kind=dp) :: probe_w,probe_t,probe_c,probe_dt,probe_dc
    real(kind=dp) :: rtmp
    integer :: k,ierr
  
    ! Calculate the vertical turbulent fluxes of the temperature and concentration
    flux_t = mean_flux(grid,state%t%p,state%w%p)
    flux_c = mean_flux(grid,state%c%p,state%w%p)
    flux_en = mean_flux(grid,state%t%p,state%w%p)

    eflux_t = mean_flux(grid,state%t%p,state%w%p,state%cp%p)
    eflux_c = mean_flux(grid,state%c%p,state%w%p,state%rho%p*state%mu%p)
    eflux_en = mean_flux(grid,state%en%p,state%w%p,state%rho%p*state%T%p)

    ! Calculate the dissipation of the temperature and concentration
    diss_t = dissipation(grid,state%t%s)
    diss_c = dissipation(grid,state%c%s)

    ! Calculate the rms of velocity, temperature (T), concentration (C), and C - T
    rms_u = sqrt(sum_squares(grid,state%u%s))
    rms_v = sqrt(sum_squares(grid,state%v%s))
    rms_w = sqrt(sum_squares(grid,state%w%s))
    rms_t = sqrt(sum_squares(grid,state%t%s))
    rms_c = sqrt(sum_squares(grid,state%c%s))    

    ke_sum = ke(grid,state)
    pe_sum = pe(grid,params,state)
    ie_sum = ie(grid,params,state)
    
    probe_w = 0.0_dp
    probe_t = 0.0_dp
    probe_c = 0.0_dp
    probe_dt = 0.0_dp
    probe_dc = 0.0_dp
    rtmp = grid%z_length / 4.0_dp
    if (grid%ly(1) == 0.0_dp) then
      if (grid%lz(1) <= rtmp .and. grid%lz(2) > rtmp) then
        probe_w = state%w%p(1,1,1)
        probe_t = state%t%p(1,1,1)
        probe_c = state%c%p(1,1,1)
        probe_dt = (state%t%p(1,1,2) - state%t%p(1,1,1)) / grid%dz
        probe_dc = (state%c%p(1,1,2) - state%c%p(1,1,1)) / grid%dz
      end if
      do k = 2,grid%z_p_count
        if (grid%lz(k) <= rtmp .and. grid%lz(k) + grid%dz > rtmp) then
          probe_w = state%w%p(1,1,k)
          probe_t = state%t%p(1,1,k)
          probe_c = state%c%p(1,1,k)
          probe_dt = (state%t%p(1,1,k) - state%t%p(1,1,k - 1)) / grid%dz
          probe_dc = (state%c%p(1,1,k) - state%c%p(1,1,k - 1)) / grid%dz
        end if
      end do
    end if

    ! Write to the diagnostics file
    if (get_id() == 0) then
      write(unit,'(i7,23(",",1E21.15))') &
        & state%timestep,state%time,state%dt, &
        & rms_u,rms_v,rms_w,rms_t,rms_c, &
        & flux_t,flux_c,flux_en,diss_t,diss_c, &
        & eflux_t,eflux_c,eflux_en, &
        & ke_sum,pe_sum,ie_sum, &
        & probe_w,probe_t,probe_c,probe_dt,probe_dc
    end if
  end subroutine

#include "output_netcdf.F90"

end module