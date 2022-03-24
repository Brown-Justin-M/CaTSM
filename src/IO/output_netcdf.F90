  subroutine open_file(grid,params,nc_ids,file_name,flags)
    !/********************************************************************|
    !| OPEN_FILE                                                          |
    !| ---------                                                          |
    !| This subroutine opens a new netCDF file, the properties of which   |
    !| are specified by flags. **WARNING** Though it is possible to       |
    !| specify a vertical spectrum with the flags, this is currently not  |
    !| implemented, and the code may fail.                                |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   nc_ids: a nc_ids_type object that contains the netCDF IDs        |
    !|   file_name: the name of the output file to open                   |
    !|   flags: the flags that define the type of output file             |
    !|     io_spectral: if specified, the output is in spectral space     |
    !|     io_vertical: if specified, the x and y dimensions are omitted  |
    !|     io_single_record: if specified, the time dimension is omitted  |
    !\********************************************************************/
    use netcdf
    use comm_module, only: get_id,abort_comm

    type(grid_type), intent(in) :: grid
    type(params_type), intent(in) :: params
    type(nc_ids_type), intent(inout) :: nc_ids
    character(len=100), intent(in) :: file_name
    integer, optional, intent(in) :: flags
    integer :: flagsin

    integer :: mode_flags
    integer, allocatable :: dimids(:)
    integer :: i

    flagsin = 0
    if (present(flags)) flagsin = flags

    nc_ids%file_name = file_name
    nc_ids%flags = flagsin

    mode_flags = ior(nf90_clobber,nf90_netcdf4)
    call check(nf90_create(trim(nc_ids%file_name),mode_flags,nc_ids%ncid))

    ! Check if the file is a single record file to be overwritten or a time series
    ! Set up the dimensions of the file accordingly
    if (iand(flagsin,io_single_record) == io_single_record) then
      nc_ids%dims = 0
      allocate(dimids(0))
    else
      nc_ids%dims = 1
      allocate(dimids(1))
      call check(nf90_def_dim(nc_ids%ncid,"time",nf90_unlimited,nc_ids%time_dimid))
      dimids = (/nc_ids%time_dimid/)
    end if

    ! Add the critical metadata variables to the file
    call check(nf90_def_var(nc_ids%ncid,"timestep",nf90_int,dimids,nc_ids%timestep))
    call check(nf90_def_var(nc_ids%ncid,"time",nf90_double,dimids,nc_ids%time))
    call check(nf90_def_var(nc_ids%ncid,"dt",nf90_double,dimids,nc_ids%dt))
    deallocate(dimids)

    mode_flags = ior(io_spectral,io_vertical)
    ! Open a 3D spectral file
    if (iand(flagsin,mode_flags) == io_spectral) then
      call open_spectral(grid,nc_ids)
    end if

    ! Open a vertical spectral file
    ! TODO: implement vertical spectral file
    if (iand(flagsin,mode_flags) == mode_flags) then
      print*,"WARNING: Vertical spectral output not yet implemented"
    end if

    ! Open a 3D physical file
    if (iand(flagsin,mode_flags) == 0) then
      call open_physical(grid,nc_ids)
    end if

    ! Open a vertical physical file
    if (iand(flagsin,mode_flags) == io_vertical) then
      call open_vertical(grid,nc_ids)
    end if

    ! Define the main variables
    allocate(dimids(nc_ids%dims))
    dimids(1:nc_ids%dims) = nc_ids%dimids(1:nc_ids%dims)
    call check(nf90_def_var(nc_ids%ncid,"RHO",nf90_double,dimids,nc_ids%rho))
    call check(nf90_def_var(nc_ids%ncid,"ETA",nf90_double,dimids,nc_ids%en))
    call check(nf90_def_var(nc_ids%ncid,"T",nf90_double,dimids,nc_ids%t))
    call check(nf90_def_var(nc_ids%ncid,"P",nf90_double,dimids,nc_ids%p))
    call check(nf90_def_var(nc_ids%ncid,"C",nf90_double,dimids,nc_ids%c))
    call check(nf90_def_var(nc_ids%ncid,"U",nf90_double,dimids,nc_ids%u))
    call check(nf90_def_var(nc_ids%ncid,"V",nf90_double,dimids,nc_ids%v))
    call check(nf90_def_var(nc_ids%ncid,"W",nf90_double,dimids,nc_ids%w))
    deallocate(dimids)

    call check(nf90_enddef(nc_ids%ncid))

    ! Write the dimension variable values
    if (iand(flagsin,io_spectral) == io_spectral) then
      call check(nf90_put_var(nc_ids%ncid,nc_ids%kx,grid%kx))
      call check(nf90_put_var(nc_ids%ncid,nc_ids%ky,grid%ky))
      call check(nf90_put_var(nc_ids%ncid,nc_ids%kz,grid%kz))
    else
      if (iand(flagsin,io_vertical) /= io_vertical) then
        call check(nf90_put_var(nc_ids%ncid,nc_ids%x,grid%x))
        call check(nf90_put_var(nc_ids%ncid,nc_ids%y,grid%y))
      end if
      call check(nf90_put_var(nc_ids%ncid,nc_ids%z,grid%z))
    end if

    call check(nf90_sync(nc_ids%ncid))

    nc_ids%step = 1
  end subroutine

  subroutine open_spectral(grid,nc_ids)
    !/********************************************************************|
    !| OPEN_SPECTRAL                                                      |
    !| -------------                                                      |
    !| This subroutine opens a new NetCDF file in spectral space.         |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   nc_ids: a nc_ids_type object that contains the netCDF IDs        |
    !\********************************************************************/
    use netcdf
    implicit none

    type(grid_type) :: grid
    type(nc_ids_type) :: nc_ids

    integer, parameter :: one = 1, two = 2
    integer, allocatable :: dimids(:)

    ! Define the wavenumber dimensions
    call check(nf90_def_dim(nc_ids%ncid,"kz",two * grid%z_modes,nc_ids%z_dimid))
    call check(nf90_def_dim(nc_ids%ncid,"kx",grid%x_modes + one,nc_ids%x_dimid))
    call check(nf90_def_dim(nc_ids%ncid,"ky",max(two * grid%y_modes,one), nc_ids%y_dimid))

    nc_ids%dims = 3 + nc_ids%dims
    nc_ids%dimids = (/nc_ids%z_dimid,nc_ids%x_dimid,nc_ids%y_dimid,nc_ids%time_dimid/)
    allocate(dimids(nc_ids%dims))
    dimids(1:nc_ids%dims) = nc_ids%dimids(1:nc_ids%dims)

    ! Define the wavenumber variables
    call check(nf90_def_var(nc_ids%ncid,"kz",nf90_double,dimids(1),nc_ids%kz))
    call check(nf90_def_var(nc_ids%ncid,"kx",nf90_double,dimids(2),nc_ids%kx))
    call check(nf90_def_var(nc_ids%ncid,"ky",nf90_double,dimids(3),nc_ids%ky))

    ! Define the imaginary components of the main variables
    call check(nf90_def_var(nc_ids%ncid,"RI",nf90_double,dimids,nc_ids%rhoaux))
    call check(nf90_def_var(nc_ids%ncid,"EI",nf90_double,dimids,nc_ids%enaux))
    call check(nf90_def_var(nc_ids%ncid,"TI",nf90_double,dimids,nc_ids%taux))
    call check(nf90_def_var(nc_ids%ncid,"PI",nf90_double,dimids,nc_ids%paux))
    call check(nf90_def_var(nc_ids%ncid,"CI",nf90_double,dimids,nc_ids%caux))
    call check(nf90_def_var(nc_ids%ncid,"UI",nf90_double,dimids,nc_ids%uaux))
    call check(nf90_def_var(nc_ids%ncid,"VI",nf90_double,dimids,nc_ids%vaux))
    call check(nf90_def_var(nc_ids%ncid,"WI",nf90_double,dimids,nc_ids%waux))
    deallocate(dimids)
  end subroutine

  subroutine open_physical(grid,nc_ids)
    !/********************************************************************|
    !| OPEN_PHYSICAL                                                      |
    !| -------------                                                      |
    !| This subroutine opens a new NetCDF file in physical space.         |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   nc_ids: a nc_ids_type object that contains the netCDF IDs        |
    !\********************************************************************/
    use netcdf
    implicit none

    type(grid_type) :: grid
    type(nc_ids_type) :: nc_ids

    integer, parameter :: one = 1
    integer, allocatable :: dimids(:)

    ! Define the spatial dimensions
    call check(nf90_def_dim(nc_ids%ncid,"x",one * grid%x_n,nc_ids%x_dimid))
    call check(nf90_def_dim(nc_ids%ncid,"y",one * grid%y_n,nc_ids%y_dimid))
    call check(nf90_def_dim(nc_ids%ncid,"z",one * grid%z_n,nc_ids%z_dimid))
    
    nc_ids%dims = 3 + nc_ids%dims
    nc_ids%dimids = (/nc_ids%x_dimid,nc_ids%y_dimid,nc_ids%z_dimid,nc_ids%time_dimid/)
    allocate(dimids(nc_ids%dims))
    dimids(1:nc_ids%dims) = nc_ids%dimids(1:nc_ids%dims)

    ! Define the spatial variables
    call check(nf90_def_var(nc_ids%ncid,"x",nf90_double,dimids(1),nc_ids%x))
    call check(nf90_def_var(nc_ids%ncid,"y",nf90_double,dimids(2),nc_ids%y))
    call check(nf90_def_var(nc_ids%ncid,"z",nf90_double,dimids(3),nc_ids%z))

    deallocate(dimids)
  end subroutine

  subroutine open_vertical(grid,nc_ids)
    !/********************************************************************|
    !| OPEN_VERTICAL                                                      |
    !| -------------                                                      |
    !| This subroutine opens a new NetCDF file for vertical profiles.     |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   nc_ids: a nc_ids_type object that contains the netCDF IDs        |
    !\********************************************************************/
    use netcdf
    implicit none

    type(grid_type) :: grid
    type(nc_ids_type) :: nc_ids

    integer, parameter :: one = 1
    integer :: dimids(2)

    ! Define the z dimension
    call check(nf90_def_dim(nc_ids%ncid,"z",one * grid%z_n,nc_ids%z_dimid))
    
    nc_ids%dims = 1 + nc_ids%dims
    nc_ids%dimids = (/nc_ids%z_dimid,nc_ids%time_dimid,0,0/)
    dimids = (/nc_ids%z_dimid,nc_ids%time_dimid/)

    call check(nf90_def_var(nc_ids%ncid,"z",nf90_double,dimids(1),nc_ids%z))

    ! Define variables to contain the horizontal means
    call check(nf90_def_var(nc_ids%ncid,"RMean",nf90_double,dimids,nc_ids%rhoaux))
    call check(nf90_def_var(nc_ids%ncid,"EMean",nf90_double,dimids,nc_ids%enaux))
    call check(nf90_def_var(nc_ids%ncid,"TMean",nf90_double,dimids,nc_ids%taux))
    call check(nf90_def_var(nc_ids%ncid,"PMean",nf90_double,dimids,nc_ids%paux))
    call check(nf90_def_var(nc_ids%ncid,"CMean",nf90_double,dimids,nc_ids%caux))
    call check(nf90_def_var(nc_ids%ncid,"UMean",nf90_double,dimids,nc_ids%uaux))
    call check(nf90_def_var(nc_ids%ncid,"VMean",nf90_double,dimids,nc_ids%vaux))
    call check(nf90_def_var(nc_ids%ncid,"WMean",nf90_double,dimids,nc_ids%waux))
  end subroutine

  subroutine close_file(nc_ids)
    !/********************************************************************|
    !| CLOSE_FILE                                                         |
    !| ----------                                                         |
    !| This subroutine closes a NetCDF file given its handle .            |
    !|   nc_ids: a nc_ids_type object that contains the netCDF IDs        |
    !\********************************************************************/
    use netcdf
    implicit none

    type(nc_ids_type) :: nc_ids

    call check(nf90_close(nc_ids%ncid))
  end subroutine

  subroutine write_file(grid,params,nc_ids,state)
    !/********************************************************************|
    !| WRITE_FILE                                                         |
    !| ----------                                                         |
    !| This subroutine writes to an opened netCDF file, using the data    |
    !| contained inside of state. This subroutine writes the metadata and |
    !| then calls the relevant subroutine for the particular type of      |
    !| file (physical, spectral, etc.).                                   |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   nc_ids: a nc_ids_type object that contains the netCDF IDs        |
    !|   state: a state_type object containing the simulation data        |
    !\********************************************************************/
    use netcdf
    use comm_module, only: get_id
    use grid_module, only: allocate_p
    implicit none

    type(grid_type) :: grid
    type(params_type) :: params
    type(nc_ids_type) :: nc_ids
    type(state_type) :: state
    
    integer :: flags
    integer, parameter :: one(1) = (/1/)
    integer :: step(1)

    ! Write the metadata
    step = (/nc_ids%step/)
    call check(nf90_put_var(nc_ids%ncid,nc_ids%timestep,state%timestep,start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%time,state%time,start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%dt,state%dt,start=step))

    ! Write to a spectral file
    flags = ior(io_spectral,io_vertical)
    if (iand(nc_ids%flags,flags) == io_spectral) then
      call write_spectral(grid,params,nc_ids,state)
    end if

    ! Write to a physical file
    if (iand(nc_ids%flags,flags) == 0) then
      call write_physical(grid,params,nc_ids,state)
    end if

    ! Write to a profile file
    if (iand(nc_ids%flags,flags) == io_vertical) then
      call write_vertical(grid,params,nc_ids,state)
    end if

    call check(nf90_sync(nc_ids%ncid))

    if (iand(nc_ids%flags,io_single_record) /= io_single_record) then
      nc_ids%step = nc_ids%step + 1
    end if
  end subroutine

  subroutine write_spectral(grid,params,nc_ids,state)
    !/********************************************************************|
    !| WRITE_SPECTRAL                                                     |
    !| --------------                                                     |
    !| This subroutine writes to an opened spectral netCDF file, using    |
    !| the data contained inside of state.                                |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   nc_ids: a nc_ids_type object that contains the netCDF IDs        |
    !|   state: a state_type object containing the simulation data        |
    !\********************************************************************/
    use netcdf
    implicit none

    type(grid_type) :: grid
    type(params_type) :: params
    type(nc_ids_type) :: nc_ids
    type(state_type) :: state

    call check(nf90_put_var(nc_ids%ncid,nc_ids%rho,real(state%rho%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%rhoaux,imag(state%rho%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%en,real(state%en%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%enaux,imag(state%en%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%t,real(state%t%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%taux,imag(state%t%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%p,real(state%p%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%paux,imag(state%p%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%c,real(state%c%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%caux,imag(state%c%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%u,real(state%u%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%uaux,imag(state%u%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%v,real(state%v%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%vaux,imag(state%v%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%w,real(state%w%s)))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%waux,imag(state%w%s)))
  end subroutine

  subroutine write_physical(grid,params,nc_ids,state)
    !/********************************************************************|
    !| WRITE_PHYSICAL                                                     |
    !| --------------                                                     |
    !| This subroutine writes to an opened physical netCDF file, using    |
    !| the data contained inside of state.                                |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   nc_ids: a nc_ids_type object that contains the netCDF IDs        |
    !|   state: a state_type object containing the simulation data        |
    !\********************************************************************/
    use netcdf
    use grid_module, only: allocate_p
    implicit none

    type(grid_type) :: grid
    type(params_type) :: params
    type(nc_ids_type) :: nc_ids
    type(state_type) :: state

    integer :: step(nc_ids%dims)
    real(kind=dp), pointer :: rwork(:,:,:)

    step = 1
    step(nc_ids%dims) = nc_ids%step

    call check(nf90_put_var(nc_ids%ncid,nc_ids%rho,state%rho%p,start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%en,state%en%p,start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%t,state%t%p,start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%p,state%p%p,start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%c,state%c%p,start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%u,state%u%p,start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%v,state%v%p,start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%w,state%w%p,start=step))
  end subroutine

  subroutine write_vertical(grid,params,nc_ids,state)
    !/********************************************************************|
    !| WRITE_VERTICAL                                                     |
    !| --------------                                                     |
    !| This subroutine writes to an opened profile netCDF file, using     |
    !| the data contained inside of state.                                |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   params: a params_type object that specifies the parameters       |
    !|   nc_ids: a nc_ids_type object that contains the netCDF IDs        |
    !|   state: a state_type object containing the simulation data        |
    !\********************************************************************/
    use netcdf
    use fields_module
    use transform_module, only: transform_xy_p_to_s,transform_xy_s_to_p
    use comm_module,only: get_id
    implicit none

    type(grid_type) :: grid
    type(params_type) :: params
    type(nc_ids_type) :: nc_ids
    type(state_type) :: state

    real(kind=dp) :: mean_profile(grid%z_n)
    integer :: step(nc_ids%dims)

    step = 1
    step(nc_ids%dims) = nc_ids%step

    ! Output the profiles at x = 0, y = 0
    call check(nf90_put_var(nc_ids%ncid,nc_ids%rho,state%rho%p(0,0,:),start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%en,state%en%p(0,0,:),start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%t,state%t%p(0,0,:),start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%p,state%p%p(0,0,:),start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%c,state%c%p(0,0,:),start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%u,state%u%p(0,0,:),start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%v,state%v%p(0,0,:),start=step))
    call check(nf90_put_var(nc_ids%ncid,nc_ids%w,state%w%p(0,0,:),start=step))

    ! Calculate the horizontal means and output them to the file
    mean_profile = sum(sum(state%rho%p,1),1) / grid%xy_n
    call check(nf90_put_var(nc_ids%ncid,nc_ids%rhoaux,mean_profile,start=step))
    mean_profile = sum(sum(state%en%p,1),1) / grid%xy_n
    call check(nf90_put_var(nc_ids%ncid,nc_ids%enaux,mean_profile,start=step))
    mean_profile = sum(sum(state%t%p,1),1) / grid%xy_n
    call check(nf90_put_var(nc_ids%ncid,nc_ids%taux,mean_profile,start=step))
    mean_profile = sum(sum(state%p%p,1),1) / grid%xy_n
    call check(nf90_put_var(nc_ids%ncid,nc_ids%paux,mean_profile,start=step))
    mean_profile = sum(sum(state%c%p,1),1) / grid%xy_n
    call check(nf90_put_var(nc_ids%ncid,nc_ids%caux,mean_profile,start=step))
    mean_profile = sum(sum(state%u%p,1),1) / grid%xy_n
    call check(nf90_put_var(nc_ids%ncid,nc_ids%uaux,mean_profile,start=step))
    mean_profile = sum(sum(state%v%p,1),1) / grid%xy_n
    call check(nf90_put_var(nc_ids%ncid,nc_ids%vaux,mean_profile,start=step))
    mean_profile = sum(sum(state%w%p,1),1) / grid%xy_n
    call check(nf90_put_var(nc_ids%ncid,nc_ids%waux,mean_profile,start=step))
  end subroutine

  subroutine check(ierror)
    !/********************************************************************|
    !| CHECK                                                              |
    !| -----                                                              |
    !| Check that the error code from a netcdf call is nf90_noerr (no     |
    !| errors). If not, print the error message and end the program.      |
    !|   ierror: the error code to check                                  |
    !\********************************************************************/
    use netcdf
    use comm_module, only: abort_comm
    implicit none

    integer,intent (in) :: ierror
      
    if (ierror /= nf90_noerr) then
      print*,"FATAL NETcdf: ",nf90_strerror(ierror)
      call abort_comm(ierror)
    end if
  end subroutine
