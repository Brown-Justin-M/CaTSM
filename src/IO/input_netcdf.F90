  !/**********************************************************************|
  !| FUNCTIONS                                                            |
  !\**********************************************************************/
  function idefault(ncid,variable,val)
    !/********************************************************************|
    !| IDEFAULT                                                           |
    !| --------                                                           |
    !| A helper function that tries to read a time-dependent variable     |
    !| from an open netCDF file. If the variable doesn't exist, a default |
    !| value (val) is used instead.                                       |
    !|   ncid: the netCDF file ID associated with the file                | 
    !|   variable: the character string representation of the variable    |
    !|   val: the default value to use on failure                         |
    !\********************************************************************/
    use mpi
    use netcdf
    implicit none

    integer :: ncid
    character(*) :: variable
    integer :: val
    integer :: idefault

    integer :: varid,dimid
    integer :: step(1)
    character(len=100) :: ctmp
    integer :: itmp(1),ierror,status

    ierror = nf90_inq_dimid(ncid,"time",dimid)

    status = nf90_inq_varid(ncid,variable,varid)
    if (status == nf90_enotvar) then
      idefault = val
    else if (status == nf90_noerr) then
      if (ierror == nf90_ebaddim) then
        call check(nf90_get_var(ncid,varid,idefault))
      else
        call check(nf90_inquire_dimension(ncid,dimid,ctmp,step(1)))
        call check(nf90_get_var(ncid,varid,idefault,start=step))
      endif
    else
      call check(status)
    endif
  end function

  function rdefault(ncid,variable,val)
    !/********************************************************************|
    !| RDEFAULT                                                           |
    !| --------                                                           |
    !| A helper function that tries to read a time-dependent variable     |
    !| from an open netCDF file. If the variable doesn't exist, a default |
    !| value (val) is used instead.                                       |
    !|   ncid: the netCDF file ID associated with the file                | 
    !|   variable: the character string representation of the variable    |
    !|   val: the default value to use on failure                         |
    !\********************************************************************/
    use mpi
    use netcdf
    implicit none

    integer :: ncid
    character(*) :: variable
    real(kind=dp) :: val
    real(kind=dp) :: rdefault

    integer :: varid,dimid
    integer :: step(1)
    real(kind=dp) :: rtmp
    character(len=100) :: ctmp
    integer :: ierror,status

    ierror = nf90_inq_dimid(ncid,"time",dimid)

    status = nf90_inq_varid(ncid,variable,varid)
    if (status == nf90_enotvar) then
      rdefault = val
    else if (status == nf90_noerr) then
      if (ierror == nf90_ebaddim) then
        call check(nf90_get_var(ncid,varid,rdefault))
      else
        call check(nf90_inquire_dimension(ncid,dimid,ctmp,step(1)))
        call check(nf90_get_var(ncid,varid,rdefault,start=step))
      endif
    else
      call check(status)
    endif
  end function

  !/**********************************************************************|
  !| SUBROUTINES                                                          |
  !\**********************************************************************/
  subroutine read_input(grid,file_name,state)
    !/********************************************************************|
    !| READ_INPUT                                                         |
    !| ----------                                                         |
    !| This subroutine reads in an input file, determines its type, and   |
    !| checks that the input file is valid for the problem. It then calls |
    !| the relevant read subroutine (read_spectral or read_data) to read  |
    !| the data into state.                                               |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   file_name: the name of the input file to read in                 |
    !|   state: the state_type object containing the simulation data      |
    !\********************************************************************/
    use mpi
    use netcdf
    use output_module
    use comm_module, only: get_id,abort_comm
    use fields_module, only: transform_all_to_s,transform_all_to_p
    use transform_module, only: transform_s_to_p,transform_p_to_s
    implicit none

    type(grid_type), intent(inout) :: grid
    character(len=100), intent(in) :: file_name
    type(state_type), intent(inout) :: state

    type(nc_ids_type) :: nc_ids
    integer :: itmp,ierror
    character(len=100) :: ctmp

    call check(nf90_open(trim(file_name),nf90_nowrite,nc_ids%ncid))

    state%timestep = idefault(nc_ids%ncid,"timestep",0)
    state%restart = state%timestep
    state%time = rdefault(nc_ids%ncid,"time",0.0_dp)
    state%dt = rdefault(nc_ids%ncid,"dt",0.0_dp)

    ! Try to open the input file in spectral space
    ierror = nf90_inq_dimid(nc_ids%ncid,"kz",nc_ids%z_dimid)
    if (ierror == nf90_noerr) then
      ! The input file is spectral
      if (get_id() == 0) then
        print*,"Interpreting input file as spectral"
      end if

      ! Check that the dimensions of the input file match the simulation
      call check(nf90_inquire_dimension(nc_ids%ncid,nc_ids%z_dimid,ctmp,itmp))
      if (grid%z_modes /= itmp / 2) then
        print*,"FATAL: Z-dimension of restart file incompatible with setup"
        call abort_comm
      end if
      call check(nf90_inq_dimid(nc_ids%ncid,"kx",nc_ids%x_dimid))
      call check(nf90_inquire_dimension(nc_ids%ncid,nc_ids%x_dimid,ctmp,itmp))
      if (grid%x_modes /= itmp - 1) then
        print*,"FATAL: X-dimension of restart file incompatible with setup"
        call abort_comm
      end if
      call check(nf90_inq_dimid(nc_ids%ncid,"ky",nc_ids%y_dimid))
      call check(nf90_inquire_dimension(nc_ids%ncid,nc_ids%y_dimid,ctmp,itmp))
      if (grid%y_modes /= itmp / 2) then
        print*,"FATAL: Y-dimension of restart file incompatible with setup"
        call abort_comm
      end if

      call read_spectral(grid,nc_ids,state)
    else if (ierror == nf90_ebaddim) then
      ! The input file is physical
      if (get_id() == 0) then
        print*,"Interpreting input file as phsyical"
      end if

      ! Check that the dimensions of the input file match the simulation
      call check(nf90_inq_dimid(nc_ids%ncid,"x",nc_ids%x_dimid))
      call check(nf90_inquire_dimension(nc_ids%ncid,nc_ids%x_dimid,ctmp,itmp))
      if (grid%x_n /= itmp) then
        print*,"FATAL: X-dimension of data file incompatible with setup"
        call abort_comm
      end if
      call check(nf90_inq_dimid(nc_ids%ncid,"y",nc_ids%y_dimid))
      call check(nf90_inquire_dimension(nc_ids%ncid,nc_ids%y_dimid,ctmp,itmp))
      if (grid%y_n /= itmp) then
        print*,"FATAL: Y-dimension of data file incompatible with setup"
        call abort_comm
      end if
      call check(nf90_inq_dimid(nc_ids%ncid,"z",nc_ids%z_dimid))
      call check(nf90_inquire_dimension(nc_ids%ncid,nc_ids%z_dimid,ctmp,itmp))
      if (grid%z_n /= itmp) then
        print*,"FATAL: Z-dimension of data file incompatible with setup"
        call abort_comm
      end if

      call read_data(grid,nc_ids,state)
    else
      call check(ierror)
    end if
    
    call check(nf90_close(nc_ids%ncid))

    call transform_s_to_p(state%u%s,state%u%p)
    call transform_s_to_p(state%v%s,state%v%p)
    call transform_s_to_p(state%w%s,state%w%p)
    call transform_s_to_p(state%rho%s,state%rho%p)

    state%rhou%p = state%rho%p * state%u%p
    state%rhov%p = state%rho%p * state%v%p
    state%rhow%p = state%rho%p * state%w%p

    call transform_p_to_s(state%rhou%p,state%rhou%s)
    call transform_p_to_s(state%rhov%p,state%rhov%s)
    call transform_p_to_s(state%rhow%p,state%rhow%s)

    ! Repopulate the physical space arrays
    call transform_all_to_p(grid,state)
  end subroutine
  
  subroutine read_spectral(grid,nc_ids,state)
    !/********************************************************************|
    !| READ_SPECTRAL                                                      |
    !| -------------                                                      |
    !| This subroutine reads in an input file in spectral coordinates. It |
    !| uses the file to populate the arrays and time information in the   |
    !| state object.                                                      |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   nc_ids: the nc_ids_type object containing the netcdf ids         |
    !|   state: the state_type object containing the simulation data      |
    !\********************************************************************/
    use netcdf
    use output_module
    use comm_module, only: get_id
    use output_module, only: check
    implicit none

    type(grid_type), intent(in) :: grid
    type(nc_ids_type), intent(inout) :: nc_ids
    type(state_type), intent(inout) :: state

    complex(kind=dp),parameter :: iu = (0.0_dp,1.0_dp)
    integer :: itmp,ierror,step(4)
    character(len=100) :: ctmp
    real(kind=dp), allocatable :: rwork(:,:,:)

    ! Determine if the input file has a time coordinate
    ! If it does, we want to read in the latest time
    itmp = 1
    ierror = nf90_inq_dimid(nc_ids%ncid,"time",nc_ids%time_dimid)
    if (ierror == nf90_noerr) then
      call check(nf90_inquire_dimension(nc_ids%ncid,nc_ids%time_dimid,ctmp,itmp))
    else if (ierror /= nf90_ebaddim) then
      call check(ierror)
    endif

    step = (/1,1,1,itmp/)

    ! Check for the real and imaginary components of the basic arrays
    call check(nf90_inq_varid(nc_ids%ncid,"RHO",nc_ids%rho))
    call check(nf90_inq_varid(nc_ids%ncid,"ETA",nc_ids%en))
    call check(nf90_inq_varid(nc_ids%ncid,"C",nc_ids%c))
    call check(nf90_inq_varid(nc_ids%ncid,"U",nc_ids%u))
    call check(nf90_inq_varid(nc_ids%ncid,"V",nc_ids%v))
    call check(nf90_inq_varid(nc_ids%ncid,"W",nc_ids%w))
    call check(nf90_inq_varid(nc_ids%ncid,"RHOI",nc_ids%rhoaux))
    call check(nf90_inq_varid(nc_ids%ncid,"EI",nc_ids%enaux))
    call check(nf90_inq_varid(nc_ids%ncid,"CI",nc_ids%caux))
    call check(nf90_inq_varid(nc_ids%ncid,"UI",nc_ids%uaux))
    call check(nf90_inq_varid(nc_ids%ncid,"VI",nc_ids%vaux))
    call check(nf90_inq_varid(nc_ids%ncid,"WI",nc_ids%waux))

    ! Read the real and imaginary components into a work array
    ! and then direct them to the real and imaginary parts of state
    allocate(rwork(0:2 * grid%z_modes - 1,0:grid%x_modes,0:max(2 * grid%y_modes - 1,0)))
    call check(nf90_get_var(nc_ids%ncid,nc_ids%rho,rwork))
    state%rho%s = rwork
    call check(nf90_get_var(nc_ids%ncid,nc_ids%rhoaux,rwork))
    state%rho%s = state%rho%s + iu * rwork
    call check(nf90_get_var(nc_ids%ncid,nc_ids%en,rwork))
    state%en%s = rwork
    call check(nf90_get_var(nc_ids%ncid,nc_ids%enaux,rwork))
    state%en%s = state%en%s + iu * rwork
    call check(nf90_get_var(nc_ids%ncid,nc_ids%c,rwork))
    state%c%s = rwork
    call check(nf90_get_var(nc_ids%ncid,nc_ids%caux,rwork))
    state%c%s = state%c%s + iu * rwork
    call check(nf90_get_var(nc_ids%ncid,nc_ids%u,rwork))
    state%u%s = rwork
    call check(nf90_get_var(nc_ids%ncid,nc_ids%uaux,rwork))
    state%u%s = state%u%s + iu * rwork
    call check(nf90_get_var(nc_ids%ncid,nc_ids%v,rwork))
    state%v%s = rwork
    call check(nf90_get_var(nc_ids%ncid,nc_ids%vaux,rwork))
    state%v%s = state%v%s + iu * rwork
    call check(nf90_get_var(nc_ids%ncid,nc_ids%w,rwork))
    state%w%s = rwork
    call check(nf90_get_var(nc_ids%ncid,nc_ids%waux,rwork))
    state%w%s = state%w%s + iu * rwork
    deallocate(rwork)
  end subroutine

  subroutine read_data(grid,nc_ids,state)
    !/********************************************************************|
    !| READ_DATA                                                          |
    !| ---------                                                          |
    !| This subroutine reads in an input file in physical coordinates. It |
    !| uses the file to populate the arrays and time information in the   |
    !| state object.                                                      |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   nc_ids: the nc_ids_type object containing the netcdf ids         |
    !|   state: the state_type object containing the simulation data      |
    !\********************************************************************/
    use netcdf
    use comm_module, only: get_id
    use fields_module, only: transform_all_to_s
    use output_module, only: check
    implicit none

    type(grid_type) :: grid
    type(nc_ids_type) :: nc_ids
    type(state_type) :: state

    integer :: itmp,ierror,step(4)
    character(len=100) :: ctmp

    ! Determine if the input file has a time coordinate
    ! If it does, we want to read in the latest time
    itmp = 1
    ierror = nf90_inq_dimid(nc_ids%ncid,"time",nc_ids%time_dimid)
    if (ierror == nf90_noerr) then
      call check(nf90_inquire_dimension(nc_ids%ncid,nc_ids%time_dimid,ctmp,itmp))
    else if (ierror /= nf90_ebaddim) then
      call check(ierror)
    endif

    step = (/1,1,1,itmp/)

    ! Check for the basic arrays and read them in directly
    call check(nf90_inq_varid(nc_ids%ncid,"RHO",nc_ids%rho))
    call check(nf90_inq_varid(nc_ids%ncid,"ETA",nc_ids%en))
    call check(nf90_inq_varid(nc_ids%ncid,"C",nc_ids%c))
    call check(nf90_inq_varid(nc_ids%ncid,"U",nc_ids%u))
    call check(nf90_inq_varid(nc_ids%ncid,"V",nc_ids%v))
    call check(nf90_inq_varid(nc_ids%ncid,"W",nc_ids%w))

    call check(nf90_get_var(nc_ids%ncid,nc_ids%rho,state%rho%p,start=step))
    call check(nf90_get_var(nc_ids%ncid,nc_ids%en,state%en%p,start=step))
    call check(nf90_get_var(nc_ids%ncid,nc_ids%c,state%c%p,start=step))
    call check(nf90_get_var(nc_ids%ncid,nc_ids%u,state%u%p,start=step))
    call check(nf90_get_var(nc_ids%ncid,nc_ids%v,state%v%p,start=step))
    call check(nf90_get_var(nc_ids%ncid,nc_ids%w,state%w%p,start=step))

    ! Transform the simulation back into spectral space
    call transform_all_to_s(grid,state)
  end subroutine
