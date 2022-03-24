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
    use pnetcdf
    use fields_module, only: transform_all_to_s,transform_all_to_p
    use transform_module, only: transform_s_to_p,transform_p_to_s
    use comm_module, only: get_id,abort_comm
    implicit none

    type(grid_type), intent(inout) :: grid
    character(len=100), intent(in) :: file_name
    type(state_type), intent(inout) :: state

    type(nc_ids_type) :: nc_ids
    integer(kind=mpi_offset_kind) :: itmp
    integer :: ierror

    call check(nfmpi_open(mpi_comm_world,trim(file_name),nf_nowrite,mpi_info_null,nc_ids%ncid))

    state%timestep = idefault(nc_ids%ncid,"timestep",0)
    state%restart = state%timestep
    state%time = rdefault(nc_ids%ncid,"time",0.0_dp)
    state%dt = rdefault(nc_ids%ncid,"dt",0.0_dp)

    ierror = nfmpi_inq_dimid(nc_ids%ncid,"kz",nc_ids%z_dimid)
    if (ierror == nf_noerr) then
      if (get_id() == 0) then
        print*,"Interpreting input file as spectral"
      end if

      call check(nfmpi_inq_dimlen(nc_ids%ncid,nc_ids%z_dimid,itmp))
      if (grid%z_modes /= itmp / 2) then
        print*,"FATAL: Z-dimension of restart file incompatible with setup"
        call abort_comm
      end if
      call check(nfmpi_inq_dimid(nc_ids%ncid,"kx",nc_ids%x_dimid))
      call check(nfmpi_inq_dimlen(nc_ids%ncid,nc_ids%x_dimid,itmp))
      if (grid%x_modes /= itmp - 1) then
        print*,"FATAL: X-dimension of restart file incompatible with setup"
        call abort_comm
      end if
      call check(nfmpi_inq_dimid(nc_ids%ncid,"ky",nc_ids%y_dimid))
      call check(nfmpi_inq_dimlen(nc_ids%ncid,nc_ids%y_dimid,itmp))
      if (grid%y_modes /= itmp / 2) then
        print*,"FATAL: Y-dimension of restart file incompatible with setup"
        call abort_comm
      end if

      call read_spectral(grid,nc_ids,state)
    else if (ierror == nf_ebaddim) then
      if (get_id() == 0) then
        print*,"Interpreting input file as physical"
      end if

      call check(nfmpi_inq_dimid(nc_ids%ncid,"x",nc_ids%x_dimid))
      call check(nfmpi_inq_dimlen(nc_ids%ncid,nc_ids%x_dimid,itmp))
      if (grid%x_n /= itmp) then
        print*,"FATAL: X-dimension of data file incompatible with setup"
        call abort_comm
      end if
      call check(nfmpi_inq_dimid(nc_ids%ncid,"y",nc_ids%y_dimid))
      call check(nfmpi_inq_dimlen(nc_ids%ncid,nc_ids%y_dimid,itmp))
      if (grid%y_n /= itmp) then
        print*,"FATAL: Y-dimension of data file incompatible with setup"
        call abort_comm
      end if
      call check(nfmpi_inq_dimid(nc_ids%ncid,"z",nc_ids%z_dimid))
      call check(nfmpi_inq_dimlen(nc_ids%ncid,nc_ids%z_dimid,itmp))
      if (grid%z_n /= itmp) then
        print*,"FATAL: Z-dimension of data file incompatible with setup"
        call abort_comm
      end if

      call read_data(grid,nc_ids,state)
      
      call transform_all_to_s(grid,state)
    else
      call check(ierror)
    end if
    call check(nfmpi_close(nc_ids%ncid))

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

    call transform_all_to_p(grid,state)
  end subroutine

  ! TODO: This fails on Bridges2
  subroutine read_spectral(grid,nc_ids,state)
    !/********************************************************************|
    !| READ_SPECTRAL                                                      |
    !| -------------                                                      |
    !| This subroutine reads in a spectral input file to read the data    |
    !| into state.                                                        |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   nc_ids: the netCDF variable ids associated with the input file   |
    !|   state: the state_type object containing the simulation data      |
    !\********************************************************************/
    use mpi
    use pnetcdf
    use output_module
    use comm_module, only: get_id
    implicit none

    type(grid_type), intent(in) :: grid
    type(nc_ids_type), intent(inout) :: nc_ids
    type(state_type), intent(inout) :: state

    complex(kind=dp),parameter :: iu = (0.0_dp,1.0_dp)
    integer(kind=mpi_offset_kind)  :: offset(4),count(4)
    integer :: itmp(1),ierror,ncid
    real(kind=dp) :: rtmp(1)
    real(kind=dp), allocatable :: rwork(:,:,:)

    ncid = nc_ids%ncid

    offset(4) = 1
    ierror = nfmpi_inq_dimid(ncid,"time",nc_ids%time_dimid)
    if (ierror == nf_noerr) then
      call check(nfmpi_inq_dimlen(ncid,nc_ids%time_dimid,offset(4)))
    else if (ierror /= nf_ebaddim) then
      call check(ierror)
    endif

    offset(1) = 1
    offset(2) = grid%x_s_offset + 1
    offset(3) = grid%y_s_offset + 1

    count(1) = 2 * grid%z_modes
    count(2) = grid%x_s_count
    count(3) = grid%y_s_count
    count(4) = 1

    call check(nfmpi_inq_varid(ncid,"RHO",nc_ids%rho))
    call check(nfmpi_inq_varid(ncid,"ETA",nc_ids%en))
    call check(nfmpi_inq_varid(ncid,"C",nc_ids%c))
    call check(nfmpi_inq_varid(ncid,"U",nc_ids%u))
    call check(nfmpi_inq_varid(ncid,"V",nc_ids%v))
    call check(nfmpi_inq_varid(ncid,"W",nc_ids%w))
    call check(nfmpi_inq_varid(ncid,"RHOI",nc_ids%rhoaux))
    call check(nfmpi_inq_varid(ncid,"ETAI",nc_ids%enaux))
    call check(nfmpi_inq_varid(ncid,"CI",nc_ids%caux))
    call check(nfmpi_inq_varid(ncid,"UI",nc_ids%uaux))
    call check(nfmpi_inq_varid(ncid,"VI",nc_ids%vaux))
    call check(nfmpi_inq_varid(ncid,"WI",nc_ids%waux))

    allocate(rwork(2 * grid%z_modes,grid%x_s_count,grid%y_s_count))
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%rho,offset,count,rwork))
    state%rho%s = rwork
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%rhoaux,offset,count,rwork))
    state%rho%s = state%rho%s + iu * rwork
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%en,offset,count,rwork))
    state%en%s = rwork
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%enaux,offset,count,rwork))
    state%en%s = state%en%s + iu * rwork
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%c,offset,count,rwork))
    state%c%s = rwork
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%caux,offset,count,rwork))
    state%c%s = state%c%s + iu * rwork
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%u,offset,count,rwork))
    state%u%s = rwork
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%uaux,offset,count,rwork))
    state%u%s = state%u%s + iu * rwork
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%v,offset,count,rwork))
    state%v%s = rwork
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%vaux,offset,count,rwork))
    state%v%s = state%v%s + iu * rwork
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%w,offset,count,rwork))
    state%w%s = rwork
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%waux,offset,count,rwork))
    state%w%s = state%w%s + iu * rwork
    deallocate(rwork)
  end subroutine

  subroutine read_data(grid,nc_ids,state)
    !/********************************************************************|
    !| READ_DATA                                                          |
    !| ---------                                                          |
    !| This subroutine reads in a physical input file to read the data    |
    !| into state.                                                        |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   nc_ids: the netCDF variable ids associated with the input file   |
    !|   state: the state_type object containing the simulation data      |
    !\********************************************************************/
    use mpi
    use pnetcdf
    use comm_module, only: get_id
    implicit none

    type(grid_type), intent(in) :: grid
    type(nc_ids_type), intent(inout) :: nc_ids
    type(state_type), intent(inout) :: state

    integer(kind=mpi_offset_kind)  :: offset(4),count(4)
    integer(kind=mpi_offset_kind) :: ltmp
    integer :: itmp(1),ierror,ncid
    real(kind=dp) :: rtmp(1)

    ncid = nc_ids%ncid

    offset(4) = 1
    ierror = nfmpi_inq_dimid(ncid,"time",nc_ids%time_dimid)
    if (ierror == nf_noerr) then
      call check(nfmpi_inq_dimlen(ncid,nc_ids%time_dimid,offset(4)))
    else if (ierror /= nf_ebaddim) then
      call check(ierror)
    endif

    offset(1) = 1
    offset(2) = grid%y_p_offset + 1
    offset(3) = grid%z_p_offset + 1

    count(1) = grid%x_n
    count(2) = grid%y_p_count
    count(3) = grid%z_p_count
    count(4) = 1

    call check(nfmpi_inq_varid(ncid,"RHO",nc_ids%rho))
    call check(nfmpi_inq_varid(ncid,"ETA",nc_ids%en))
    call check(nfmpi_inq_varid(ncid,"C",nc_ids%c))
    call check(nfmpi_inq_varid(ncid,"U",nc_ids%u))
    call check(nfmpi_inq_varid(ncid,"V",nc_ids%v))
    call check(nfmpi_inq_varid(ncid,"W",nc_ids%w))

    call check(nfmpi_get_vara_double_all(ncid,nc_ids%rho,offset,count,state%rho%p))
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%en,offset,count,state%en%p))
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%c,offset,count,state%c%p))
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%u,offset,count,state%u%p))
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%v,offset,count,state%v%p))
    call check(nfmpi_get_vara_double_all(ncid,nc_ids%w,offset,count,state%w%p))
  end subroutine

  !/**********************************************************************|
  !| FUNCTIONS                                                            |
  !\**********************************************************************/
  function idefault(ncid,variable,val)
    !/********************************************************************|
    !| IDEFAULT                                                           |
    !| --------                                                           |
    !| A helper function that tries to read a variable from an open       |
    !| netCDF file. If the variable doesn't exist, a default value (val)  |
    !| is used instead.                                                   |
    !|   ncid: the netCDF file ID associated with the file                | 
    !|   variable: the character string representation of the variable    |
    !|   val: the default value to use on failure                         |
    !\********************************************************************/
    use mpi
    use pnetcdf
    implicit none

    integer :: ncid
    character(*) :: variable
    integer :: val
    integer :: idefault

    integer :: varid,dimid
    integer(kind=mpi_offset_kind) :: offset(1),count(1),step
    integer :: itmp(1),ierror,status

    ierror = nfmpi_inq_dimid(ncid,"time",dimid)

    status = nfmpi_inq_varid(ncid,variable,varid)
    if (status == nf_enotvar) then
      idefault = val
    else if (status == nf_noerr) then
      if (ierror == nf_ebaddim) then
        call check(nfmpi_get_var_int_all(ncid,varid,itmp))
      else
        call check(nfmpi_inq_dimlen(ncid,dimid,step))
  
        offset(1) = step
        count = (/1/)

        call check(nfmpi_get_vara_int_all(ncid,varid,offset,count,itmp))
      endif
      idefault = itmp(1)
    else
      call check(status)
    endif
  end function

  function rdefault(ncid,variable,val)
    !/********************************************************************|
    !| RDEFAULT                                                           |
    !| --------                                                           |
    !| A helper function that tries to read a variable from an open       |
    !| netCDF file. If the variable doesn't exist, a default value (val)  |
    !| is used instead.                                                   |
    !|   ncid: the netCDF file ID associated with the file                | 
    !|   variable: the character string representation of the variable    |
    !|   val: the default value to use on failure                         |
    !\********************************************************************/
    use mpi
    use pnetcdf
    implicit none

    integer :: ncid
    character(*) :: variable
    real(kind=dp) :: val
    real(kind=dp) :: rdefault

    integer :: varid,dimid
    integer(kind=mpi_offset_kind) :: offset(1),count(1),step
    real(kind=dp) :: rtmp(1)
    integer :: ierror,status

    ierror = nfmpi_inq_dimid(ncid,"time",dimid)

    status = nfmpi_inq_varid(ncid,variable,varid)
    if (status == nf_enotvar) then
      rdefault = val
    else if (status == nf_noerr) then
      if (ierror == nf_ebaddim) then
        call check(nfmpi_get_var_double_all(ncid,varid,rtmp))
      else
        call check(nfmpi_inq_dimlen(ncid,dimid,step))
  
        offset(1) = step
        count = (/1/)

        call check(nfmpi_get_vara_double_all(ncid,varid,offset,count,rtmp))
      endif
      rdefault = rtmp(1)
    else
      call check(status)
    endif
  end function
