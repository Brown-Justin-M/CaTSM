!/**********************************************************************|
!|            CCCCCC             TTTTTTTT   SSSSSS   MM    MM           |
!|           CCC  CCC   aaaaaa      TT     SSS  SSS  MMM  MMM           |
!|           CCC            aaa     TT      SSS      MMMMMMMM           |
!|           CCC        aaaaaaa     TT         SSS   MM MM MM           |
!|           CCC  CCC  aaa  aaa     TT     SSS  SSS  MM MM MM           |
!|            CCCCCC    aaaaa a     TT      SSSSSS   MM MM MM           |
!|           ================================================           |
!|           ================================================           |
!|                                                                      |
!| Written by: Justin Brown                                             |
!| Contact: jmbrown2@nps.edu                                            |
!|                                                                      |
!| CaTSM is a fluid code designed to solve double-diffusive and         |
!| microscale fluid problems using compressible physics. The equation   |
!| state permits linear and ideal gas forms. The code uses the          |
!| pseudo-spectral method and is periodic horizontally. Initial steps   |
!| are carried out using a Runge--Kutta/Crank--Nicolson method, but the |
!| bulk of the code uses the Adams--Bashforth/Backwards Differentiation |
!|  Formula.                                                            |
!\**********************************************************************/
program rome
  use comm_module, only: get_id,get_wtime,init_comm,end_comm
  use grid_module, only: dp,grid_type,init_grid,free_grid
  use params_module, only: params_type,read_params
  use fields_module, only: state_type,rhs_type,init_records,free_records
  use schemes_module, only: update
  use equations_module, only: update_eos
  use output_module, only: controllers_type
  use output_module, only: open_files,close_files,write_output_files
  use input_module, only: open_inputs,read_input
  implicit none

  ! These are the main global variables for the calculation
  type(params_type) :: params ! a container for the simulation parameters
  type(grid_type) :: grid ! a container for the grid information
  type(controllers_type) :: controls ! a container for the IO preferences

  type(state_type) :: state ! a container for the arrays of the simulation
  type(state_type), pointer :: prior(:) ! the prior state containers from previous times
  type(rhs_type), pointer :: rhs(:) ! the current and prior RHS for each equation

  real(kind=dp) :: rtmp = 0.0_dp
  integer :: i

  ! Initialize MPI, grid, and record data; read parameters from file
  call init_comm

  call read_params(grid,params)
  call init_grid(grid)

  call open_inputs(grid,params,controls)

  call init_records(grid,state,prior,rhs)
  ! If an input file was specified, read it into the record data
  if (controls%file_input /= "") then
    call read_input(grid,controls%file_input,state)
    call update_eos(grid,params,state)
  endif

  call open_files(grid,params,controls)

  if (params%max_steps .le. 0) then 
    print*,"ERROR: the number of steps is less than 1. Exiting."
    call end_comm
  end if

  ! Write the first outputs using the initial data
  call write_output_files(grid,params,controls,state)

  ! Main simulation loop
  do i = 1,params%max_steps
    if (state%time .ge. params%max_time) then
      if (get_id() == 0) then
        print*, "Reached max_time. Ending."
      end if
      exit
    end if

    rtmp = get_wtime()

    call update(grid,params,state,prior,rhs)

    rtmp = get_wtime() - rtmp
    if (get_id() == 0) then
      print*,state%timestep,", dt :",state%dt,", wallclock:",rtmp,", CS%:",state%dt/grid%dz*params%cs
    end if

    call write_output_files(grid,params,controls,state)
  end do

  ! Clean up opened files, allocated memory, and MPI
  call close_files(controls)

  call free_records(state,prior,rhs)
  call free_grid(grid)

  call end_comm
end program rome
