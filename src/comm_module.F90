module comm_module
  !/**********************************************************************|
  !| COMM MODULE                                                          |
  !| -----------                                                          |
  !| This module contains the subroutines and variables used for inter-   |
  !| process communication. If the parallel transform library is not      |
  !| available, it provides routines that work in the absence of MPI.     |
  !\**********************************************************************/
  implicit none
  save
  private
  public get_id,get_tasks,get_comm,get_wtime
  public init_comm,end_comm,abort_comm
  public collect_max,collect_array_max,collect_sum

  integer :: id
  integer :: tasks
  integer :: comm

contains
  !/**********************************************************************|
  !| FUNCTIONS                                                            |
  !\**********************************************************************/
  function get_id()
    !/********************************************************************|
    !| GET_ID                                                             |
    !| ------                                                             |
    !| A getter function that returns the current global process ID of    |
    !| the task. This permits other modules to read this variable without |
    !| providing write access.                                            |
    !|   returns: the current process id (0 for serial code)              |
    !\********************************************************************/
    integer :: get_id

    get_id = id
  end function

  function get_tasks()
    !/********************************************************************|
    !| GET_TASKS                                                          |
    !| ---------                                                          |
    !| A getter function that returns the current number of tasks running |
    !| This permits other modules to read this variable without providing |
    !| write access.                                                      |
    !|   returns: the number of running tasks (1 for serial code)         |
    !\********************************************************************/
    integer :: get_tasks

    get_tasks = tasks
  end function

  function get_comm()
    !/********************************************************************|
    !| GET_COMM                                                           |
    !| --------                                                           |
    !| A getter function that returns the current global MPI communicator.|
    !|   returns: the MPI communicator (0 for serial code)                |
    !\********************************************************************/
    integer :: get_comm

    get_comm = comm
  end function

  function get_wtime()
    !/********************************************************************|
    !| GET_WTIME                                                          |
    !| ---------                                                          |
    !| This function calculates the current time for profiling.           |
    !|   returns: the current time as a floating point number             |
    !\********************************************************************/
    real :: get_wtime

    call cpu_time(get_wtime)
  end function

  function collect_max(in)
    !/********************************************************************|
    !| COLLECT_MAX                                                        |
    !| -----------                                                        |
    !| Calculate the max of in across the processes and return the max to |
    !| the root process. The output of this function is undefined for     |
    !| other processes.                                                   |
    !|   in: the floating point number from each process to compare       |
    !|   returns: the maximum of in from all processes                    |
    !\********************************************************************/
    implicit none

    real(kind=8), intent(in) :: in
    real(kind=8) :: collect_max

    integer :: ierror

    collect_max = in
  end function

  function collect_array_max(in)
    !/********************************************************************|
    !| COLLECT_ARRAY_MAX                                                  |
    !| -----------------                                                  |
    !| Calculate the max of in across the processes and return the max to |
    !| the root process. The output of this function is undefined for     |
    !| other processes.                                                   |
    !|   in: the array from each process from which the max is calculated |
    !|   returns: the maximum of in from all processes                    |
    !\********************************************************************/
    implicit none

    real(kind=8), intent(in) :: in(:,:,:)
    real(kind=8) :: collect_array_max

    real(kind=8) :: local

    local = maxval(in)
    collect_array_max = collect_max(local)
  end function

  function collect_sum(in)
    !/********************************************************************|
    !| COLLECT_SUM                                                        |
    !| -----------                                                        |
    !| Calculate the sum of in across the processes and return the sum to |
    !| the root process. The output of this function is undefined for     |
    !| other processes.                                                   |
    !|   in: the floating point number from each process to sum           |
    !|   returns: the sum of in from all processes                        |
    !\********************************************************************/
    implicit none

    real(kind=8), intent(in) :: in
    real(kind=8) :: collect_sum

    integer :: ierror

    collect_sum = in
  end function

  !/**********************************************************************|
  !| SUBROUTINES                                                          |
  !\**********************************************************************/
  subroutine init_comm
    !/********************************************************************|
    !| INIT_COMM                                                          |
    !| ---------                                                          |
    !| Initialize message passing if MPI is available; otherwise, set the |
    !| communicator variables to their serial defaults. This will raise   |
    !| an error if something goes wrong while setting up MPI.             |
    !\********************************************************************/
    implicit none

    integer :: ierror

    comm = 0
    id = 0
    tasks = 1
  end subroutine

  subroutine end_comm
    !/********************************************************************|
    !| END_COMM                                                           |
    !| --------                                                           |
    !| Safely finalize MPI communication.                                 |
    !\********************************************************************/
    implicit none

    integer :: ierror
  end subroutine

  subroutine abort_comm(errcode)
    !/********************************************************************|
    !| ABORT_COMM                                                         |
    !| ----------                                                         |
    !| If an unrecoverable error has arisen, kill all processes and send  |
    !| the appropriate errorcode to the user.                             |
    !|   errcode: the errorcode to send                                   |
    !\********************************************************************/
    implicit none

    integer, optional :: errcode

    integer :: errcodein,ierror

    errcodein = 1
    if (present(errcode)) errcodein = errcode
    stop
  end subroutine

end module comm_module
