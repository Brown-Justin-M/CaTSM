module fields_module
  !/**********************************************************************|
  !| FIELDS MODULE                                                        |
  !| -------------                                                        |
  !| This module contains the state_type and rhs_type derived type, which |
  !| contain the state information and records of the simulation. It also |
  !| Includes a number of utilities for copying, updating, and setting up |
  !| the fields contained in these objects.                               |
  !\**********************************************************************/
  use grid_module, only: dp,li,grid_type
  implicit none
  private
  public field_type,vector_type,state_type,rhs_type
  public sum_squares,sum_squares_vector,sum_product_p
  public init_state,free_state,init_rhs,free_rhs,init_records,free_records
  public step_time,copy_state
  public transform_all_to_p,transform_all_to_s

  !/**********************************************************************|
  !| DERIVED TYPES                                                        |
  !\**********************************************************************/
  type field_type
    !/********************************************************************|
    !| FIELD_TYPE                                                         |
    !| ----------                                                         |
    !| The field type contains the pointers for a field's arrays, with    |
    !| p pointer pointing towards the field in physical space, and the s  |
    !| pointer pointer towards the field in spectral space. This is a     |
    !| helper type typically contained within a state object.             |
    !\********************************************************************/
    complex(kind=dp), pointer :: s(:,:,:) ! the pointer to the spectral array
    real(kind=dp), pointer :: p(:,:,:) ! the pointer to the physical array
  end type field_type

  type vector_type
    !/********************************************************************|
    !| VECTOR_TYPE                                                        |
    !| -----------                                                        |
    !| The vector type contains the pointers for a vector's arrays, with  |
    !| p pointer pointing towards the vector in physical space, and the s |
    !| pointer pointer towards the vector in spectral space. This is a    |
    !| helper type typically contained within a state object.             |
    !\********************************************************************/
    complex(kind=dp), pointer :: s(:,:,:,:) ! the pointer to the spectral array
    real(kind=dp), pointer :: p(:,:,:,:) ! the pointer to the physical array
  end type vector_type

  type state_type
    !/********************************************************************|
    !| STATE_TYPE                                                         |
    !| ----------                                                         |
    !| The state type contains everything pertaining to a single step.    |
    !| This derived type should be set up using the init_state            |
    !| subroutine. The number of records can be specified before this     |
    !| subroutine is called but should not be changed afterwards. The     |
    !| state type should be freed with free_state.                        |
    !\********************************************************************/
    complex(kind=dp), pointer :: all_s(:,:,:,:) ! all major variables in spectral space
    real(kind=dp), pointer :: all_p(:,:,:,:) ! all major variables in physical space

    type(vector_type) :: vel ! the velocity vector
    type(vector_type) :: rhovel ! the momentum vector
    type(field_type) :: u ! the x-component of the velocity vector
    type(field_type) :: v ! the y-component of the velocity vector
    type(field_type) :: w ! the z-component of the velocity vector
    type(field_type) :: rhou ! the x-component of the momentum vector
    type(field_type) :: rhov ! the y-component of the momentum vector
    type(field_type) :: rhow ! the z-component of the momentum vector

    type(field_type) :: t ! the temperature field
    type(field_type) :: c ! the concentration field

    type(field_type) :: en ! the entropy field
    type(field_type) :: rho ! the density field
    type(field_type) :: mu ! the potential temperature field
    type(field_type) :: p ! the pressure field
    type(field_type) :: cp ! the specific heat field

    complex(kind=dp), pointer :: init_en(:)
    logical :: init_en_setup = .true.
    complex(kind=dp), pointer :: init_c(:)
    logical :: init_c_setup = .true.

    integer :: timestep ! the integer timestep
    real(kind=dp) :: time ! the floating point time
    real(kind=dp) :: dt = 0.0_dp ! the timestep duration
    integer :: restart = 0 ! the last recorded restart timestep
    integer :: n_records = 3 ! the number of records to store
    integer :: n_vars ! the number of variables in the state
  end type state_type

  type rhs_type
    !/********************************************************************|
    !| RHS_TYPE                                                           |
    !| --------                                                           |
    !| The rhs type contains the various right hand sides of the          |
    !| equations for the variables in the state type. These are all in    |
    !| spectral space. This derived type should be constructed using the  |
    !| init_rhs subroutine and freed with free_rhs.                       |
    !\********************************************************************/
    complex(kind=dp), pointer :: all(:,:,:,:) ! the array of all RHS
    
    complex(kind=dp), pointer :: vel(:,:,:,:) ! the velocity RHS
    complex(kind=dp), pointer :: rhovel(:,:,:,:) ! the momentum RHS
    complex(kind=dp), pointer :: u(:,:,:) ! the x-component of the velocity RHS
    complex(kind=dp), pointer :: v(:,:,:) ! the y-component of the velocity RHS
    complex(kind=dp), pointer :: w(:,:,:) ! the z-component of the velocity RHS
    complex(kind=dp), pointer :: rhou(:,:,:) ! the x-component of the momentum RHS
    complex(kind=dp), pointer :: rhov(:,:,:) ! the y-component of the momentum RHS
    complex(kind=dp), pointer :: rhow(:,:,:) ! the z-component of the momentum RHS

    complex(kind=dp), pointer :: en(:,:,:) ! the entropy RHS
    complex(kind=dp), pointer :: rho(:,:,:) ! the density RHS
    complex(kind=dp), pointer :: c(:,:,:) ! the concentration RHS
    complex(kind=dp), pointer :: t(:,:,:) ! the temperature RHS
    complex(kind=dp), pointer :: mu(:,:,:) ! the chemical potential RHS
    complex(kind=dp), pointer :: cp(:,:,:) ! the specific heat RHS
    complex(kind=dp), pointer :: p(:,:,:) ! the pressure RHS

    integer :: n_vars ! the number of variables recorded
  end type rhs_type

contains

  !/**********************************************************************|
  !| FUNCTIONS                                                            |
  !\**********************************************************************/
  function sum_squares(grid,in)
    !/********************************************************************|
    !| SUM_SQUARES                                                        |
    !| -----------                                                        |
    !| Calculate the sum of the squares for in, which should be in        |
    !| spectral space. Note that this is scaled to be equal to the volume |
    !| mean of the relevant physical field.                               |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   in: a 3D array to be averaged in spectral space                  |
    !|   returns: the volume average of the field/sum of the squares      |
    !\********************************************************************/
    use comm_module, only: get_comm,collect_sum
    implicit none

    type(grid_type), intent(in) :: grid
    complex(kind=dp), intent(in) :: in(:,:,:)
    real(kind=dp) :: sum_squares

    integer ::  ierror,i,j

    sum_squares = 0.0_dp
    do j = 1,grid%y_s_count
      do i = 1,grid%x_s_count
        ! Because the spectral representation only contains half the modes
        ! non-zero modes must be doubled (due to rfft symmetry)
        if (grid%lkx(i) == 0.0_dp) then
          sum_squares = sum_squares + 1.0_dp * sum(abs(in(:,i,j))**2)
        else
          sum_squares = sum_squares + 2.0_dp * sum(abs(in(:,i,j))**2)
        end if
      end do
    end do

    sum_squares = collect_sum(sum_squares)
  end function

  function sum_squares_vector(grid,in)
    !/********************************************************************|
    !| SUM_SQUARES_VECTOR                                                 |
    !| -----------                                                        |
    !| Calculate the sum of the squares for in, which should be in        |
    !| spectral space. Note that this is scaled to be equal to the volume |
    !| mean of the relevant physical field.                               |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   in: a 4D array to be averaged in spectral space                  |
    !|   returns: the volume average of the field/sum of the squares      |
    !\********************************************************************/
    use comm_module, only: get_comm,collect_sum
    implicit none

    type(grid_type), intent(in) :: grid
    complex(kind=dp), intent(in) :: in(:,:,:,:)
    real(kind=dp) :: sum_squares_vector
    
    real(kind=dp) :: local
    integer :: ierror,i,j,l

    local = 0.0_dp
    do l = 1,3
      do j = 1,grid%y_s_count
        do i = 1,grid%x_s_count
          ! Because the spectral representation only contains half the modes
          ! non-zero modes must be doubled (due to rfft symmetry)
          if (grid%lkx(i) == 0.0_dp) then
            local = local + 1.0_dp * sum(abs(in(:,i,j,l))**2)
          else
            local = local + 2.0_dp * sum(abs(in(:,i,j,l))**2)
          end if
        end do
      end do
    end do

    sum_squares_vector = collect_sum(local)
  end function

  function sum_product_p(grid,a,b)
    !/********************************************************************|
    !| SUM_PRODUCT_P                                                      |
    !| -------------                                                      |
    !| Calculate the sum of the product of two quantities a and b, which  |
    !| both be in spectral space.                                         |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   a: a 3D array to be multiplied by and summed in physical space   |
    !|   b: a 3D array to be multiplied by and summed in physical space   |
    !|   returns: the sum of the product a * b                            |
    !\********************************************************************/
    use comm_module, only: get_comm,collect_sum
    implicit none

    type(grid_type), intent(in) :: grid
    real(kind=dp), intent(in) :: a(:,:,:)
    real(kind=dp), intent(in) :: b(:,:,:)
    real(kind=dp) :: sum_product_p
    
    real(kind=dp) :: local
    integer :: ierror,k,l

    local = sum(a * b)

    sum_product_p = collect_sum(local)
  end function

  !/**********************************************************************|
  !| SUBROUTINES                                                          |
  !\**********************************************************************/
  subroutine init_state(grid,state)
    !/********************************************************************|
    !| INIT_STATE                                                         |
    !| ----------                                                         |
    !| Initialize a state_type object. This involves allocating the all_p |
    !| and all_s arrays and assigning the various fields and vectors      |
    !| within the state_type object to sub arrays of those arrays.        |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   state: a state_type object containing the simulation data        |
    !\********************************************************************/
    use grid_module, only: allocate_p,allocate_s
    use grid_module, only: allocate_p_vector,allocate_s_vector
    implicit none

    type(grid_type), intent(in) :: grid
    type(state_type), intent(inout) :: state

    state%n_vars = 13

    ! Allocate the variables
    call allocate_p_vector(grid,state%all_p,state%n_vars)
    call allocate_s_vector(grid,state%all_s,state%n_vars)
    allocate (state%init_en(grid%z_n))
    allocate (state%init_c(grid%z_n))

    ! Point the main variables to their locations in all_s/p
    state%vel%s => state%all_s(:,:,:,1:3)
    state%rhovel%s => state%all_s(:,:,:,4:6)
    state%t%s => state%all_s(:,:,:,7)
    state%c%s => state%all_s(:,:,:,8)
    state%en%s => state%all_s(:,:,:,9)
    state%mu%s => state%all_s(:,:,:,10)
    state%p%s => state%all_s(:,:,:,11)
    state%rho%s => state%all_s(:,:,:,12)
    state%cp%s => state%all_s(:,:,:,13)

    state%vel%p => state%all_p(:,:,:,1:3)
    state%rhovel%p => state%all_p(:,:,:,4:6)
    state%t%p => state%all_p(:,:,:,7)
    state%c%p => state%all_p(:,:,:,8)
    state%en%p => state%all_p(:,:,:,9)
    state%mu%p => state%all_p(:,:,:,10)
    state%p%p => state%all_p(:,:,:,11)
    state%rho%p => state%all_p(:,:,:,12)
    state%cp%p => state%all_p(:,:,:,13)

    ! Point the vector components to their locations in the vector arrays
    state%u%s => state%vel%s(:,:,:,1)
    state%v%s => state%vel%s(:,:,:,2)
    state%w%s => state%vel%s(:,:,:,3)

    state%u%p => state%vel%p(:,:,:,1)
    state%v%p => state%vel%p(:,:,:,2)
    state%w%p => state%vel%p(:,:,:,3)

    ! Point the vector components to their locations in the vector arrays
    state%rhou%s => state%rhovel%s(:,:,:,1)
    state%rhov%s => state%rhovel%s(:,:,:,2)
    state%rhow%s => state%rhovel%s(:,:,:,3)

    state%rhou%p => state%rhovel%p(:,:,:,1)
    state%rhov%p => state%rhovel%p(:,:,:,2)
    state%rhow%p => state%rhovel%p(:,:,:,3)

    ! Initialize all fields to 0
    state%all_s = 0.0_dp
    state%all_p = 0.0_dp
  end subroutine

  subroutine free_state(state)
    !/********************************************************************|
    !| FREE_STATE                                                         |
    !| ----------                                                         |
    !| Free a state_type object by deallocating all the relevant arrays.  | 
    !|   state: a state_type object containing the simulation data        |
    !\********************************************************************/
    implicit none

    type(state_type), intent(inout) :: state

    deallocate(state%all_p,state%all_s,state%init_en,state%init_c)
  end subroutine

  subroutine init_rhs(grid,rhs)
    !/********************************************************************|
    !| INIT_RHS                                                           |
    !| --------                                                           |
    !| Initialize an rhs_type object. This involves allocating the all    |
    !| array and assigning the various fields and vectors within the      |
    !| object to sub arrays of that array.                                |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   rhs: an rhs_type object containing the RHS of the simulation     |
    !\********************************************************************/
    use grid_module, only: allocate_s_vector
    implicit none

    type(grid_type), intent(in) :: grid
    type(rhs_type), intent(out) :: rhs

    rhs%n_vars = 13

    call allocate_s_vector(grid,rhs%all,rhs%n_vars)

    ! Point the variables to their locations inside all
    rhs%vel => rhs%all(:,:,:,1:3)
    rhs%rhovel => rhs%all(:,:,:,4:6)
    rhs%t => rhs%all(:,:,:,7)
    rhs%c => rhs%all(:,:,:,8)
    rhs%en => rhs%all(:,:,:,9)
    rhs%mu => rhs%all(:,:,:,10)
    rhs%p => rhs%all(:,:,:,11)
    rhs%rho => rhs%all(:,:,:,12)
    rhs%cp => rhs%all(:,:,:,13)

    rhs%u => rhs%vel(:,:,:,1)
    rhs%v => rhs%vel(:,:,:,2)
    rhs%w => rhs%vel(:,:,:,3)

    rhs%rhou => rhs%rhovel(:,:,:,1)
    rhs%rhov => rhs%rhovel(:,:,:,2)
    rhs%rhow => rhs%rhovel(:,:,:,3)

    ! Initialize the RHS to 0
    rhs%all = (0.0_dp,0.0_dp)
  end subroutine

  subroutine free_rhs(rhs)
    !/********************************************************************|
    !| FREE_RHS                                                           |
    !| --------                                                           |
    !| Free an rhs_type object by deallocating all the relevant arrays.   | 
    !|   rhs: an rhs_type object containing the RHS of the simulation     |
    !\********************************************************************/
    implicit none

    type(rhs_type), intent(inout) :: rhs

    deallocate(rhs%all)
  end subroutine

  subroutine allocate_records(grid,state,prior,rhs)
    !/********************************************************************|
    !| ALLOCATE_RECORDS                                                   |
    !| ----------------                                                   |
    !| This subroutine allocates all the necessary state and rhs records  |
    !| for a simulation. It allocates a number of prior states and rhs    |
    !| history equal to state%n_records.                                  |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   state: a state_type object containing the current state          |
    !|   prior: an array of state_type objects containing previous states |
    !|   rhs: an array of rhs_type objects containing previous RHS        |
    !\********************************************************************/
    implicit none

    type(grid_type), intent(in) :: grid
    type(state_type), intent(inout) :: state
    type(state_type), pointer, intent(inout) :: prior(:)
    type(rhs_type), pointer, intent(inout) :: rhs(:)

    integer :: i

    call init_state(grid,state)

    allocate(prior(state%n_records))
    allocate(rhs(state%n_records))
    do i = 1,state%n_records
      call init_state(grid,prior(i))
      call init_rhs(grid,rhs(i))
    end do
  end subroutine

  subroutine free_records(state,prior,rhs)
    !/********************************************************************|
    !| FREE_RECORDS                                                       |
    !| ----------------                                                   |
    !| Free all the record objects by calling the appropriate free        |
    !| routines.                                                          | 
    !|   state: a state_type object containing the current state          |
    !|   prior: an array of state_type objects containing previous states |
    !|   rhs: an array of rhs_type objects containing previous RHS        |
    !\********************************************************************/
    implicit none

    type(state_type), intent(inout) :: state
    type(state_type), pointer, intent(inout) :: prior(:)
    type(rhs_type), pointer, intent(inout) :: rhs(:)

    integer :: i

    call free_state(state)
    do i = 1,state%n_records
      call free_state(prior(i))
      call free_rhs(rhs(i))
    end do
    deallocate(prior)
    deallocate(rhs)
  end subroutine

  subroutine copy_state(in,out)
    !/********************************************************************|
    !| COPY_STATE                                                         |
    !| ----------                                                         |
    !| Copy the contents of one state_type object into another existing   |
    !| state_type object. This is a deep copy operation.                  |
    !|   in: the state_type object containing the state to copy from      |
    !|   out: the state_type object containing the state to copy to       |
    !\********************************************************************/
    implicit none

    type(state_type), intent(in) :: in
    type(state_type), intent(out) :: out

    out%all_p = in%all_p
    out%all_s = in%all_s

    out%timestep = in%timestep
    out%time = in%time
    out%dt = in%dt

    out%init_en = in%init_en
    out%init_c = in%init_c
  end subroutine

  subroutine step_time(state,prior,rhs)
    !/********************************************************************|
    !| STEP_TIME                                                          |
    !| ---------                                                          |
    !| This performs the necessary, backend shuffle associated with       |
    !| stepping forward in time. It cycles the arrays in prior and rhs    |
    !| and copies the current state to the most recent prior state. This  |
    !| is done to limit the number of expensive deep copy operations.     |
    !|   state: a state_type object containing the current state          |
    !|   prior: an array of state_type objects containing previous states |
    !|   rhs: an array of rhs_type objects containing previous RHS        |
    !\********************************************************************/
    implicit none

    type(state_type), intent(in) :: state
    type(state_type), intent(inout) :: prior(:)
    type(rhs_type), intent(inout) :: rhs(:)
    
    prior = cshift(prior,-1)
    rhs = cshift(rhs,-1)
    call copy_state(state,prior(1))
  end subroutine

  subroutine init_records(grid,state,prior,rhs)
    !/********************************************************************|
    !| INIT_RECORDS                                                       |
    !| ------------                                                       |
    !| This subroutine initializes the records for the simulation,        |
    !| calling allocate_records as it runs. It uses random numbers to     |
    !| generate the initial fields and ensures that the initial velocity  |
    !| field is divergence free and consistent with the curl. It also     |
    !| converts the initial conditions to spectral space such that the    |
    !| resulting initial fields are entirely self-consistent.             |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   state: a state_type object containing the current state          |
    !|   prior: an array of state_type objects containing previous states |
    !|   rhs: an array of rhs_type objects containing previous RHS        |
    !\********************************************************************/
    use grid_ops_module, only: curl,remove_div
    use transform_module, only: init_plans,transform_s_to_p
    implicit none

    type(grid_type), intent(in) :: grid
    type(state_type), intent(inout) :: state
    type(state_type), pointer, intent(inout) :: prior(:)
    type(rhs_type), pointer, intent(inout) :: rhs(:)

    integer :: ierror

    call allocate_records(grid,state,prior,rhs)

    ! Initialize the fields with random perturbations
    state%vel%p = 0.0_dp
    call init_rn_p(grid,state%t)
    call init_rn_p(grid,state%c)

    ! Initialize the timestep to 0
    state%timestep = 0
    state%restart = 0
    state%time = 0.0_dp
    state%dt = 0.0_dp

    ! Update the rest of the fields to be self consistent
    call transform_all_to_s(grid,state)

    call remove_div(grid,state%time,state%vel%s,state%vel%s)

    call transform_all_to_p(grid,state)
  end subroutine

  subroutine transform_all_to_s(grid,state)
    !/********************************************************************|
    !| TRANSFORM_ALL_TO_S                                                 |
    !| ------------------                                                 |
    !| Transform every field to spectral space. This also sets the mean   |
    !| velocity field to be that of the background flow (typically 0).    |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   state: a state_type object containing the current state          |
    !\********************************************************************/
    use comm_module, only: get_id
    use transform_module, only: transform_p_to_s,transform_s_to_p
    implicit none

    type(grid_type), intent(in) :: grid
    type(state_type), intent(inout) :: state

    integer :: l
    
    do l = 1,state%n_vars
      call transform_p_to_s(state%all_p(:,:,:,l),state%all_s(:,:,:,l))
    end do
  end subroutine

  subroutine transform_all_to_p(grid,state)
    !/********************************************************************|
    !| TRANSFORM_ALL_TO_P                                                 |
    !| ------------------                                                 |
    !| Transform every field to physical space. This also sets the mean   |
    !| velocity field to be that of the background flow (typically 0).    |
    !| The curl is also updated at this stage.                            |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   state: a state_type object containing the current state          |
    !\********************************************************************/
    use grid_ops_module, only: curl,deriv_x,deriv_z
    use transform_module, only: transform_s_to_p
    use comm_module, only: get_id
    implicit none
  
    type(grid_type), intent(in) :: grid
    type(state_type), intent(inout) :: state
  
    integer :: l
  
    do l = 1,state%n_vars
      call transform_s_to_p(state%all_s(:,:,:,l),state%all_p(:,:,:,l))
    end do
  end subroutine

  subroutine init_rn_p(grid,x)
    !/********************************************************************|
    !| INIT_RN_P                                                          |
    !| ---------                                                          |
    !| Initialize a field in physical space using random numbers between  |
    !| -5.e-4 and 5e.-4.                                                  |
    !|   grid: a grid_type object that describes the grid                 | 
    !|   x: the field_type object to initialize                           |
    !\********************************************************************/
    use comm_module, only: get_id
    implicit none

    type(grid_type), intent(in) :: grid
    type(field_type), intent(out) :: x

    integer :: n = 33,seed(33) = 0
    save seed

    ! Ensure that the seed is different on each process but consistent from
    ! simulation to simulation
    seed = seed + 1 + get_id()

    call random_seed(n)
    call random_seed(put=seed)
    call random_number(x%p)
    x%p = (x%p - 0.5_dp) * 1.e-3
  end subroutine

endmodule fields_module
