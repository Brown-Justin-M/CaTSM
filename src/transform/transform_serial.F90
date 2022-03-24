  subroutine init_plans(x_modes,y_modes,z_modes,x_n,y_n,z_n,y_pencils,z_pencils)
    !/********************************************************************|
    !| INIT_PLANS                                                         |
    !| ----------                                                         |
    !| This subroutine uses the information provided about the grid to    |
    !| construct the FFTW plan objects and stores them in a               |
    !| transform_type object.                                             |
    !|   x_modes: the number of Fourier modes in the x direction          |
    !|   y_modes: the number of Fourier modes in the y direction          |
    !|   z_modes: the number of Fourier modes in the z direction          |
    !|   x_n: the number of grid points in the x direction                |
    !|   y_n: the number of grid points in the y direction                |
    !|   z_n: the number of grid points in the z direction                |
    !|   y_pencils: unused (for MPI implementation)                       |
    !|   z_pencils: unused (for MPI implementation)                       |
    !\********************************************************************/
    implicit none
    include "fftw3.f"

    integer, intent(in) :: x_modes,y_modes,z_modes,x_n,y_n,z_n
    integer, optional :: y_pencils,z_pencils

    complex(kind=dp), pointer :: cwork(:,:,:)
    real(kind=dp), pointer :: rwork(:,:,:)
    integer :: howmany ! the total number of transforms
    integer :: nembed ! the size of the individual arrays
    integer :: inembed2(2),onembed2(2) ! the size of the individual arrays
    integer :: stride ! the distance in memory between adjacent elements
    integer :: dist ! the distance in memory between adjacent arrays
    integer :: flags

    transform%x_modes = x_modes
    transform%y_modes = y_modes
    transform%z_modes = z_modes
    transform%x_n = x_n
    transform%z_n = z_n
    transform%x_s_count = x_modes + 1
    transform%y_s_count = max(2 * y_modes,1)
    transform%z_s_count = max(2 * z_modes,1)
    transform%x_p_count = x_n
    transform%y_p_count = y_n
    transform%z_p_count = z_n
    flags = ior(fftw_patient,fftw_preserve_input)

    ! Set up the parameters for the z transform
    howmany = transform%x_s_count * transform%y_s_count
    nembed = z_n
    stride = 1
    dist = z_n

    ! Construct the z transform plans
    allocate(cwork(z_n,transform%x_s_count,transform%y_s_count))
    call dfftw_plan_many_dft(transform%z_s_to_p,1,z_n,howmany, &
      & cwork,nembed,stride,dist,cwork,nembed,stride,dist,fftw_backward,flags)
    call dfftw_plan_many_dft(transform%z_p_to_s,1,z_n,howmany, &
      & cwork,nembed,stride,dist,cwork,nembed,stride,dist,fftw_forward,flags)
    deallocate(cwork)

    ! Set up the parameters for the xy transform
    inembed2 = (/x_n / 2 + 1,y_n/)
    onembed2 = (/x_n,y_n/)
    stride = z_n
    dist = 1
    flags = ior(fftw_patient,fftw_destroy_input)

    ! Construct the xy transform plans
    allocate(rwork(transform%z_p_count,transform%x_n,transform%y_p_count))
    allocate(cwork(transform%z_p_count,transform%x_n / 2 + 1,transform%y_p_count))
    call dfftw_plan_many_dft_c2r(transform%xy_s_to_p,2,(/x_n,y_n/),z_n, &
      & cwork,inembed2,stride,dist,rwork,onembed2,stride,dist,flags)
    call dfftw_plan_many_dft_r2c(transform%xy_p_to_s,2,(/x_n,y_n/),z_n, &
      & rwork,onembed2,stride,dist,cwork,inembed2,stride,dist,flags)
    deallocate(cwork,rwork)
  end subroutine

  subroutine free_plans
    !/********************************************************************|
    !| FREE_PLANS                                                         |
    !| ----------                                                         |
    !| Clean up everything associated with the transforms.                |
    !\********************************************************************/
    implicit none
    include "fftw3.f"
    
    call dfftw_cleanup
  end subroutine

  ! TODO: If the arrays are the same size (almost never), this could be done in place
  subroutine transform_p_to_s(in,out)
    !/********************************************************************|
    !| TRANSFORM_P_TO_S                                                   |
    !| ----------------                                                   |
    !| Given an input in physical space, transform it to spectral space   |
    !| and place the output into out. The arrays are expected to have the |
    !| following dimensions:                                              |
    !|   in(transform%x_p_count,transform%y_p_count,transform%z_p_count)  |
    !|   out(transform%z_s_count,transform%x_s_count,transform%y_s_count) |
    !\********************************************************************/
    implicit none
    real(kind=dp), intent(in) :: in(:,:,:)
    complex(kind=dp), intent(out) :: out(:,:,:)

    complex(kind=dp), allocatable :: cwork(:,:,:)
    integer :: z_modes

    z_modes = transform%z_modes

    allocate(cwork(transform%z_n,transform%x_s_count,transform%y_s_count))

    call transform_xy_p_to_s(in,cwork)

    call dfftw_execute_dft(transform%z_p_to_s,cwork,cwork)

    ! Remove excess dealiased modes and set the Nyquist frequency to 0
    out(:z_modes,:,:) = cwork(:z_modes,:,:)
    out(z_modes + 1,:,:) = 0.0_dp
    out(z_modes + 2:,:,:) = cwork(transform%z_n - z_modes + 2:,:,:)
    deallocate(cwork)

    out = out / transform%z_n
  end subroutine

  ! TODO: If the arrays are the same size (almost never), this could be done in place
  subroutine transform_s_to_p(in,out)
    !/********************************************************************|
    !| TRANSFORM_S_TO_P                                                   |
    !| ----------------                                                   |
    !| Given an input in spectral space, transform it to physical space   |
    !| and place the output into out. The arrays are expected to have the |
    !| following dimensions:                                              |
    !|   in(transform%z_s_count,transform%x_s_count,transform%y_s_count)  |
    !|   out(transform%x_p_count,transform%y_p_count,transform%z_p_count) |
    !\********************************************************************/
    implicit none
    include "fftw3.f"

    complex(kind=dp), intent(in) :: in(:,:,:)
    real(kind=dp), intent(out) :: out(:,:,:)

    complex(kind=dp), allocatable :: cwork(:,:,:)
    integer :: z_modes

    z_modes = transform%z_modes

    allocate(cwork(transform%z_p_count,transform%x_s_count,transform%y_s_count))

    ! Expand the arrays with zero padding for the physical grid
    cwork(:z_modes,:,:) = in(:z_modes,:,:)
    cwork(z_modes + 1:transform%z_p_count - z_modes,:,:) = 0.0_dp
    cwork(transform%z_p_count - z_modes + 1:,:,:) = in(z_modes + 1:,:,:)

    call dfftw_execute_dft(transform%z_s_to_p,cwork,cwork)

    call transform_xy_s_to_p(cwork,out)
    deallocate(cwork)
  end subroutine

  subroutine transform_xy_p_to_s(in,out)
    !/********************************************************************|
    !| TRANSFORM_XY_P_TO_S                                                |
    !| -------------------                                                |
    !| Given an input in physical space, transform it to be spectral in   |
    !| only x and y. This produces an array that's in physical space in   |
    !| z. This is useful for various operations that depend on the z      |
    !| coordinate. Note that the dimensions are reordered. The arrays are |
    !| expected to have the following dimensions:                         |
    !|   in(transform%x_p_count,transform%y_p_count,transform%z_p_count)  |
    !|   out(transform%z_n,transform%x_s_count,transform%y_s_count)       |
    !\********************************************************************/
    implicit none
    real(kind=dp), intent(in) :: in(:,:,:)
    complex(kind=dp), intent(out) :: out(:,:,:)

    integer :: x_modes,y_modes
    complex(kind=dp), allocatable :: cwork(:,:,:)
    real(kind=dp), allocatable :: rwork(:,:,:)
    integer :: i,j,k

    x_modes = transform%x_modes
    y_modes = transform%y_modes
  
    allocate(rwork(transform%z_p_count,transform%x_n,transform%y_p_count))

    ! Transpose the array
    do k = 1,transform%z_p_count
      do j = 1,transform%y_p_count
        do i = 1,transform%x_n
          rwork(k,i,j) = in(i,j,k)
        end do
      end do
    end do

    allocate(cwork(transform%z_p_count,transform%x_n / 2 + 1,transform%y_p_count))
    call dfftw_execute_dft_r2c(transform%xy_p_to_s,rwork,cwork)
    deallocate(rwork)

    ! Remove excess dealiased modes and set the Nyquist frequency to 0
    if (y_modes /= 0) then
      out(:,:,:y_modes) = cwork(:,:x_modes + 1,:y_modes)
      out(:,:,y_modes + 2) = 0.0_dp
      out(:,:,y_modes + 2:) = cwork(:,:x_modes + 1,transform%y_p_count - y_modes + 2:)
    else
      out(:,:,:) = cwork(:,:x_modes + 1,:)
    end if
    out(:,x_modes + 1,:) = 0.0_dp
    deallocate(cwork)

    out = out / (transform%x_n * transform%y_p_count)
  end subroutine

  subroutine transform_xy_s_to_p(in,out)
    !/********************************************************************|
    !| TRANSFORM_XY_S_TO_P                                                |
    !| -------------------                                                |
    !| Given an input in physical space in z but spectral space in x and  |
    !| y, transform it to be in physical space. This is useful for        |
    !| various operations that depend on the z coordinate. Note that the  |
    !| dimensions are reordered. The arrays are expected to have the      |
    !| following dimensions:                                              |
    !|   in(transform%z_n,transform%x_s_count,transform%y_s_count)        |
    !|   out(transform%x_p_count,transform%y_p_count,transform%z_p_count) |
    !\********************************************************************/
    implicit none
    complex(kind=dp), intent(in) :: in(:,:,:)
    real(kind=dp), intent(out) :: out(:,:,:)
  
    integer :: x_modes,y_modes
    complex(kind=dp), allocatable :: cwork(:,:,:)
    real(kind=dp), allocatable :: rwork(:,:,:)
    integer :: i,j,k

    x_modes = transform%x_modes
    y_modes = transform%y_modes

    allocate(cwork(transform%z_p_count,transform%x_n / 2 + 1,transform%y_p_count))

    ! Expand the arrays with zero padding for the physical grid
    cwork(:,x_modes + 2:,:) = 0.0_dp
    if (y_modes /= 0) then
      cwork(:,:x_modes + 1,:y_modes) = in(:,:,:y_modes)
      cwork(:,:,y_modes + 1:transform%y_p_count - y_modes) = 0.0_dp
      cwork(:,:x_modes + 1,transform%y_p_count - y_modes + 1:) = in(:,:,y_modes + 1:)
    else
      cwork(:,:x_modes + 1,:) = in(:,:,:)
    end if

    allocate(rwork(transform%z_p_count,transform%x_p_count,transform%y_p_count))
    call dfftw_execute_dft_c2r(transform%xy_s_to_p,cwork,rwork)
    deallocate(cwork)

    ! Transpose the array
    do k = 1,transform%z_p_count
      do j = 1,transform%y_p_count
        do i = 1,transform%x_p_count
          out(i,j,k) = rwork(k,i,j)
        end do
      end do
    end do
    deallocate(rwork)
  end subroutine

  subroutine transform_z_p_to_s(in,out)
    !/********************************************************************|
    !| TRANSFORM_Z_P_TO_S                                                 |
    !| ------------------                                                 |
    !| Given an input in physical space in z and spectral space in x and  |
    !| y, transform it to spectral space. The arrays are expected to have |
    !| the following dimensions:                                          |
    !|   in(transform%z_n,transform%x_s_count,transform%y_s_count)        |
    !|   out(transform%z_s_count,transform%x_s_count,transform%y_s_count) |
    !\********************************************************************/
    implicit none
    complex(kind=dp), intent(in) :: in(:,:,:)
    complex(kind=dp), intent(out) :: out(:,:,:)

    integer :: z_modes
    complex(kind=dp), allocatable :: cwork(:,:,:)

    z_modes = transform%z_modes

    allocate(cwork(transform%z_n,transform%x_s_count,transform%y_s_count))

    call dfftw_execute_dft(transform%z_p_to_s,in,cwork)

    ! Remove excess dealiased modes and set the Nyquist frequency to 0
    out(:z_modes,:,:) = cwork(:z_modes,:,:)
    out(z_modes + 1,:,:) = 0.0_dp
    out(z_modes + 2:,:,:) = cwork(transform%z_n - z_modes + 2:,:,:)
    deallocate(cwork)

    out = out / transform%z_n
  end subroutine

  subroutine transform_z_s_to_p(in,out)
    !/********************************************************************|
    !| TRANSFORM_Z_S_TO_P                                                 |
    !| ------------------                                                 |
    !| Given an input in spectral space, transform it to physical space   |
    !| in z and spectral space in x and y. The arrays are expected to     |
    !| have the following dimensions:                                     |
    !|   in(transform%z_s_count,transform%x_s_count,transform%y_s_count)  |
    !|   out(transform%z_n,transform%x_s_count,transform%y_s_count)       |
    !\********************************************************************/
    implicit none
    complex(kind=dp) :: in(:,:,:)
    complex(kind=dp) :: out(:,:,:)

    integer :: z_modes
    complex(kind=dp), allocatable :: cwork(:,:,:)

    z_modes = transform%z_modes

    allocate(cwork(transform%z_n,transform%x_s_count,transform%y_s_count))

    ! Expand the arrays with zero padding for the physical grid
    cwork(:z_modes,:,:) = in(:z_modes,:,:)
    cwork(z_modes + 1:transform%z_n - z_modes,:,:) = 0.0_dp
    cwork(transform%z_n - z_modes + 1:,:,:) = in(z_modes + 1:,:,:)

    call dfftw_execute_dft(transform%z_s_to_p,cwork,out)
    deallocate(cwork)
  end subroutine