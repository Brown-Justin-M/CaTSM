module transform_module
  !/**********************************************************************|
  !| TRANSFORM MODULE                                                     |
  !| ----------------                                                     |
  !| This module contains the interface for the transformation to and     |
  !| from physical space for the variables contained within a state       |
  !| object.                                                              |
  !\**********************************************************************/
  implicit none

  !/**********************************************************************|
  !| PARAMETERS                                                           |
  !\**********************************************************************/
  integer, parameter :: dp = selected_real_kind(14)
  integer, parameter :: li = selected_int_kind(18)

  !/**********************************************************************|
  !| DERIVED TYPES                                                        |
  !\**********************************************************************/
  type transform_type
    !/********************************************************************|
    !| TRANSFORM_TYPE                                                     |
    !| --------------                                                     |
    !| This derived type contains the interface to the various fftw plans |
    !| used during the simulation. It also contains the sizes of the      |
    !| arrays for convenience. This type is constructed and handled by    |
    !| init_plans and should not be accessed by the user.                 |
    !\********************************************************************/
    integer(kind=li) :: x_p_to_s,x_s_to_p ! the FFTW plans for x transforms
    integer(kind=li) :: y_p_to_s,y_s_to_p ! the FFTW plans for y transforms
    integer(kind=li) :: xy_p_to_s,xy_s_to_p ! the FFTW plans for xy transforms
    integer(kind=li) :: z_p_to_s,z_s_to_p ! the FFTW plans for z transforms
    integer :: x_modes,y_modes,z_modes,x_n,z_n  ! the global sizes of the arrays
    integer :: x_s_count,y_s_count,z_s_count ! the local sizes of the arrays
    integer :: x_p_count,y_p_count,z_p_count ! the local sizes of the arrays
  end type

  ! TODO: make this a local variable to be read into the subroutines
  type(transform_type) :: transform ! A global transform object

contains
  
  !/**********************************************************************|
  !| SUBROUTINES                                                          |
  !\**********************************************************************/
#include "transform_serial.F90"
  
end module