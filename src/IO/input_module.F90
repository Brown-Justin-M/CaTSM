module input_module
  !/**********************************************************************|
  !| INPUT MODULE                                                         |
  !| ------------                                                         |
  !| This module contains the interface to read an input file into a      |
  !| state_type derived type. The input module supports inputs in both    |
  !| physical and spectral space.                                         |
  !\**********************************************************************/
  use grid_module, only: dp,li,grid_type
  use params_module, only: params_type
  use fields_module, only: state_type
  use output_module, only: nc_ids_type,check,controllers_type
  implicit none
  private
  public read_input,open_inputs

contains
  subroutine open_inputs(grid,params,controls)
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

    integer :: ierror

    namelist /io/ controls

    read(params%unit,nml=io,iostat=ierror)  

    if (ierror /= 0) then 
      print*,"WARNING: While reading IO namelist, encountered error (",ierror,")"
      print*,"WARNING: Using defaults for IO."
    end if
    close(params%unit)
  end subroutine

#include "input_netcdf.F90"

end module