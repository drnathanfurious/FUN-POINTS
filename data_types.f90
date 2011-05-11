
module data_types_module
implicit none

 integer, parameter :: D = 2 ! dimensionality of the problem
 integer, parameter :: N = 6 ! number of points

 ! derived type for representing points
  type point_type
    integer :: id ! some identifying value
    real, dimension(D) :: position        ! position of a point
  end type point_type

  type opt_function_type
    integer :: opt ! optimization guy
    integer :: ires ! sucess or failure of nlopt 
    real :: f_data ! some data you can pass to nlopt
    real :: minf ! minimum fitness found as a result of nlopt's routines
  end type opt_function_type

  ! derived type for data passed to the nlopt function (possible future use)
  type nlo_fdata_type
    type(point_type) :: points(N) ! points
    real, allocatable :: adjacency_matrix(:,:) ! adj. matrix
  end type nlo_fdata_type


end module data_types_module
