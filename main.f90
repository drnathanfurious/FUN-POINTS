! FUN-POINTS
! 
! Solution to this Erdos thing.  You could win $50!

program main

  ! this just makes things pretty and keeps all the derived types in one place
  use points_module
  use optimization_module
  use sort_module
  implicit none


  integer,parameter :: number_of_points = 5
  integer,parameter :: dimensions = 2


  !!!!!
  ! a flattened version of all the positions of the points
  ! fortran matrices are just arrays... but in >>> COLUMN-MAJOR ORDER <<<
  ! let's just stick with arrays...
  real, target :: points(number_of_points,dimensions)  
  ! using the idea of an upper triangular adjaceny matrix to track
  ! distances between points...
  real,target :: adjacency_matrix (number_of_points,number_of_points)

  type(opt_function_type) :: opt_function
  
  integer :: i,j,k
  real :: tmp

  opt_function%ires=0

  ! randomly select points in the plane
  call init_random_seed()
  points = InitializePoints(number_of_points,dimensions)

  ! assign the points and adjacency matrix to the opt_function structure
  opt_function%f_data%points => points
  opt_function%f_data%adjacency_matrix => adjacency_matrix

  ! initialize the optimization routine
  call InitOptimizationFunction(opt_function)

  !!!!! start optimization
  call nlo_optimize(opt_function%ires, opt_function%opt, points, opt_function%minf)

  if (opt_function%ires.lt.0) then 
    write(*,*) "nlopt failed! This reflects poorly on you."
  else
    !write(*,*) "found min at "
    call PrintResults(points)
    !do i=1,number_of_points
    !  write(*,*) points(i,:)
    !end do
    write(*,*) "# minimum ", opt_function%minf
  end if

  !!!!! clean up
  call nlo_destroy(opt_function%opt)

  contains

  ! pretty print the results as pairs of points (lines)
  subroutine PrintResults(points) 
    real, intent(in), dimension(:,:) :: points
    integer :: i,j,k,N
    real :: sticks(size(points(:,1))*(size(points(:,1))-1)/2)
    integer :: natural_order(size(points(:,1))*(size(points(:,1))-1)/2)
    N = size(points(:,1))
    ! calculate lengths of all the sticks
    k=1
    do i=1,N-1
      do j=i+1,N
        sticks(k) = distance(points(i,:), points(j,:))
        write(*,*) points(i,:), sticks(k)
        write(*,*) points(j,:), sticks(k)
        k=k+1
      end do
    end do

  end subroutine PrintResults


end program main

