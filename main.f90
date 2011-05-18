! FUN-POINTS

! 
! Solution to this Erdos thing.  You could win $50!

program main

  ! this just makes things pretty and keeps all the derived types in one place
  use points_module
  use optimization_module

  implicit none

  integer,parameter :: number_of_points = 7
  integer,parameter :: dimensions = 2

  !!!!!
  ! a flattened version of all the positions of the points
  ! fortran matrices are just arrays... but in >>> COLUMN-MAJOR ORDER <<<
  ! let's just stick with arrays...
  real, target :: points(number_of_points,dimensions)  

  integer, target :: groupings(number_of_points*(number_of_points-1)/2)

  ! The optimization function data to be passed around during the minimization
  type(opt_function_type) :: opt_function

  real :: prev_minf



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Main routine for optimization
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! this random seed initializer has to be called from MAIN PROGRAM !!
  call init_random_seed()

  ! randomly select positions for the points in the plane
  points = InitializePoints(number_of_points,dimensions)

  opt_function%f_data%number_of_points = number_of_points
  opt_function%f_data%dimensions = dimensions

  ! group the pairs of points into their respective length-groups
  groupings = InitPairGroupings(number_of_points)
  opt_function%f_data%groupings => groupings


  prev_minf = 1000.0

  !call PrintPoints(points)
  call OptimizeGroupings (prev_minf, points, opt_function)
  !call PrintPoints(points)
  call PrintResults(points)










contains

  ! pretty print the results as pairs of points (lines)
  subroutine PrintResults(points) 
    real, intent(in), dimension(:,:) :: points
    integer :: i,j,N
    !real :: sticks(size(points(:,1))*(size(points(:,1))-1)/2)
    !integer :: natural_order(size(points(:,1))*(size(points(:,1))-1)/2)
    N = size(points(:,1))
    ! calculate lengths of all the sticks
    do i=1,N-1
    do j=i+1,N

    write(*,200) &
    i, j,    &
    points(i,:), points(j,:),   &
    distance(points(i,:), points(j,:))

    200 format('(',2I3') :: ', '(',2f6.3,')', '(',2f6.3,')', " ---> ", f9.4)

    end do
    end do

  end subroutine PrintResults


end program main

