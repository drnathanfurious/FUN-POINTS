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
  real :: prev_points(number_of_points, dimensions)

  integer, target :: groupings(number_of_points*(number_of_points-1)/2), i
  integer :: prev_groupings(size(groupings))

  ! The optimization function data to be passed around during the minimization
  type(opt_function_type) :: opt_function

  real :: prev_minf
  integer :: accepted



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

  ! set up to restore the previous state if the move is rejected
  prev_minf = 1000.0
  prev_groupings = groupings
  prev_points = points

  accepted = 0
  do while (accepted < 10)
    call OptimizePoints (points, opt_function)

    ! if the previous state was a lower/better fitness, then restore the state,
    ! and try altering the groupings again
    if (prev_minf < opt_function%minf) then
      points = prev_points
      groupings = prev_groupings
      opt_function%minf = prev_minf
    ! otherwise, update the state, and try a new set of groupings
    else
      prev_points = points
      prev_groupings = groupings
      prev_minf = opt_function%minf
      accepted = accepted + 1
    end if

    ! try a new state
    opt_function%f_data%groupings = RandomlySwapGrouping (opt_function%f_data%groupings)
    write (*,*) opt_function%minf

  end do


  call PrintPoints(points)
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

    200 format('(',2I3') :: ', '(',2f6.3,')', '(',2f6.3,')', " ---> ", f25.15)

    end do
    end do

  end subroutine PrintResults


end program main

