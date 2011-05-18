! FUN-POINTS
! 
! Solution to this Erdos thing.  You could win $50!

program main

  ! this just makes things pretty and keeps all the derived types in one place
  use points_module
  use optimization_module

  implicit none

  integer,parameter :: number_of_points = 4
  integer,parameter :: dimensions = 2

  !!!!!
  ! a flattened version of all the positions of the points
  ! fortran matrices are just arrays... but in >>> COLUMN-MAJOR ORDER <<<
  ! let's just stick with arrays...
  real, target :: points(number_of_points,dimensions)  
  integer, target :: groupings(number_of_points*(number_of_points-1)/2,3)

  ! The optimization function data to be passed around during the minimization
  type(opt_function_type) :: opt_function


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

  ! run twice for good measure
  call RunOptimization (opt_function, points)
  call RunOptimization (opt_function, points)
  call PrintResults(points, groupings)







  contains

  ! pretty print the results as pairs of points (lines)
  subroutine PrintResults(points, groupings) 
    real, intent(in), dimension(:,:) :: points
  integer, target :: groupings(size(points(:,1))*(size(points(:,1))-1)/2,3)
    integer :: i,j,k,l,N
    N = size(points(:,1))

    ! write indicies, groupings, lengths
    do i=1, size(groupings(:,1))
      j = groupings(i,1) ! index of first point
      k = groupings(i,2) ! index of second point
      l = groupings(i,3) ! group this point is in
      write(*,200) &
        j, k, &
        points(j,:), points(k,:), &
        distance(points(j,:), points(k,:)), &
        l
      200 format('(',2I3') :: ', '(',2f6.3,')', '(',2f6.3,')', " ---> ", f9.4 ,2I3 )
    end do

  end subroutine PrintResults


end program main

