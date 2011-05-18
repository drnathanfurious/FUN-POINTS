
module optimization_module
  use points_module
  use sort_module

  implicit none

 include 'nlopt.f' ! use nlopt - non-linear optimization routines

 !external :: f_opt ! optimization function 

  ! derived type for data passed to the nlopt function (possible future use)
  type nlo_fdata_type
    integer :: number_of_points
    integer :: dimensions
    ! list of groups to which each pair of points belongs. this should be
    ! N*(N-1)/2 x 3, where each row contains the index to pointA, pointB, and
    ! the group in which it belongs
    ! 
    ! i.e. 
    ! [[1,2,1], [1,3,2], [1,4,2], [2,3,3]]
    ! means the connection from point 1 to 2 is in group 1,
    !   the connection from point 1 to 3 is in group 2, etc.
    !
    ! (1,2) = group 1
    ! (1,3) = group 2
    ! (1,4) = group 2 
    ! (2,3) = group 3
    integer,pointer :: groupings(:)
  end type nlo_fdata_type

  type opt_function_type
    integer :: opt ! optimization guy
    integer :: ires ! sucess or failure of nlopt 
    type(nlo_fdata_type) :: f_data ! some data you can pass to nlopt
    real :: minf ! minimum fitness found as a result of nlopt's routines
  end type opt_function_type

  ! subroutines made avaliable to alles
  contains


  ! initialize optimization function
  subroutine InitOptimizationFunction (opt_func)
    type(opt_function_type), intent(inout) :: opt_func

    ! create optimization function
    opt_func%opt = 0

    call nlo_create(opt_func%opt, NLOPT_LN_NELDERMEAD,  &
        opt_func%f_data%number_of_points * opt_func%f_data%dimensions)

    ! lower bounds
    call nlo_set_lower_bounds1(opt_func%ires, opt_func%opt, 0.0d0)
    if(opt_func%ires.lt.0) then
      write(*,*) "set_lower_bounds failed"
    end if

    ! upper bounds
    call nlo_set_upper_bounds1(opt_func%ires, opt_func%opt, 1.0d0)
    if(opt_func%ires.lt.0) then
      write(*,*) "set_upper_bounds failed"
    end if

    ! want to minimize the objective function (as opposed to maximize)
    call nlo_set_min_objective(opt_func%ires, opt_func%opt, f_opt, opt_func%f_data)
    if(opt_func%ires.lt.0) then
      write(*,*) "set_min_objective failed"
    end if

    !call nlo_set_maxeval(opt_func%ires, opt_func%opt, 10000000)

    ! stopping criteria... 
    ! not sure whether ftol or xtol works better
    call nlo_set_xtol_abs1(opt_func%ires, opt_func%opt, 1.0E-14)
    if(opt_func%ires.lt.0) then
      write(*,*) "set_xtol_abs1 failed"
    end if

    !call nlo_set_ftol_abs(opt_func%ires, opt_func%opt, 0.00000000001d0)
    !if(opt_func%ires.lt.0) then
    !  write(*,*) "set_xtol_abs1 failed"
    !end if

  end subroutine InitOptimizationFunction




  ! The main routine for running an optimization on a set of points to get a
  ! local minimum
  subroutine RunOptimization (opt_function, points)
    type(opt_function_type) :: opt_function
    real :: points(:,:)  

    opt_function%ires=0

    ! initialize the optimization routine
    call InitOptimizationFunction(opt_function)

    !!!!!!!!!!!!!!!!!!!!!! 
    ! start optimization !
    call nlo_optimize(opt_function%ires, opt_function%opt, points, opt_function%minf)

    if (opt_function%ires.lt.0) then 
      write(*,*) "nlopt failed! This reflects poorly on you."
    else
      !write(*,*) "found min at "
      !call PrintResults(points)
      !  call PrintSticks(sticks)
      !do i=1,number_of_points
      !  write(*,*) points(i,:)
      !end do
      write(*,*) "# minimum ", opt_function%minf
    end if

    !!!!! clean up
    call nlo_destroy(opt_function%opt)

  end subroutine RunOptimization




  ! fitness function for nlopt
  subroutine f_opt(fitness, number_of_points, points, grad, need_gradient, f_data)
    type(nlo_fdata_type),intent(in) :: f_data
    real, intent(out) :: fitness ! fitness
    integer,intent(in) :: number_of_points
    real,intent(in) :: points(f_data%number_of_points,f_data%dimensions)
    real,intent(in) :: grad(f_data%number_of_points)
    integer,intent(in) :: need_gradient

    ! distances is a list of all the distances
    real :: distances(f_data%number_of_points * (f_data%number_of_points-1) / 2)

    integer :: natural_order(f_data%number_of_points * (f_data%number_of_points-1) / 2)

    distances = CalculateDistances(points)

    ! Sort the distances.  The variable natural_order is returned by the sort
    ! command, it tells how the rows were changed.
    call Qsort(distances,natural_order)

    fitness = CalcFitness(f_data%number_of_points,distances,f_data%groupings)

  end subroutine f_opt






  function CalcFitness (N,distances,groupings) result (fitness)
      real :: fitness
      integer :: N
      real :: distances(N*(N-1)/2)
      integer :: groupings(N*(N-1)/2)

      real :: avg_distances(N-1)
      integer :: i,group
      integer :: tally(N-1)

      tally(:) = 0

      do i=1,size(distances)
        group = groupings(i)
        avg_distances(group) = avg_distances(group) + distances(i)
        tally(group) = tally(group) + 1
      end do
      avg_distances(:) = avg_distances(:) / tally(:)
  

      ! based upon this average, the fitness can be calculated
      fitness=0.0
      do i=1,size(distances)
        group = groupings(i)
        fitness = fitness + abs((distances(i) - avg_distances(group))/(0.5*(distances(i) + avg_distances(group))) )
      end do
  
      ! add a penalty for two groups being close to each other in length?
      do i=1,N-2
        fitness = fitness + abs(avg_distances(i)-avg_distances(i+1))
      end do
        fitness = fitness + abs(avg_distances(1)-avg_distances(N-1))

  end function CalcFitness

  
end module optimization_module
