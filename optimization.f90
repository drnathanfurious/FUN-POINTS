
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
    ! groupings are ordered by point pairs. The first element is for the pair
    ! (1,2), the 2nd element for (1,3)... (1,N)
    ! then start (2,3), (2,4) ..., (2,N)....
    ! finally the last elements are ..., (N-3,N), (N-2,N), (N-1,N)
    ! each element contains the integer id of the group that the pair belongs to
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


  subroutine OptimizePoints (points, opt_function)
    type(opt_function_type) :: opt_function
    real :: points(:,:), prev_points(size(points(:,1)),size(points(1,:)))
    real :: prev_minf
    real :: run_min
    logical :: accepted ! if the move was accepted or not
    real,parameter :: tol = 1.0E-14 ! tolerance
    integer,parameter :: maxtrys = 250   ! max number of attempts to make before bailing
    integer :: attempt


    opt_function%minf = 1000.0  ! set this high to start out

    accepted = .true.

    do  while (accepted) 
      prev_minf = opt_function%minf
      prev_points = points


      run_min = opt_function%minf + 1000.0  ! guarantee that we loop through at least once!

      attempt = 1
      do while (abs(run_min - opt_function%minf) > tol .and. attempt < maxtrys )
        call RunOptimization (opt_function, points)
        run_min = opt_function%minf
        attempt = attempt + 1
        !write (*,100) attempt, opt_function%minf
        !100 format(1I, 2e11.3)
      end do
      !call OptimizePoints (points, opt_function)

      ! If the optimized fitness/energy is not as good as what we started with
      ! then restore the previous state
      if (prev_minf - opt_function%minf < tol) then
        opt_function%minf = prev_minf
        points = prev_points
        !write (*,*) "reject: ", prev_minf - opt_function%minf
        accepted = .false.

        ! otherwise, always accept a lower energy/fitness
      else 
        !write (*,*) "accept: ", prev_minf - opt_function%minf
        accepted = .true.
      end if

    end do


  end subroutine OptimizePoints




  subroutine MonteCarlo (points, opt_function)
    real :: points(:,:)  
    type(opt_function_type) :: opt_function

    real :: prev_points(size(points(:,1)), size(points(1,:)))
    integer :: prev_groupings(size(points(:,1))*(size(points(:,1))-1)/2)
    real :: prev_minf
    real :: monte, rnum
    real, parameter :: temp = 50.0

    integer :: accepted


    ! set up to restore the previous state if the move is rejected
    prev_minf = 1000.0
    prev_groupings = opt_function%f_data%groupings
    prev_points = points

    accepted = 0
    do while (accepted < 250)
    call OptimizePoints (points, opt_function)

    ! if the previous state was a lower/better fitness, then restore the state,
    ! and try altering the groupings again
    if (prev_minf < opt_function%minf) then

      !monte = prev_minf - opt_function%minf
      !monte = exp(monte/temp)
      !rnum = RandomReal(1.0)
      !write (*,*) monte, rnum

      ! accept the move with probability exp(-E/T)
      !if (rnum > monte) then
        !accepted = accepted + 1
        !prev_points = points
        !prev_groupings = opt_function%f_data%groupings
        !prev_minf = opt_function%minf
        !write (*,*) opt_function%minf
        !write (*,*) "accepted the monte!"
      ! otherwise reject the move
      !else
        points = prev_points
        opt_function%f_data%groupings = prev_groupings
        opt_function%minf = prev_minf
        !write (*,*) "rejected the monte..."
      !end if

      ! otherwise, accept the move, update the state, and try a new set of groupings
    else
      prev_points = points
      prev_groupings = opt_function%f_data%groupings
      prev_minf = opt_function%minf
      accepted = accepted + 1
      write (*,*) opt_function%minf
    end if

    ! try a new state
    opt_function%f_data%groupings = RandomlySwapGrouping (opt_function%f_data%groupings)

    end do

  end subroutine MonteCarlo



end module optimization_module
