
module optimization_module
  use points_module

  implicit none

 include 'nlopt.f' ! use nlopt - non-linear optimization routines

 !external :: f_opt ! optimization function 

  ! derived type for data passed to the nlopt function (possible future use)
  type nlo_fdata_type
    real,pointer :: points(:,:) ! positions in the system
    real,pointer :: adjacency_matrix(:,:) ! adj. matrix of inter-point distances
  end type nlo_fdata_type

  type opt_function_type
    integer :: opt ! optimization guy
    integer :: ires ! sucess or failure of nlopt 
    type(nlo_fdata_type) :: f_data ! some data you can pass to nlopt
    real :: minf ! minimum fitness found as a result of nlopt's routines
  end type opt_function_type

  ! subroutines made avaliable to alles
  contains

  ! dumby function
  ! Print the results to a file such that you can generate a pretty colored 
  ! plot.  Gnuplot will work for this.
  !subroutine PrintResults(points)
  ! 
  !end subroutine PrintResults

  ! initialize optimization function
  subroutine InitOptimizationFunction (opt_func)
    type(opt_function_type), intent(inout) :: opt_func
    integer :: N, D
    N = size(opt_func%f_data%points(:,1))
    D = size(opt_func%f_data%points(1,:))

    ! create optimization function
    opt_func%opt = 0

    call nlo_create(opt_func%opt, NLOPT_LN_NELDERMEAD, D*N)

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
    call nlo_set_xtol_abs1(opt_func%ires, opt_func%opt, 1.0E-12)
    if(opt_func%ires.lt.0) then
      write(*,*) "set_xtol_abs1 failed"
    end if

    !call nlo_set_ftol_abs(opt_func%ires, opt_func%opt, 0.00000000001d0)
    !if(opt_func%ires.lt.0) then
    !  write(*,*) "set_xtol_abs1 failed"
    !end if

  end subroutine InitOptimizationFunction

  ! fitness function
  function CalcFitness(points, matrix) result (error)
    real :: points(:,:)
    integer :: i,j,N
    real :: matrix(:,:)
    real :: error
    real :: distances(size(points(:,1))-1)  !!! see how N is set... below

    N = size(points(:,1)) !!! check this... is the extent to N or N-1?

    distances = CalculateDistanceAverages(matrix)

    error = 0
    ! first contribution: how close the line distance is to other members
    ! in its group
    do i=1, N-1
      do j=i+1, N
        error = error + abs(matrix(i,j)-distances(i))**2
      end do
    end do

    ! the above tends to converge on results where all the points are
    ! exactly the same, so introduce a reward for not being too close to
    ! other distances
    do i=1, N-2
       error = error - abs(distances(i)-distances(i+1))
    end do
    error = error - abs(distances(N-2)-distances(1))

    ! penalty for very small distances
    do i=1, N-1
       if(distances(i) .lt. 0.001d0) error = error+10000
    end do

  end function CalcFitness 



  ! fitness function for nlopt
  subroutine f_opt(res, en, x, grad, need_gradient, f_data)
    type(nlo_fdata_type) :: f_data
    integer :: en,need_gradient,i,j,k
    real :: x(en), grad(en)
    real, intent(out) :: res
  
    ! convert this linear array of points back into a points structure
    ! (seems wasteful)
    ! it's not really very wasteful... but it is unecessary ~es
    f_data%adjacency_matrix = UpdateAdjacencyMatrix(f_data%points)
    ! All the important fitness stuff is calculated in CalcFitness
    res = CalcFitness(f_data%points, f_data%adjacency_matrix)
  end subroutine f_opt
  
end module optimization_module
