! FUN-POINTS
! 
! Solution to this Erdos thing.  You could win $50!

program main

  ! this just makes things pretty and keeps all the derived types in one place
  use data_types_module

 implicit none

 include 'nlopt.f' ! use nlopt - non-linear optimization routines

 type(point_type) :: points(N) ! set of points in the system

 ! using the idea of an upper triangular adjaceny matrix to track
 ! distances between points...
 real, dimension(N,N) :: adjacency_matrix

 ! a flattened version of all the positions of the points
 real :: x(N*D)
 integer :: i

 type(opt_function_type) :: opt_function

 external :: f_opt ! optimization function 

 opt_function%ires=0

 ! randomly select points in the plane
 call InitializePoints(points)

 ! generate the connectivity matrix filled with distances between points
 !call UpdateAdjacencyMatrix(points, adjacency_matrix)
   
 ! pretty print out the inter-point distances
 !i = PrintDistances(adjacency_matrix)

 ! here are the average distances for each point
 !write(*,*) CalculateDistanceAverages(adjacency_matrix)

 ! convert these initial random points to a linear array
 call PointsToLinearArray(points,x)

 ! initialize the optimization routine
 call InitOptimizationFunction (opt_function)

 ! assign the points and adjacency matrix to the opt_function structure
 opt_function%f_data%points = points
 opt_function%f_data%adjacency_matrix = adjacency_matrix

 !!!!!!!!!!!!!!!!!!!! some diagnostics !!!!!!!!!!!!!!!!!!
 write (*,*) "ires: ", opt_function%ires
 write (*,*) "opt: ", opt_function%opt       !! this should not be zero by now... if it's
                                    !! if it is, there's a problem.
 write (*,*) "flat position array (x): ", x
 !write (*,*) "minf: ", opt_function%minf
 write (*,*) "the points structure: ", points
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! start optimization
 call nlo_optimize(opt_function%ires, opt_function%opt, x, opt_function%minf)

 if (opt_function%ires.lt.0) then 
   write(*,*) "nlopt failed! This reflects poorly on you."
 else
   write(*,*) "found min at ", x
   write(*,*) "minimum ", opt_function%minf
 end if

 ! clean up
 call nlo_destroy(opt_function%opt)
 call LinearArrayToPoints(points,x)
 call UpdateAdjacencyMatrix(points, adjacency_matrix)
 !call PrintResults(points)
 i = PrintDistances(adjacency_matrix)

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

  ! pretty printer for inter-point (scalar) distances
  function PrintDistances (matrix) result (stat)
    real :: matrix(:,:)
    integer :: i,j, stat

    write (*,200), N
    write (*,210), D
    200 format("Number of points in the system: ",I2)
    210 format("Dimensionality of the sytem: ",I2)

    do i=1,N-1
      do j=i+1, N
        write (*,100) i,j,matrix(i,j)
      end do
    end do

    100 format("  (",I2,"--",I2,")  ",f9.5)

    stat = 1
  end function PrintDistances

end program main

! fitness function for nlopt
subroutine f_opt(res, en, x, grad, need_gradient, f_data)
  use data_types_module
  type(nlo_fdata_type) :: f_data
  integer :: en,need_gradient,i,j,k
  real :: x(en), grad(en)
  real, intent(out) :: res

  ! convert this linear array of points back into a points structure
  ! (seems wasteful)
  call LinearArrayToPoints(f_data%points,x)
  call UpdateAdjacencyMatrix(f_data%points, f_data%adjacency_matrix)
  ! All the important fitness stuff is calculated in CalcFitness
  res = CalcFitness(f_data%points, f_data%adjacency_matrix)
end subroutine f_opt
