! FUN-POINTS
! 
! Solution to this Erdos thing.  You could win $50!

program main
 implicit none
 include 'nlopt.f' ! use nlopt
 integer, parameter :: D = 2 ! dimensionality of the problem
 integer, parameter :: N = 6 ! number of points

 ! derived type for representing points
  type point_type
    integer :: id ! some identifying value
    real, dimension(D) :: position ! position
  end type point_type

  ! derived type for data passed to the nlopt function (possible future use)
  type nlo_fdata_type
    type(point_type), dimension(N) :: points ! points
    real, dimension(:,:), allocatable :: adjacency_matrix ! adj. matrix
  end type nlo_fdata_type

  integer :: i,j,k,l,m ! variables to iterate over

  type(point_type), dimension(N) :: points ! set of points
  type(point_type), dimension(N) :: best_points ! set of best points

  ! using the idea of an upper triangular adjaceny matrix to track
  ! distances between points...
  real, dimension(:,:), allocatable :: adjacency_matrix
  real, dimension(:), allocatable :: x

  integer :: opt ! optimization guy
  integer :: ires ! sucess or failure of nlopt 
  real :: minf ! minimum fitness found as a result of nlopt's routines
  real :: f_data ! some data you can pass to nlopt
  external :: f_opt ! optimization function 
  ires=0


  call InitializePoints(points)
  ! randomly select points in the plane

  ! generate the connectivity matrix filled with distances between points
  call UpdateAdjacencyMatrix(points, adjacency_matrix)
    
  ! pretty print out the inter-point distances
  i = PrintDistances(adjacency_matrix)

  ! here are the average distances for each point
  write(*,*) CalculateDistanceAverages(adjacency_matrix)

  ! convert these points to a linear array
  call PointsToLinearArray(x,points)
  write(*,*) x

  ! creat optimization function
  opt = 0
  call nlo_create(opt,NLOPT_GN_DIRECT_L,D*N)
  ! lower bounds
  call nlo_set_lower_bounds1(ires, opt, -100.0d0)
  if(ires.lt.0) then
    write(*,*) "set_lower_bounds failed"
  end if
  ! upper bounds
  call nlo_set_upper_bounds1(ires, opt, 100.0d0)
  if(ires.lt.0) then
    write(*,*) "set_upper_bounds failed"
  end if
  ! want to minimize the objective function (as opposed to maximize)
  call nlo_set_min_objective(ires, opt, f_opt, f_data)
  if(ires.lt.0) then
    write(*,*) "set_min_objective failed"
  end if
  ! stopping criteria... stop when the function changes less than the below
  ! amount 
  call nlo_set_ftol_rel(ires, opt, 0.001d0)
  if(ires.lt.0) then
    write(*,*) "set_xtol_abs1 failed"
  end if

  ! start optimization
  call nlo_optimize(ires, opt, x, minf)
  if (ires.lt.0) then 
    write(*,*) "nlopt failed! This reflects poorly on you."
  else
    write(*,*) "found min at ", x
    write(*,*) "minimum ", minf
  end if

  ! clean up
  call nlo_destroy(opt)


 contains 

  ! converts point datatypes to a 1d array
  subroutine PointsToLinearArray(x,points)
    type(point_type), dimension(:), intent(inout) :: points
    integer :: i,j,num_points
    real, dimension(:), intent(out), allocatable :: x

    num_points = size(points)
    allocate(x(N*D))

    do i=1, num_points
      do j=1, D
        x((i-1)*D+j)=points(i)%position(j)
        write(*,*) (i-1)*D+j,points(i)%position(j)
      end do
    end do
  end subroutine PointsToLinearArray
 

  ! gives each point in the system a random position
  subroutine InitializePoints (points)
    type(point_type), dimension(:), intent(inout) :: points
    real :: rnum  ! random number
    integer :: i,j

    ! seed some random points
    call random_seed
    do i=1, size(points)
      points(i)%id = i
      ! give random a random coordinate in D-dimensions
      do j=1, D
        call random_number(rnum)
        points(i)%position(j)=rnum  ! coordinate
      end do
    end do
  end subroutine InitializePoints



   ! fills out the adjaceny matrix's upper half (not the diagonal) with 
   ! inter-point distances
  subroutine UpdateAdjacencyMatrix (points, matrix)
    type(point_type), dimension(:), intent(in) :: points
    integer :: num_points
    integer :: i,j
    real, dimension(:,:), intent(out), allocatable :: matrix

    num_points = size(points)
    allocate (matrix(num_points,num_points))
    matrix(:,:) = 0.0d0 ! first null out the connectivities/distances

    do i=1,num_points-1
      do j=i+1,num_points ! this avoids the matrix diagonal, 
                          !! and the lower half of the matrix.
        matrix(i,j) = distance (points(i), points(j))
      end do
    end do
  end subroutine UpdateAdjacencyMatrix



  ! given an adjacency matrix with the upper triangle filled out,
  ! here we calculate the averages of distances for each *row*
  ! in the matrix. 
  function CalculateDistanceAverages (matrix) result (distances)
    real, dimension(:,:) :: matrix
    integer :: num_points, row, col
    real, dimension(size(matrix(1,:))-1) :: distances
    num_points = size(distances)

    ! summation and averaging along each row
    do row=1,num_points
      distances(row) = sum(matrix(row,:))/(num_points-row+1)
    end do

  end function CalculateDistanceAverages
    

  ! calculate the distance between two points in D dimensions
  function distance (point1, point2)
    real :: distance
    real, dimension(D) :: temp ! the distance vector between the two points
    type(point_type) :: point1, point2

     temp = point2%position(:) - point1%position(:)
     temp = temp * temp 
     distance = sqrt(sum(temp)) ! sqrt of the sum of the squares... the norm
  end function distance



  ! pretty printer for inter-point (scalar) distances
  function PrintDistances (matrix) result (stat)
    real, dimension(:,:) :: matrix
    integer :: i,j, num_points, stat

    num_points = size(matrix(1,:))
    write (*,200), num_points
    write (*,210), D
    200 format("Number of points in the system: ",I2)
    210 format("Dimensionality of the sytem: ",I2)

    do i=1,num_points-1
      do j=i+1, num_points
        write (*,100) i,j,matrix(i,j)
      end do
    end do

    100 format("  (",I2,"--",I2,")  ",f9.5)

    stat = 1
  end function PrintDistances

  ! incomplete routine to calculate the fitness
  function CalcFitness(points, matrix) result (error)
    type(point_type), dimension(:), intent(in) :: points
    integer :: num_points
    integer :: i,j
    real, dimension(:,:), intent(in) :: matrix
    real :: error

    real, dimension(size(matrix(1,:))-1) :: distances
    distances = CalculateDistanceAverages(matrix)

    num_points = size(points)

    error = 0
    do i=1,num_points-1
      do j=i+1, num_points
        error = error + abs(matrix(i,j)-distances(i))
      end do
    end do
  end function CalcFitness 

  
end program main

! fitness function for nlopt
! dummy routine for now... just tries to make all values 0.5
subroutine f_opt(res, n, x, grad, need_gradient, f_data)
  real :: f_data
  integer :: n,need_gradient,i
  real :: x(n), grad(n)
  real, intent(out) :: res
  res = 0
  do i=1,n
    res = res + abs(x(i)-0.5d0)
  end do

end subroutine f_opt
