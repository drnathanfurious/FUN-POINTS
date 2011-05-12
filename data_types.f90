
module data_types_module
implicit none

 integer, parameter :: D = 2 ! dimensionality of the problem
 integer, parameter :: N = 6 ! number of points

 ! derived type for representing points
  type point_type
    integer :: id ! some identifying value
    real, dimension(D) :: position        ! position of a point
  end type point_type

  ! derived type for data passed to the nlopt function (possible future use)
  type nlo_fdata_type
    type(point_type) :: points(N) ! points
    real, dimension(N,N) :: adjacency_matrix ! adj. matrix
  end type nlo_fdata_type

  type opt_function_type
    integer :: opt ! optimization guy
    integer :: ires ! sucess or failure of nlopt 
    type(nlo_fdata_type) :: f_data ! some data you can pass to nlopt
    real :: minf ! minimum fitness found as a result of nlopt's routines
  end type opt_function_type

  ! subroutines made avaliable to alles
  contains

  ! converts point datatypes to a 1d array
  subroutine PointsToLinearArray(points,x)
    type(point_type) :: points(:)
    integer :: i,j
    real :: x(N*D)

    do i=1, N
      do j=1, D
        x((i-1)*D+j)=points(i)%position(j)
        !write(*,*) (i-1)*D+j,points(i)%position(j)
      end do
    end do
  end subroutine PointsToLinearArray
 
  ! converts linear array to point type 
  ! to improve performance, might want to scrap this in favor of using a
  ! linear array of points everywere
  subroutine LinearArrayToPoints(points,x)
    type(point_type) :: points(N)
    integer :: i,j
    real, intent(in) :: x(N*D)
    do i=1, N
      do j=1, D
        points(i)%position(j)=x((i-1)*D+j)
      end do
    end do
  end subroutine LinearArrayToPoints

  ! gives each point in the system a random position
  subroutine InitializePoints (points)
    type(point_type), intent(inout) :: points(:)
    real :: rnum  ! random number
    integer :: i,j

    ! seed some random points
    call random_seed
    do i=1, N
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
    type(point_type), intent(in) :: points(N)
    integer :: i,j
    real, intent(out) :: matrix(N,N)

    matrix(:,:) = 0.0d0 ! first null out the connectivities/distances

    do i=1,N-1
      do j=i+1,N ! this avoids the matrix diagonal, 
                          !! and the lower half of the matrix.
        matrix(i,j) = distance (points(i), points(j))
      end do
    end do
  end subroutine UpdateAdjacencyMatrix

  ! given an adjacency matrix with the upper triangle filled out,
  ! here we calculate the averages of distances for each *row*
  ! in the matrix. 
  function CalculateDistanceAverages (matrix) result (distances)
    real :: matrix(:,:)
    integer :: row, col
    real :: distances(N-1)

    ! summation and averaging along each row
    do row=1,N
      distances(row) = sum(matrix(row,:))/(N-row+1)
    end do

  end function CalculateDistanceAverages
    
  ! calculate the distance between two points in D dimensions
  function distance (point1, point2)
    real :: distance
    real :: temp(D) ! the distance vector between the two points
    type(point_type) :: point1, point2

     temp = point2%position(:) - point1%position(:)
     temp = temp * temp 
     distance = sqrt(sum(temp)) ! sqrt of the sum of the squares... the norm
  end function distance

  ! fitness function
  function CalcFitness(points, matrix) result (error)
    type(point_type), intent(in) :: points(N)
    integer :: i,j
    real, intent(in) :: matrix(N,N)
    real :: error

    real :: distances(N-1)
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


end module data_types_module
