module points_module

  implicit none

contains

  ! Initialize the random seed based on system clock.  Grabbed this off of
  ! some forum somewhere, but seems to be the standard way of doing things.
  subroutine init_random_seed()
    integer :: i, n, clock 
    integer, dimension(:), allocatable :: seed 
    call random_seed(size = n) 
    allocate(seed(n)) 
    call system_clock(count=clock) 
    seed = clock + 37 * (/ (i - 1, i = 1, n) /) 
    call random_seed(put = seed) 
    deallocate(seed) 
  end subroutine init_random_seed

  ! gives each point in the system a random position
  function InitializePoints (num_points,dimensions) result (points)
    integer :: i,j, num_points,dimensions
    real :: points(num_points,dimensions)
    real :: rnum  ! random number

    ! seed some random points
    do i=1,num_points
      ! give random a random coordinate in N-dimensions of the problem
      do j=1,dimensions
        call random_number(rnum)
        points(i,j) = rnum
      end do
    end do
  end function InitializePoints


  ! calculate the distance between two points in D dimensions
  function distance (point1, point2)
    real :: point1(:), point2(:)
    real :: distance
    real :: temp(size(point1)) ! the distance vector between the two points

     temp = point2(:) - point1(:)
     temp = temp * temp 
     distance = sqrt(sum(temp)) ! sqrt of the sum of the squares... the norm
  end function distance


  function CalculateDistances (points) result (sticks)
    real :: points(:,:)
    real :: sticks(size(points(:,1))*(size(points(:,1))-1))
    integer :: i,j,k, N

    N = size(points(:,1))

    ! first calculate the lengths of all the sticks
    k=1
    do i=1,N-1
      do j=i+1,N
        sticks(k) = distance(points(i,:), points(j,:))
        k=k+1
      end do
    end do

  end function CalculateDistances

  ! given an adjacency matrix with the upper triangle filled out,
  ! here we calculate the averages of distances for each *row*
  ! in the matrix. 
  function CalculateDistanceAverages (matrix) result (distances)
    real :: matrix(:,:) ! this is a square matrix
    integer :: row, col, N
    ! for N points, there's N-1 distances
    real :: distances(size(matrix(:,1))-1)    

    N = size(matrix(:,1))-1

    ! summation and averaging along each row
    do row=1,N
      distances(row) = sum(matrix(row,:))/(N-row+1)
    end do

  end function CalculateDistanceAverages


  ! fills out the adjaceny matrix's upper half (not the diagonal) with 
  ! inter-point distances
  function UpdateAdjacencyMatrix (points) result(matrix)
    real :: points(:,:)
    integer :: i,j, N
    real :: matrix (size(points(:,1)),size(points(:,1)))

    N = size(points(:,1))

    matrix(:,:) = 0.0d0 ! first null out the connectivities/distances

    do i=1,N-1
      do j=i+1,N ! this avoids the matrix diagonal, 
                 !! and the lower half of the matrix.
        matrix(i,j) = distance (points(i,:), points(j,:))
      end do
    end do
  end function UpdateAdjacencyMatrix


  ! pretty printer for a single point
  subroutine PrintPoint (point)
    real :: point(:)
    integer :: i

    write (*,"(A)",advance="no") "("
    do i=1,size(point)
      write (*,"(f6.3)",advance="no") point(i)
    end do
    write (*,*) ")"
  end subroutine PrintPoint


  ! print all the points of a group
  subroutine PrintPoints (points)
    real :: points (:,:)
    integer :: i

    do i=1,size(points(:,1))
      call PrintPoint (points(i,:))
    end do
  end subroutine PrintPoints


  ! pretty printer for inter-point (scalar) distances
  function PrintDistances (matrix) result (stat)
    real :: matrix(:,:)
    integer :: i,j, stat, N

    N = size(matrix(:,1))

    write (*,200), N
    200 format("Number of points in the system: ",I2)

    do i=1,N-1
      do j=i+1, N
        write (*,100) i,j,matrix(i,j)
      end do
    end do

    100 format("  (",I2,"--",I2,")  ",f9.5)

    stat = 1
  end function PrintDistances



end module points_module
