! FUN-POINTS
! 
! Solution to this Erdos thing.  You could win $50!
program main
 implicit none

 ! derived type for representing points in the 2-d plane
  type point_type
    integer :: id ! some identifying value
    real*8, dimension(2) :: position ! 2d (x,y) position
  end type point_type

  integer, parameter :: N = 4 ! number of points
  integer, parameter :: D = 2 ! dimensionality of the problem
  integer :: i,j,k,l,m ! variables to iterate over

  type(point_type), dimension(N) :: points ! all the points

  ! using the idea of an upper triangular adjaceny matrix to track
  ! distances between points...
  real*8, dimension(:,:), allocatable :: adjacency_matrix

  ! randomly select points in the plane
  call InitializePoints(points)

  ! generate the connectivity matrix filled with distances between points
  call UpdateAdjacencyMatrix(points, adjacency_matrix)

  i = PrintDistances (adjacency_matrix)

 !!!!! >>> I can't understand what this this is doing... ~ess !!!!!!

 ! N-i is the number of partners in your group
 ! coordinates i,j are the pairs of points your segment is attached to
 ! l is the segment number currently working on
 ! m is the number of iterations
 !l=1
 !do i=1, N-1
  !do k=1, N-i
   !! go over all other segments in your group
   !do j=i+1, N 
    !if(j-i .eq. k .and. i .ne. N-1) cycle ! don't choose yourself
    !write(*,*) i,j,k,l    !!!!! >>> what are i, j, k, and l?
   !end do
   !l=l+1
  !end do
 !end do

 contains 
  


  ! gives each point in the system a random position
  subroutine InitializePoints (points)
    type(point_type), dimension(:), intent(inout) :: points
    real*8 :: rnum  ! random number

    ! seed some random points
    call random_seed
    do i=1, size(points)
      points(i)%id = i
      call random_number(rnum)
      points(i)%position(1)=rnum  ! x coordinate
      call random_number(rnum)
      points(i)%position(2)=rnum  ! x coordinate
    end do
  end subroutine InitializePoints



   ! fills out the adjaceny matrix's upper half (not the diagonal) with 
   ! inter-point distances
  subroutine UpdateAdjacencyMatrix (points, matrix)
    type(point_type), dimension(:), intent(in) :: points
    integer :: num_points
    integer :: i,j
    real*8, dimension(:,:), intent(out), allocatable :: matrix

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


  ! calculate the distance between two points in the 2D plane
  function distance (point1, point2)
    real*8 :: distance
    real*8, dimension(2) :: temp ! the distance vector between the two points
    type(point_type) :: point1, point2

     temp = point2%position(:) - point1%position(:)
     temp = temp * temp 
     distance = sqrt(sum(temp)) ! sqrt of the sum of the squares... the norm
  end function distance

  function PrintDistances (matrix) result (stat)
    real*8, dimension(:,:) :: matrix
    integer :: i,j, num_points, stat

    num_points = size(matrix(1,:))
    write (*,200), num_points
    200 format("Number of points in the system: ",I2)

    do i=1,num_points-1
      do j=i+1, num_points
        write (*,100) i,j,matrix(i,j)
      end do
    end do

    100 format("  (",I2,"--",I2,")  ",f5.2)

    stat = 1
  end function PrintDistances


    

end program main
