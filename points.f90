module points_module

  use sort_module
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
        points(i,j) = ceiling(rnum*10)
      end do
    end do
  end function InitializePoints

  ! initialize the groupings of pair of points
  function InitPairGroupings (N) result (groupings)
    integer :: groupings(N*(N-1)/2,3)
    integer :: N, i,j,pair,tmp
    real :: rnum

    pair = 1
    do i=1,N-1
      do j = i+1,N
        groupings(pair,1) = i
        groupings(pair,2) = j
        groupings(pair,3) = i
        pair = pair + 1
      end do
    end do

    ! Shuffle the groupings
    do i=1, size(groupings(:,1))
      call random_number(rnum)
      rnum = rnum * size(groupings(:,1))
      j = ceiling(rnum)
      tmp = groupings(j,3)
      groupings(j,3) = groupings(i,3)
      groupings(i,3) = tmp
    end do

  end function InitPairGroupings


  ! grab the pairs of points for a particular group
  ! the returned list from this function is a list of pairs of points
  ! ... rather, it's a (N-group) x 2 list, where each row contains the index
  ! to point1 and point2 for a pair in a given group
  function GrabGroupPairs (group,N,groupings) result (pairs)
    integer :: group ! the group that we need to isolate
    integer :: N     ! total number of points
    integer :: pairs (N-group,2)
    integer :: groupings (N*(N-1)/2, 3)
    integer :: i,pair

    write (*,*) 

    pair = 1
    do i=1,size(groupings(:,1))
      if (groupings(i,3) == group) then
        pairs(pair,:) = groupings(i,1:2)
        pair = pair + 1
      end if
    end do

  end function 


  ! calculate the distance between two points in D dimensions
  function distance (point1, point2)
    real :: point1(:), point2(:)
    real :: distance
    real :: temp(size(point1)) ! the distance vector between the two points

     temp = point2(:) - point1(:)
     temp = temp * temp 
     distance = sqrt(sum(temp)) ! sqrt of the sum of the squares... the norm
  end function distance

  ! checks for more than two colinear points in 2D
  ! returns 1 if no more than two colinear points
  ! returns 0 if three or more colinear points
  ! this function only works in 2D!
  function CheckCollinear(points) result (is)
    real :: points(:,:)
    integer :: N,i,j,k,is
    N = size(points(:,1))
    is = 1
    ! same formula to calculate area of triangle in 2D.  If two points are
    ! collinear the area of the triangle should be zero. Here i,j,and k
    ! represent groupings of three.
    do i=1,N-1
      do j=i+1,N 
        do k=j+1,N
          if(points(i,1)*(points(j,2)-points(k,2))+&
             points(j,1)*(points(k,2)-points(i,2))+&
             points(k,1)*(points(i,2)-points(j,2)).eq.0) then
            is = 0
            return
          end if
        end do
      end do
    end do
  end function CheckCollinear

  ! Check whether there are more than three points on any circle.  This
  ! function only works in 2D!  Return 0 if the condition is NOT
  ! satisfied, and 1 if the condition is satisfied
  function CheckCircle(points) result (is)
    real :: points(:,:)
    integer :: i,j,k,l,N,is
    real :: d,c(2),r
    N = size(points(:,1))
    is = 1

    ! unique groupings of three points
    do i=1,N-1
      do j=i+1,N 
        do k=j+1,N
          ! area of the triangle whose points are the circumcircle
          d=2*(&
            points(i,1)*(points(j,2)-points(k,2))+&
            points(j,1)*(points(k,2)-points(i,2))+&
            points(k,1)*(points(i,2)-points(j,2)))
          ! center of circle defined by three points - x
          c(1)=(&
            (points(i,2)**2+points(i,1)**2)*(points(j,2)-points(k,2))+&
            (points(j,2)**2+points(j,1)**2)*(points(k,2)-points(i,2))+&
            (points(k,2)**2+points(k,1)**2)*(points(i,2)-points(j,2))&
          )/d
          ! center of circle defined by three points - y
          c(2)=(&
            (points(i,2)**2+points(i,1)**2)*(points(k,2)-points(j,2))+&
            (points(j,2)**2+points(j,1)**2)*(points(i,2)-points(k,2))+&
            (points(k,2)**2+points(k,1)**2)*(points(j,2)-points(i,2))&
          )/d
          ! radius of the circle
          r=distance(c,points(i,:))
          ! if the radius of the extra points is equal to the radius of the
          ! circle, then they lie on that circle and we have to regect this
          ! solution
          do l=k+1,N
            if(abs(distance(c,points(l,:))-r).lt.epsilon(c(1))) then
              is = 0
              return
            end if

          end do
        end do
      end do
    end do

  end function CheckCircle


  ! check the grouping condition.  For N points we want N-1 seperate
  ! groups.  Return 0 if the condition is NOT satisfied, return 1 if it is.
  function CheckGroupings(points) result (is)
    real :: points(:,:)
    real :: distances(size(points(:,1))*(size(points(:,1))-1)/2)
    integer :: check(size(points(:,1))*(size(points(:,1))-1)/2)
    integer :: i,j,k,is,N,Ndist
    real :: tmp1
    N = size(points(:,1))

    distances =  CalculateDistances(points)
    check(:) = 0
    is = 1
    ! count the number of grouped pairs and put them in the check array
    do i=1,size(distances)
      j = count(distances.eq.distances(i))
      check(j) = check(j) + 1
    end do

    ! check the check...
    do i=1,N-1
      if(check(i).ne.i) then
        is = 0
        return
      end if
    end do

  end function CheckGroupings


  function CalculateDistances (points) result (distances)
    real :: points(:,:)
    real :: distances(size(points(:,1))*(size(points(:,1))-1)/2)
    integer :: i,j,k, N

    N = size(points(:,1))
    distances(:) = 0.0d0

    ! first calculate the lengths of all the sticks
    k=1
    do i=1,N-1
      do j=i+1,N
        distances(k) = distance(points(i,:), points(j,:))
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

    !write (*,"(A)",advance="no") "("
    do i=1,size(point)
      write (*,"(f6.3)",advance="no") point(i)
    end do
    write (*,*)
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
