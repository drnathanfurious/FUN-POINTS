! FUN-POINTS
! 
! Solution to this Erdos thing.  You could win $50!
program main
 implicit none
 integer, parameter :: N = 4 ! number of points
 integer, parameter :: D = 2 ! dimensionality of the problem
 integer :: i,j,k,l,m ! variables to iterate over
 real, dimension(N-1,D) :: point ! all the points
 real :: rnum
 
 ! seed some random points
 call random_seed
 do i=1, N
  call random_number(rnum)
  point(i,1)=rnum  ! x coordinate
  call random_number(rnum)
  point(i,2)=rnum  ! y coordinate
 end do

 ! N-i is the number of partners in your group
 ! coordinates i,j are the pairs of points your segment is attached to
 ! l is the segment number currently working on
 ! m is the number of iterations
 l=1
 do i=1, N-1
  do k=1, N-i
   ! go over all other segments in your group
   do j=i+1, N 
    if(j-i .eq. k .and. i .ne. N-1) cycle ! don't choose yourself
    write(*,*) i,j,k,l
   end do
   l=l+1
  end do
 end do

end program main
