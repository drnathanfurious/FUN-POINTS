module stick_module

  use points_module, only : distance

  implicit none

  type stick_type
    real,pointer :: point1(:), point2(:)  ! the 2 points that make up the stick
    real :: length  ! the length of the stick
    integer :: group  ! the group identifier that the stick is in
  end type stick_type

contains


  function MakeStick (point1, point2) result (stick)
    type(stick_type) :: stick
    real,target :: point1(:), point2(:)

    stick%point1 => point1
    stick%point2 => point2
    stick%length = distance(point1,point2)
    stick%group = 0 ! null this out for now...
  end function MakeStick


  ! returns a set of sticks created from the points given
  function BundleSticksFromPoints (points) result (bundle)
    real :: points(:,:)
    type(stick_type) :: bundle(size(points(:,1)) * (size(points(:,1)) - 1) / 2)
    integer :: N, i, j, k

    N = size(points(:,1))

    k = 1
    do i=1,N-1
      do j=i+1,N
        bundle(k) = MakeStick(points(i,:), points(j,:))
        k = k + 1
      end do
    end do
        
  end function BundleSticksFromPoints


  ! Naively assign sticks to be in different groupings by lengths
  subroutine InitializeStickGroups (sticks)
    type(stick_type) :: sticks(:)
    integer :: num_groups, i, stick, group

    num_groups = size(sticks) - 1

    ! this little trick goes through and groups the sticks into N-1 length-
    ! groups. This should be good enough for an initial guess.
    stick = 1
    do group=1,num_groups
      do i = 1,group
        sticks(stick)%group = group
        stick = stick + 1
      end do
    end do
  end subroutine InitializeStickGroups


  ! pretty printer for sticks
  subroutine PrintStick (stick)
    type(stick_type) :: stick
    write (*,100) stick%point1, stick%point2, stick%length, stick%group
    100 format ("(",2f6.3,")","(",2f6.3,")",  &
                "  length: ", f6.3, "  group: ", i3)
  end subroutine PrintStick



  subroutine PrintSticks (sticks)
    type(stick_type) :: sticks(:)
    integer :: N,i
    N = size(sticks)

    do i=1,N
      call PrintStick(sticks(i))
    end do
  end subroutine PrintSticks

end module stick_module
