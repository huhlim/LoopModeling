!-------------------------------------------------------------------------------
MODULE RANDOM
!-------------------------------------------------------------------------------
use globals, only: me, dp

implicit none
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_random()
!-------------------------------------------------------------------------------
integer :: i, n, clock
integer, allocatable :: seed(:)

call random_seed(size=n)
allocate(seed(n))

call system_clock(count=clock)

seed = clock + 37 * (/ (i-1, i=1, n) /) + me
call random_seed(put=seed)

deallocate(seed)

end subroutine initialize_random
!-------------------------------------------------------------------------------
function rand()
!-------------------------------------------------------------------------------
real(dp) :: rand

call random_number(rand)

end function rand
!-------------------------------------------------------------------------------
END MODULE RANDOM
!-------------------------------------------------------------------------------
