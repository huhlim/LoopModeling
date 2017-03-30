!-------------------------------------------------------------------------------
MODULE LOGGER
!-------------------------------------------------------------------------------
use globals, only: log_level, king

implicit none
private

public :: write_log
public :: write_divider
public :: terminate_with_error

character(len=80) :: divider = &
    '--------------------------------------------------------------------------------'

!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
subroutine write_log(msg, level, me)
!-------------------------------------------------------------------------------
character(len=*), intent(in) :: msg
integer, intent(in), optional :: level
integer, intent(in), optional :: me

if (present(level)) then
    if (level > log_level) return
else
    if (log_level < 9) return
end if

if (present(me)) then
    if (me /= king) return
end if

write(*,"(A)") msg

end subroutine write_log
!-------------------------------------------------------------------------------
subroutine write_divider(level, me)
!-------------------------------------------------------------------------------
integer, intent(in), optional :: level
integer, intent(in), optional :: me

if (present(level)) then
    if (level > log_level) return
else
    if (log_level < 9) return
end if

if (present(me)) then
    if (me /= king) return
end if

write(*,"(A)") divider

end subroutine write_divider
!-------------------------------------------------------------------------------
subroutine terminate_with_error(msg)
!-------------------------------------------------------------------------------
character(len=*), intent(in) :: msg

write(*,'(A)') msg
stop

end subroutine terminate_with_error
!-------------------------------------------------------------------------------
END MODULE LOGGER
!-------------------------------------------------------------------------------
