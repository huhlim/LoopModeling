!-------------------------------------------------------------------------------
MODULE LOOP_MODELING
!-------------------------------------------------------------------------------
use globals
use random, only: rand
use mathfunction, only: combinations
use geometry, only: bond_length, bond_angle, torsion_angle, rotate_torsion, &
                    internal2cartesian, internal2cartesian_reverse
use loop_closure, only: max_soln, solve_tripep_closure

implicit none
private

public :: close_loop

CONTAINS
!-------------------------------------------------------------------------------
subroutine close_loop(protein, res_i, res_j)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein
integer, intent(in) :: res_i, res_j
integer :: n_closing, i_close
integer, allocatable :: closing_s(:,:)
real(dp) :: rv(3,0:4,3), b_len(6), b_ang(7), t_ang(2), r_anchor(3,4)
integer :: i_soln, n_soln, i_res, res_no, atm_no
real(dp), allocatable :: soln(:,:,:,:), s_ang(:,:,:)

call get_closing_residue_list(protein, res_i, res_j, n_closing, closing_s)
i_close = int(rand()*n_closing)+1
write(*,'(A6,2x,A,4(2x,I4))') 'REMARK', 'CLOSING', i_close, closing_s(:,i_close)

call get_tripep_geometry(protein, closing_s(:,i_close), rv, b_len, b_ang, t_ang)
r_anchor(:,1:2) = rv(:,1:2,1)
r_anchor(:,3:4) = rv(:,2:3,3)

call solve_tripep_closure(b_len, b_ang, t_ang, r_anchor, n_soln, soln)
if (n_soln == 0) then
    return
end if

allocate(s_ang(2,3,n_soln))
call get_torsion_angles_from_the_solutions(n_soln, soln, rv, s_ang)

do i_soln = 1, n_soln
    do i_res = 1, 3
        ! rotate phi
        res_no = closing_s(i_res,i_close)
        atm_no = 3
        call rotate_torsion(protein, res_no, atm_no, s_ang(1,i_res,i_soln))
        !
        ! rotate psi
        res_no = res_no + 1
        atm_no = 1
        call rotate_torsion(protein, res_no, atm_no, s_ang(2,i_res,i_soln))
    end do
end do

call internal2cartesian(protein, res_i, res_j)

if (allocated(soln)) deallocate(soln)
deallocate(s_ang)

end subroutine close_loop
!-------------------------------------------------------------------------------
subroutine get_closing_residue_list(protein, res_i, res_j, n_closing, closing_s)
!-------------------------------------------------------------------------------
type(protein_type), intent(in) :: protein
integer, intent(in) :: res_i, res_j
integer, intent(out) :: n_closing
integer, intent(out), allocatable :: closing_s(:,:)
integer :: i_res, i, n_cand
integer, allocatable :: cand(:), comb(:,:)

n_cand = 0
do i_res = res_i, res_j
    if (protein%residue(i_res)%res_name /= 'PRO') then
        n_cand = n_cand + 1
    end if
end do

i = 0
allocate(cand(n_cand))
do i_res = res_i, res_j
    if (protein%residue(i_res)%res_name /= 'PRO') then
        i = i + 1
        cand(i) = i_res
    end if
end do

call combinations(n_cand, 3, n_closing, comb)

allocate(closing_s(3, n_closing))
do i = 1, n_closing
    closing_s(:,i) = cand(comb(:,i))
end do

deallocate(cand)

end subroutine get_closing_residue_list
!-------------------------------------------------------------------------------
subroutine get_tripep_geometry(protein, closing, rv, b_len, b_ang, t_ang)
!-------------------------------------------------------------------------------
type(protein_type), intent(in) :: protein
integer, intent(in) :: closing(3)
real(dp), intent(out) :: rv(3,0:4,3)
real(dp), intent(out) :: b_len(6), b_ang(7), t_ang(2)
real(dp) :: rs(3,9)
integer :: i_res, ic, i

do i_res = 1, 3
    ic = closing(i_res)
    rv(:,1:3,i_res) = protein%residue(ic)%R(:,1:3)
    rv(:,0,i_res) = protein%residue(ic-1)%R(:,3)
    rv(:,4,i_res) = protein%residue(ic+1)%R(:,1)
end do
rs = reshape(rv(:,1:3,:), (/3,9/))

! virtual bond lengths
do i = 1, 6
    b_len(i) = bond_length(rs(:,i+1:i+2))
end do

! virtual bond angles
do i = 1, 7
    b_ang(i) = bond_angle(rs(:,i:i+2))
end do

! virtual torsion angles
t_ang(1) = torsion_angle(rs(:,2:5))
t_ang(2) = torsion_angle(rs(:,5:8))

end subroutine get_tripep_geometry
!-------------------------------------------------------------------------------
subroutine get_torsion_angles_from_the_solutions(n_soln, soln, rv, t_ang)
!-------------------------------------------------------------------------------
integer, intent(in) :: n_soln
real(dp), intent(in) :: soln(:,:,:,:), rv(3,0:4,3)
real(dp), intent(out) :: t_ang(:,:,:)
integer :: i_soln, i, j
real(dp) :: t_ang0(2,2,3), r_tor(3,4)

! compute the initial angles
do i = 1, 3
    do j = 1, 2
        r_tor = rv(:,j-1:j+2, i)
        t_ang0(1,j,i) = torsion_angle(r_tor)
        !
        if (i /= 1 .and. j == 1) then
            r_tor(:,1) = rv(:,3,i-1)
            r_tor(:,2:4) = rv(:,1:3,i)
        else if (i /= 3 .and. j == 2) then
            r_tor(:,1:3) = rv(:,1:3,i)
            r_tor(:,4) = rv(:,1,i+1)
        else
            cycle
        end if
        t_ang0(2,j,i) = torsion_angle(r_tor)
    end do
end do

do i_soln = 1, n_soln
    do i = 1, 3
        do j = 1, 2
            if (i == 1 .and. j == 1) then
                r_tor(:,1) = rv(:,0,i)
                r_tor(:,2:4) = soln(:,1:3,i,i_soln)
                t_ang(j,i,i_soln) = torsion_angle(r_tor)
                cycle
            else if (i == 3 .and. j == 2) then
                r_tor(:,1:3) = soln(:,1:3,i,i_soln)
                r_tor(:,4) = rv(:,4,i)
                t_ang(j,i,i_soln) = torsion_angle(r_tor)
                cycle
            end if
            !
            if (j == 1) then
                r_tor(:,1) = soln(:,3,i-1,i_soln)
                r_tor(:,2:4) = soln(:,1:3,i,i_soln)
            else
                r_tor(:,1:3) = soln(:,1:3,i,i_soln)
                r_tor(:,4) = soln(:,1,i+1,i_soln)
            end if
            t_ang(j,i,i_soln) = t_ang0(1,j,i) + (torsion_angle(r_tor)-t_ang0(2,j,i))
        end do
    end do
end do

end subroutine get_torsion_angles_from_the_solutions
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END MODULE LOOP_MODELING
!-------------------------------------------------------------------------------
