!-------------------------------------------------------------------------------
MODULE LOOP_MODELING
!-------------------------------------------------------------------------------
use globals
use ran, only: random
use mathfunction, only: combinations
use geometry, only: bond_length, bond_angle, torsion_angle
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
real(dp) :: rv(3,3,3), b_len(6), b_ang(7), t_ang(2)
integer :: n_soln

call get_closing_residue_list(protein, res_i, res_j, n_closing, closing_s)
i_close = int(random()*n_closing)+1

call get_tripep_geometry(protein, closing_s(:,i_close), rv, b_len, b_ang, t_ang)

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
real(dp), intent(out) :: rv(3,3,3)
real(dp), intent(out) :: b_len(6), b_ang(7), t_ang(2)
real(dp) :: rs(3,9)
integer :: i_res, ic, i

do i_res = 1, 3
    ic = closing(i_res)
    rv(:,:,i_res) = protein%residue(ic)%R(:,1:3)
end do
rs = reshape(rv, (/3,9/))

! virtual bond lengths
do i = 1, 6
    b_len(i) = bond_length(rs(:,i+1:i+2))
end do

! virtual bond angles
do i = 1, 7
    b_ang(i) = pi - bond_angle(rs(:,i:i+2))
end do

! virtual torsion angles
t_ang(1) = torsion_angle(rs(:,2:5))
t_ang(2) = torsion_angle(rs(:,5:8))

end subroutine get_tripep_geometry
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END MODULE LOOP_MODELING
!-------------------------------------------------------------------------------
