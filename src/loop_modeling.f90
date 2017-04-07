!-------------------------------------------------------------------------------
MODULE LOOP_MODELING
!-------------------------------------------------------------------------------
use globals
use mathfunction, only: combinations
use geometry, only: bond_length, bond_angle, torsion_angle, rotate_torsion, &
                    cartesian2internal, internal2cartesian, &
                    internal2cartesian_reverse, bound_angle, &
                    report_protein_geometry
use loop_closure, only: max_soln, solve_tripep_closure
use in_out, only: write_pdb

implicit none
private

public :: close_loop
public :: close_loop_complete

CONTAINS
!-------------------------------------------------------------------------------
subroutine close_loop(protein0, res_i, res_j, unperturbed)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein0
integer, intent(in) :: res_i, res_j
type(protein_type), intent(in), optional :: unperturbed
!
type(protein_type) :: protein
integer :: n_closing, i_close
integer, allocatable :: closing_s(:,:)
real(dp) :: rv(3,0:4,3), b_len(6), b_ang(7), t_ang(2), r_anchor(3,4)
integer :: i_soln, n_soln, i_res, res_no, atm_no
real(dp), allocatable :: soln(:,:,:,:), s_ang(:,:,:)
integer :: closing_min(3)
real(dp) :: s_ang_min(2,3), d_min, d

if (present(unperturbed)) then
    call link_broken_bond(protein0, res_j+1, unperturbed=unperturbed)
else
    call link_broken_bond(protein0, res_j+1)
end if

d_min = 999999999.d0
call get_closing_residue_list(protein0, res_i, res_j, n_closing, closing_s)
do i_close = 1, n_closing
    protein = protein0
    call internal2cartesian_reverse(protein, closing_s(3,i_close), res_j+1)
    !
    call get_tripep_geometry(protein0, closing_s(:,i_close), rv, b_len, b_ang, t_ang)
    r_anchor(:,1:2) = protein%residue(closing_s(1,i_close))%R(:,1:2)
    r_anchor(:,3:4) = protein%residue(closing_s(3,i_close))%R(:,2:3)
    rv(:,4,3) = protein%residue(closing_s(3,i_close)+1)%R(:,1)

    call solve_tripep_closure(b_len, b_ang, t_ang, r_anchor, n_soln, soln)
    if (n_soln == 0) cycle

    allocate(s_ang(2,3,n_soln))
    call get_torsion_angles_from_the_solutions(n_soln, soln, rv, s_ang)

    ! pick up the least changed solution in torsions
    do i_soln = 1, n_soln
        d = 0.0d0
        do i_res = 1, 3
            ! delta-phi
            res_no = closing_s(i_res,i_close)
            atm_no = 3
            d = d + abs(bound_angle(protein%residue(res_no)%t_ang(atm_no) - s_ang(1,i_res,i_soln)))
            !
            ! delta-psi
            res_no = res_no + 1
            atm_no = 1
            d = d + abs(bound_angle(protein%residue(res_no)%t_ang(atm_no) - s_ang(2,i_res,i_soln)))
        end do
        !
        if (d < d_min) then
            d_min = d
            closing_min = closing_s(:,i_close)
            s_ang_min = s_ang(:,:,i_soln)
        end if
    end do

    deallocate(soln)
    deallocate(s_ang)
end do

call internal2cartesian_reverse(protein0, closing_min(3), res_j+1)
do i_res = 1, 3
    ! rotate phi
    res_no = closing_min(i_res)
    atm_no = 3
    call rotate_torsion(protein0, res_no, atm_no, s_ang_min(1,i_res))
    !
    ! rotate psi
    res_no = res_no + 1
    atm_no = 1
    call rotate_torsion(protein0, res_no, atm_no, s_ang_min(2,i_res))
end do

call internal2cartesian(protein0, res_i, closing_min(3))
call cartesian2internal(protein0)

end subroutine close_loop
!-------------------------------------------------------------------------------
subroutine close_loop_complete(protein_in, res_i, res_j, n_output, output, unperturbed)
!-------------------------------------------------------------------------------
type(protein_type), intent(in) :: protein_in
integer, intent(in) :: res_i, res_j
integer, intent(out) :: n_output
type(protein_type), intent(out), allocatable :: output(:)
type(protein_type), intent(in), optional :: unperturbed
type(protein_type), allocatable :: protein_tmp(:)
type(protein_type) :: tmp, protein, protein0
!
integer :: n_closing, i_close
integer, allocatable :: closing_s(:,:)
real(dp) :: rv(3,0:4,3), b_len(6), b_ang(7), t_ang(2), r_anchor(3,4)
integer :: i_soln, n_soln, i_res, res_no, atm_no
real(dp), allocatable :: soln(:,:,:,:), s_ang(:,:,:)

protein0 = protein_in
if (present(unperturbed)) then
    call link_broken_bond(protein0, res_j+1, unperturbed=unperturbed)
else
    call link_broken_bond(protein0, res_j+1)
end if

n_output = 0
allocate(output(n_output))

call get_closing_residue_list(protein0, res_i, res_j, n_closing, closing_s)
do i_close = 1, n_closing
    protein = protein0
    call internal2cartesian_reverse(protein, closing_s(3,i_close), res_j+1)
    !
    call get_tripep_geometry(protein0, closing_s(:,i_close), rv, b_len, b_ang, t_ang)
    r_anchor(:,1:2) = protein%residue(closing_s(1,i_close))%R(:,1:2)
    r_anchor(:,3:4) = protein%residue(closing_s(3,i_close))%R(:,2:3)
    rv(:,4,3) = protein%residue(closing_s(3,i_close)+1)%R(:,1)

    call solve_tripep_closure(b_len, b_ang, t_ang, r_anchor, n_soln, soln)
    if (n_soln == 0) cycle

    allocate(s_ang(2,3,n_soln))
    call get_torsion_angles_from_the_solutions(n_soln, soln, rv, s_ang)

    allocate(protein_tmp(n_output+n_soln))
    protein_tmp(1:n_output) = output(1:n_output)
    call move_alloc(protein_tmp, output)

    do i_soln = 1, n_soln
        tmp = protein
        !
        do i_res = 1, 3
            ! rotate phi
            res_no = closing_s(i_res,i_close)
            atm_no = 3
            call rotate_torsion(tmp, res_no, atm_no, s_ang(1,i_res,i_soln))
            !
            ! rotate psi
            res_no = res_no + 1
            atm_no = 1
            call rotate_torsion(tmp, res_no, atm_no, s_ang(2,i_res,i_soln))
        end do
        call internal2cartesian(tmp, res_i, res_j)
        call cartesian2internal(tmp)
        !
        output(n_output+i_soln) = tmp
    end do

    n_output = n_output + n_soln

    if (allocated(soln)) deallocate(soln)
    deallocate(s_ang)
end do

end subroutine close_loop_complete
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
subroutine link_broken_bond(protein, res_no, unperturbed)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein
integer, intent(in) :: res_no
type(protein_type), intent(in), optional :: unperturbed
integer :: i_res
real(dp), parameter :: b_len0 = 1.335d0
real(dp), parameter :: b_ang0(2) = (/116.6d0, 121.9d0/)*deg2rad

if (present(unperturbed)) then
    protein%residue(res_no)%b_len(1) = unperturbed%residue(res_no)%b_len(1)
    protein%residue(res_no)%b_ang(1:2) = unperturbed%residue(res_no)%b_ang(1:2)
    protein%residue(res_no)%t_ang(1:3) = unperturbed%residue(res_no)%t_ang(1:3)
else
    protein%residue(res_no)%b_len(1) = b_len0
    protein%residue(res_no)%b_ang(1:2) = b_ang0
end if

end subroutine link_broken_bond
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END MODULE LOOP_MODELING
!-------------------------------------------------------------------------------
