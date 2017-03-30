!-------------------------------------------------------------------------------
MODULE GEOMETRY
!-------------------------------------------------------------------------------
use globals
use mathfunction

implicit none
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine cartesian2internal(protein, res_i, res_j)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein
integer, intent(in) :: res_i, res_j
integer :: i_res, i_atm, ia, j_res, j_atm
real(dp) :: R(3,4)
logical :: defined(3)

do i_res = res_i, res_j
    do i_atm = 1, protein%residue(i_res)%n_atom
        R(:,1) = protein%residue(i_res)%R(:,i_atm)
        !
        defined = .false.
        do ia = 1, 3
            j_res = protein%residue(i_res)%atom_prev(1,ia,i_atm)
            j_atm = protein%residue(i_res)%atom_prev(2,ia,i_atm)
            if (j_res > 0) then
                R(:,ia+1) = protein%residue(j_res)%R(:,j_atm)
                defined(ia) = .true.
            end if
        end do
        !
        if (.not. defined(1)) cycle
        protein%residue(i_res)%b_len(i_atm) = bond_length(R(:,1:2))
        !
        if (.not. defined(2)) cycle
        protein%residue(i_res)%b_ang(i_atm) = bond_angle(R(:,1:3))
        !
        if (.not. defined(3)) cycle
        protein%residue(i_res)%t_ang(i_atm) = torsion_angle(R(:,1:4))
    end do
end do

end subroutine cartesian2internal
!-------------------------------------------------------------------------------
subroutine internal2cartesian(residue, res_i, res_j)
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: residue(:)
integer, intent(in) :: res_i, res_j

end subroutine internal2cartesian
!-------------------------------------------------------------------------------
function bond_length(R)
!-------------------------------------------------------------------------------
real(dp) :: bond_length
real(dp), intent(in) :: R(3,2)

bond_length = v_size(R(:,2) - R(:,1))

end function bond_length
!-------------------------------------------------------------------------------
function bond_angle(R)
!-------------------------------------------------------------------------------
real(dp) :: bond_angle
real(dp), intent(in) :: R(3,3)
real(dp) :: dR(3,2)
real(dp) :: r12(3), r32(3)

dR(:,1) = v_norm(R(:,1) - R(:,2))
dR(:,2) = v_norm(R(:,3) - R(:,2))

bond_angle = dacos(dot_product(dR(:,1),dR(:,2)))

end function bond_angle
!-------------------------------------------------------------------------------
function torsion_angle(R)
!-------------------------------------------------------------------------------
real(dp) :: torsion_angle
real(dp), intent(in) :: R(3,4)
real(dp) :: dR(3,3), p(3), q(3), s(3)
real(dp) :: arg

dR(:,:) = R(:,2:4) - R(:,1:3)

p = v_norm(cross(dR(:,1), dR(:,2)))
q = v_norm(cross(dR(:,2), dR(:,3)))
s = cross(dR(:,3), dR(:,1))

arg = dot_product(p,q)
arg = min(1.0d0, max(-1.0d0, arg))

torsion_angle = sign(dacos(arg), dot_product(s, dR(:,2)))

end function torsion_angle
!-------------------------------------------------------------------------------
subroutine rotate_torsion(protein, i_res, i_atm, t_ang)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein
integer, intent(in) :: i_res, i_atm
real(dp), intent(in) :: t_ang

end subroutine rotate_torsion
!-------------------------------------------------------------------------------
subroutine initialize_protein_geometry(protein)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein
integer :: i_res, i_atm, ia, j_res, j_atm, i_dep
integer :: bb_geom(2,3,4)
logical :: has_same_prev_atoms

bb_geom = reshape( (/ -1,3, -1,2, -1,1, &
                       0,1, -1,3, -1,2, &
                       0,2,  0,1, -1,3, &
                       0,3,  0,2,  0,1 /), (/2,3,4/))

! find indices for internal coordinate definiations
do i_res = 1, protein%n_res
    do i_atm = 1, protein%residue(i_res)%n_atom
        do ia = 1, 3
            protein%residue(i_res)%atom_prev(1,ia,i_atm) = i_res + bb_geom(1,ia,i_atm)
            protein%residue(i_res)%atom_prev(2,ia,i_atm) = bb_geom(2,ia,i_atm)
        end do
    end do
end do

! find indices for torsion angle dependencies
do i_res = 1, protein%n_res
    do i_atm = 1, protein%residue(i_res)%n_atom
        i_dep = 0
        do j_res = max(1, i_res-1), min(protein%n_res, i_res+1)
            do j_atm = 1, protein%residue(j_res)%n_atom
                has_same_prev_atoms = .true.
                do ia = 1, 3
                    if (protein%residue(i_res)%atom_prev(1,ia,i_atm) /= &
                        protein%residue(j_res)%atom_prev(1,ia,j_atm)) then
                        has_same_prev_atoms = .false.
                        exit
                    end if
                    if (protein%residue(i_res)%atom_prev(2,ia,i_atm) /= &
                        protein%residue(j_res)%atom_prev(2,ia,j_atm)) then
                        has_same_prev_atoms = .false.
                        exit
                    end if
                end do
                if (has_same_prev_atoms) then
                    i_dep = i_dep + 1
                    protein%residue(i_res)%torsion_dep(1:2,i_dep,i_atm) = (/j_res,j_atm/)
                end if
            end do
        end do
        !
        protein%residue(i_res)%n_torsion_dep(i_atm) = i_dep
    end do
end do

do i_res = 1, protein%n_res
    protein%residue(i_res)%b_len(:) = 0.0d0
    protein%residue(i_res)%b_ang(:) = 0.0d0
    protein%residue(i_res)%t_ang(:) = 0.0d0
end do

end subroutine initialize_protein_geometry
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END MODULE GEOMETRY
!-------------------------------------------------------------------------------
