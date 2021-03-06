!-------------------------------------------------------------------------------
MODULE GEOMETRY
!-------------------------------------------------------------------------------
use globals
use mathfunction

implicit none
public

! indices
integer, parameter :: index_phi = 3 
integer, parameter :: index_psi = 1
integer, parameter :: index_omg = 2

CONTAINS
!-------------------------------------------------------------------------------
subroutine cartesian2internal(protein, res_i, res_j)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein
integer, intent(in), optional :: res_i, res_j
integer :: res_start, res_end
integer :: i_res, i_atm, ia, j_res, j_atm
real(dp) :: R(3,4)
logical :: defined(3)

if (present(res_i)) then
    res_start = res_i
else
    res_start = 1
end if

if (present(res_j)) then
    res_end = res_j
else
    res_end = protein%n_res
end if

do i_res = res_start, res_end
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
subroutine internal2cartesian(protein, res_i, res_j)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein
integer, intent(in), optional :: res_i, res_j
integer :: res_start, res_end
integer :: i_res, i_atm, ia, j_res, j_atm
real(dp) :: R(3,3), b_len, b_ang, t_ang
logical :: defined

if (present(res_i)) then
    res_start = res_i
else
    res_start = 1
end if

if (present(res_j)) then
    res_end = res_j
else
    res_end = protein%n_res
end if

do i_res = res_start, res_end
    do i_atm = 1, protein%residue(i_res)%n_atom
        defined = .true.
        do ia = 1, 3
            j_res = protein%residue(i_res)%atom_prev(1,ia,i_atm)
            j_atm = protein%residue(i_res)%atom_prev(2,ia,i_atm)
            if (j_res > 0) then
                R(:,ia) = protein%residue(j_res)%R(:,j_atm)
            else
                defined = .false.
                exit
            end if
        end do
        if (.not. defined) cycle
        !
        b_len = protein%residue(i_res)%b_len(i_atm)
        b_ang = protein%residue(i_res)%b_ang(i_atm)
        t_ang = protein%residue(i_res)%t_ang(i_atm)
        !
        protein%residue(i_res)%R(:,i_atm) = place_an_atom(R, b_len, b_ang, t_ang)
    end do
end do

end subroutine internal2cartesian
!-------------------------------------------------------------------------------
subroutine internal2cartesian_reverse(protein, res_i, res_j)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein
integer, intent(in), optional :: res_i, res_j
integer :: res_start, res_end
integer :: i_res, i_atm, ia, j_res, j_atm, atom_prev(2,3)
real(dp) :: R(3,3), par(3), b_len, b_ang, t_ang
logical :: defined, localized

if (present(res_i)) then
    res_start = res_i
else
    res_start = 1
end if

if (present(res_j)) then
    res_end = res_j
else
    res_end = protein%n_res
end if

do i_res = res_end, res_start, -1
    do i_atm = protein%residue(i_res)%n_atom, 1, -1
        localized = .true.
        do ia = 1, 3
            if (protein%residue(i_res)%atom_prev(1,ia,i_atm) /= i_res) then
                localized = .false.
                exit
            end if
        end do
        if (localized) cycle
        !
        call find_atom_prev_reverse(protein, i_res, i_atm, atom_prev, par)
        !
        defined = .true.
        do ia = 1, 3
            j_res = atom_prev(1,ia)
            j_atm = atom_prev(2,ia)
            if (j_res > 0) then
                R(:,ia) = protein%residue(j_res)%R(:,j_atm)
            else
                defined = .false.
                exit
            end if
        end do
        if (.not. defined) cycle
        !
        b_len = par(1)
        b_ang = par(2)
        t_ang = par(3)
        !
        protein%residue(i_res)%R(:,i_atm) = place_an_atom(R, b_len, b_ang, t_ang)
    end do
    !
    do i_atm = 1, protein%residue(i_res)%n_atom
        localized = .true.
        do ia = 1, 3
            if (protein%residue(i_res)%atom_prev(1,ia,i_atm) /= i_res) then
                localized = .false.
                exit
            end if
        end do
        if (.not. localized) cycle
        !
        defined = .true.
        do ia = 1, 3
            j_res = protein%residue(i_res)%atom_prev(1,ia,i_atm)
            j_atm = protein%residue(i_res)%atom_prev(2,ia,i_atm)
            if (j_res > 0) then
                R(:,ia) = protein%residue(j_res)%R(:,j_atm)
            else
                defined = .false.
                exit
            end if
        end do
        if (.not. defined) cycle
        !
        b_len = protein%residue(i_res)%b_len(i_atm)
        b_ang = protein%residue(i_res)%b_ang(i_atm)
        t_ang = protein%residue(i_res)%t_ang(i_atm)
        !
        protein%residue(i_res)%R(:,i_atm) = place_an_atom(R, b_len, b_ang, t_ang)
    end do
end do

end subroutine internal2cartesian_reverse
!-------------------------------------------------------------------------------
subroutine report_protein_geometry(protein, res_i, res_j)
!-------------------------------------------------------------------------------
type(protein_type), intent(in) :: protein
integer, intent(in), optional :: res_i, res_j
integer :: i_res, i_atm, res_start, res_end

if (present(res_i)) then
    res_start = res_i
else
    res_start = 1
end if
if (present(res_j)) then
    res_end = res_j
else
    res_end = protein%n_res
end if

do i_res = res_start, res_end
    do i_atm = 1, protein%residue(i_res)%n_atom
        write(*,'(A,2x,I4,1x,I2,2x,3F8.3,2x,F8.3,2(1x,F7.1))') &
            'GEOMETRY', i_res, i_atm, &
            protein%residue(i_res)%R(:,i_atm), &
            protein%residue(i_res)%b_len(i_atm), &
            protein%residue(i_res)%b_ang(i_atm)*rad2deg, &
            protein%residue(i_res)%t_ang(i_atm)*rad2deg
    end do
end do

end subroutine report_protein_geometry
!-------------------------------------------------------------------------------
function calc_bond(u, v)
!-------------------------------------------------------------------------------
real(dp) :: calc_bond
real(dp), intent(in) :: u(3), v(3)

calc_bond = v_size(v-u)

end function calc_bond
!-------------------------------------------------------------------------------
function bond_length(R)
!-------------------------------------------------------------------------------
real(dp) :: bond_length
real(dp), intent(in) :: R(3,2)

bond_length = calc_bond(R(:,1), R(:,2))

end function bond_length
!-------------------------------------------------------------------------------
function calc_angle(u, v)
!-------------------------------------------------------------------------------
real(dp) :: calc_angle
real(dp), intent(in) :: u(3), v(3)

calc_angle = dacos(dot_product(v_norm(u), v_norm(v)))

end function calc_angle
!-------------------------------------------------------------------------------
function bond_angle(R)
!-------------------------------------------------------------------------------
real(dp) :: bond_angle
real(dp), intent(in) :: R(3,3)

bond_angle = calc_angle(R(:,1)-R(:,2), R(:,3)-R(:,2))

end function bond_angle
!-------------------------------------------------------------------------------
function calc_torsion(u, v, w)
!-------------------------------------------------------------------------------
real(dp) :: calc_torsion
real(dp), intent(in) :: u(3), v(3), w(3)
real(dp) :: p(3), q(3), s(3), arg

p = v_norm(cross(u, v))
q = v_norm(cross(v, w))
s = cross(w, u)

arg = dot_product(p, q)
arg = min(1.0d0, max(-1.0d0, arg))

calc_torsion = sign(dacos(arg), dot_product(s, v))

end function calc_torsion
!-------------------------------------------------------------------------------
function torsion_angle(R)
!-------------------------------------------------------------------------------
real(dp) :: torsion_angle
real(dp), intent(in) :: R(3,4)
real(dp) :: dR(3,3)

dR(:,:) = R(:,2:4) - R(:,1:3)
torsion_angle = calc_torsion(dR(:,1), dR(:,2), dR(:,3))

end function torsion_angle
!-------------------------------------------------------------------------------
subroutine rotate_torsion(protein, i_res, i_atm, t_ang)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein
integer, intent(in) :: i_res, i_atm
real(dp), intent(in) :: t_ang
real(dp) :: d_ang, ang
integer :: i_dep, j_res, j_atm

d_ang = t_ang - protein%residue(i_res)%t_ang(i_atm)

do i_dep = 1, protein%residue(i_res)%n_torsion_dep(i_atm)
    j_res = protein%residue(i_res)%torsion_dep(1, i_dep, i_atm)
    j_atm = protein%residue(i_res)%torsion_dep(2, i_dep, i_atm)
    ang = protein%residue(j_res)%t_ang(j_atm) 
    protein%residue(j_res)%t_ang(j_atm) = &
        bound_angle(protein%residue(j_res)%t_ang(j_atm) + d_ang)
    !write(*,'(A,3(1x,I4), 3(2x,F6.1))') 'REMARK', i_dep, j_res, j_atm, ang*rad2deg, &
    !        protein%residue(j_res)%t_ang(j_atm) *rad2deg
end do

end subroutine rotate_torsion
!-------------------------------------------------------------------------------
function place_an_atom(R, b_len, b_ang, t_ang)
!-------------------------------------------------------------------------------
real(dp) :: place_an_atom(3)
real(dp), intent(in) :: R(3,3)
real(dp), intent(in) :: b_len, b_ang, t_ang
real(dp) :: b(3), q1(4), q2(4), q(4), axis(3), U(3,3)

b = R(:,1) - R(:,2)

axis = cross(b, R(:,3)-R(:,2))
q1 = quaternion(axis, pi-b_ang)
q2 = quaternion(b, t_ang)
q = q_product(q2, q1)

U = rotation_matrix(q)
place_an_atom = R(:,1) + rotate(v_norm(b), U) * b_len

end function place_an_atom
!-------------------------------------------------------------------------------
subroutine find_atom_prev_reverse(protein, i_res, i_atm, atom_prev, par)
!-------------------------------------------------------------------------------
type(protein_type), intent(in) :: protein
integer, intent(in) :: i_res, i_atm
integer, intent(out) :: atom_prev(2,3)
real(dp), intent(out) :: par(3)
integer :: j_res, j_atm

atom_prev = 0
par = 0.0d0

do j_res = i_res, min(i_res+1, protein%n_res)
    do j_atm = 1, protein%residue(j_res)%n_atom
        if (protein%residue(j_res)%atom_prev(1,3,j_atm) == i_res .and. &
            protein%residue(j_res)%atom_prev(2,3,j_atm) == i_atm) then
            atom_prev(:,1) = protein%residue(j_res)%atom_prev(:,2,j_atm)
            atom_prev(:,2) = protein%residue(j_res)%atom_prev(:,1,j_atm)
            atom_prev(:,3) = (/j_res, j_atm/)
            !
            par(1) = protein%residue(atom_prev(1,1))%b_len(atom_prev(2,1))
            par(2) = protein%residue(atom_prev(1,2))%b_ang(atom_prev(2,2))
            par(3) = protein%residue(atom_prev(1,3))%t_ang(atom_prev(2,3))
        end if
    end do
end do

end subroutine find_atom_prev_reverse
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
function calc_rmsd(ref, protein, res_i, res_j)
!-------------------------------------------------------------------------------
real(dp) :: calc_rmsd
type(protein_type), intent(in) :: ref, protein
integer, intent(in) :: res_i, res_j
integer :: i_res, i_atm
real(dp) :: dr(3,4)

calc_rmsd = 0.0d0
do i_res = res_i, res_j
    dr = protein%residue(i_res)%R - ref%residue(i_res)%R
    do i_atm = 1, 4
        calc_rmsd = calc_rmsd + dot_product(dr(:,i_atm),dr(:,i_atm))
    end do
end do

calc_rmsd = dsqrt(calc_rmsd / (4.0d0 * (res_j-res_i+1)))

end function calc_rmsd
!-------------------------------------------------------------------------------
function calc_rmsd_CA(ref, protein, res_i, res_j)
!-------------------------------------------------------------------------------
real(dp) :: calc_rmsd_CA
type(protein_type), intent(in) :: ref, protein
integer, intent(in) :: res_i, res_j
integer :: i_res
real(dp) :: dr(3)

calc_rmsd_CA = 0.0d0
do i_res = res_i, res_j
    dr = protein%residue(i_res)%R(:,2) - ref%residue(i_res)%R(:,2)
    calc_rmsd_CA = calc_rmsd_CA + dot_product(dr,dr)
end do

calc_rmsd_CA = dsqrt(calc_rmsd_CA / (res_j-res_i+1))

end function calc_rmsd_CA
!-------------------------------------------------------------------------------
END MODULE GEOMETRY
!-------------------------------------------------------------------------------
