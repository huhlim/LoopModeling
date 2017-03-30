!-------------------------------------------------------------------------------
MODULE GLOBALS
!-------------------------------------------------------------------------------

implicit none
public

integer, parameter :: dp = kind(1.0d0)
real(dp), parameter :: small_real = 1.0d-10
real(dp), parameter :: pi = 3.141592653589793238462643383279502884197d0

integer :: me
integer, parameter :: king = 0

integer, parameter :: len_fname = 1024
integer, parameter :: len_read = 120
integer, parameter :: len_write = 120

integer :: log_level
character(len=len_write) :: log_msg

integer, parameter :: max_chain = 100
integer, parameter :: max_res = 1000

integer, parameter :: n_stdres = 20
character(len=3), parameter :: stdres(n_stdres) = &
    (/ 'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', &
       'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR' /)
integer, parameter :: n_backbone_atom = 4
character(len=4), parameter :: backbone(n_backbone_atom) = &
    (/ ' N  ', ' CA ', ' C  ', ' O  ' /)

!-------------------------------------------------------------------------------
! Types for protein
!-------------------------------------------------------------------------------
type residue_type
!-------------------------------------------------------------------------------
character(len=1) :: chain_id
character(len=1) :: ter_type
!
integer :: res_index
character(len=3) :: res_name
character(len=5) :: res_no_char
!
integer :: n_atom
character(len=4), allocatable :: atom_name(:)
logical, allocatable :: atom_fixed(:)
!
real(dp), allocatable :: R(:,:)
real(dp), allocatable :: b_len(:)
real(dp), allocatable :: b_ang(:)
real(dp), allocatable :: t_ang(:)
!
integer, allocatable :: atom_prev(:,:,:)
integer, allocatable :: n_torsion_dep(:)
integer, allocatable :: torsion_dep(:,:,:)
!-------------------------------------------------------------------------------
end type residue_type
!-------------------------------------------------------------------------------
type protein_type
!-------------------------------------------------------------------------------
integer :: n_chain
integer, allocatable :: chain_res(:,:)
!
integer :: n_res
type(residue_type), allocatable :: residue(:)
!-------------------------------------------------------------------------------
end type protein_type
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
subroutine allocate_protein_type(protein, n_res, n_chain)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein
integer, intent(in), optional :: n_res
integer, intent(in), optional :: n_chain

if (.not. allocated(protein%residue)) then
    if (present(n_res)) then
        allocate(protein%residue(n_res))
        protein%n_res = n_res
    else
        allocate(protein%residue(max_res))
        protein%n_res = max_res
    end if
end if

if (.not. allocated(protein%chain_res)) then
    if (present(n_chain)) then
        allocate(protein%chain_res(2, n_chain))
        protein%n_chain = n_chain
    else
        allocate(protein%chain_res(2, max_chain))
        protein%n_chain = max_chain
    end if
end if

end subroutine allocate_protein_type
!-------------------------------------------------------------------------------
subroutine deallocate_protein_type(protein)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein

deallocate(protein%residue)
deallocate(protein%chain_res)

protein%n_chain = 0
protein%n_res = 0

end subroutine deallocate_protein_type
!-------------------------------------------------------------------------------
subroutine reallocate_protein_type(protein, n_res, n_chain)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein
integer, intent(in), optional :: n_res
integer, intent(in), optional :: n_chain
!
type(residue_type), allocatable :: residue(:)
integer, allocatable :: chain_res(:,:)

if (present(n_res)) then
    allocate(residue(n_res))
    residue(1:min(n_res, protein%n_res)) = protein%residue(1:min(n_res, protein%n_res))
    call move_alloc(residue, protein%residue)
    protein%n_res = n_res
end if

if (present(n_chain)) then
    allocate(chain_res(2, n_chain))
    chain_res(:,1:min(n_chain, protein%n_chain)) = protein%chain_res(:,1:min(n_chain, protein%n_chain))
    call move_alloc(chain_res, protein%chain_res)
    protein%n_chain = n_chain
end if

end subroutine reallocate_protein_type
!-------------------------------------------------------------------------------
subroutine allocate_residue_type(residue, n_atom)
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: residue
integer, intent(in) :: n_atom

if (allocated(residue%R)) return

allocate(residue%atom_name(n_atom))
allocate(residue%atom_fixed(n_atom))
residue%atom_fixed = .false.

allocate(residue%R(3, n_atom))
allocate(residue%b_len(n_atom))
allocate(residue%b_ang(n_atom))
allocate(residue%t_ang(n_atom))

allocate(residue%atom_prev(2, 3, n_atom))
allocate(residue%n_torsion_dep(n_atom))
allocate(residue%torsion_dep(2, 3, n_atom))

residue%n_atom = 0

end subroutine allocate_residue_type
!-------------------------------------------------------------------------------
subroutine deallocate_residue_type(residue)
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: residue

deallocate(residue%atom_name)
deallocate(residue%atom_fixed)

deallocate(residue%R)
deallocate(residue%b_len)
deallocate(residue%b_ang)
deallocate(residue%t_ang)

deallocate(residue%atom_prev)
deallocate(residue%n_torsion_dep)
deallocate(residue%torsion_dep)

end subroutine deallocate_residue_type
!-------------------------------------------------------------------------------
subroutine reallocate_residue_type(residue, n_atom)
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: residue
integer, intent(in) :: n_atom
!
integer :: n_min
character(len=4), allocatable :: atom_name(:)
logical, allocatable :: atom_fixed(:)
real(dp), allocatable :: R(:,:)
real(dp), allocatable :: b_len(:)
real(dp), allocatable :: b_ang(:)
real(dp), allocatable :: t_ang(:)
integer, allocatable :: atom_prev(:,:,:)
integer, allocatable :: n_torsion_dep(:), torsion_dep(:,:,:)

allocate(atom_name(n_atom))
allocate(atom_fixed(n_atom))
allocate(R(3, n_atom))
allocate(b_len(n_atom))
allocate(b_ang(n_atom))
allocate(t_ang(n_atom))
allocate(atom_prev(2, 3, n_atom))
allocate(n_torsion_dep(n_atom))
allocate(torsion_dep(2, 3, n_atom))

n_min = min(residue%n_atom, n_atom)
atom_name(1:n_min) = residue%atom_name(1:n_min)
atom_fixed(1:n_min) = residue%atom_fixed(1:n_min)
R(:,1:n_min) = residue%R(:,1:n_min)
b_len(1:n_min) = residue%b_len(1:n_min)
b_ang(1:n_min) = residue%b_ang(1:n_min)
t_ang(1:n_min) = residue%t_ang(1:n_min)
atom_prev(:,:,1:n_min) = residue%atom_prev(:,:,1:n_min)
n_torsion_dep(1:n_min) = residue%n_torsion_dep(1:n_min)
torsion_dep(:,:,1:n_min) = residue%torsion_dep(:,:,1:n_min)

call move_alloc(atom_name, residue%atom_name)
call move_alloc(atom_fixed, residue%atom_fixed)
call move_alloc(R, residue%R)
call move_alloc(b_len, residue%b_len)
call move_alloc(b_ang, residue%b_ang)
call move_alloc(t_ang, residue%t_ang)
call move_alloc(atom_prev, residue%atom_prev)
call move_alloc(n_torsion_dep, residue%n_torsion_dep)
call move_alloc(torsion_dep, residue%torsion_dep)

residue%n_atom = n_atom

end subroutine reallocate_residue_type
!-------------------------------------------------------------------------------
END MODULE GLOBALS
!-------------------------------------------------------------------------------
