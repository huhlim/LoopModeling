!-------------------------------------------------------------------------------
MODULE IN_OUT
!-------------------------------------------------------------------------------
use globals
use logger
use geometry, only: initialize_protein_geometry

implicit none
private

public :: read_pdb

public :: open_write_pdb
public :: close_write_pdb
public :: write_pdb

CONTAINS
!-------------------------------------------------------------------------------
subroutine read_pdb(infile_pdb, protein)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_pdb
type(protein_type), intent(out) :: protein
!
integer, parameter :: f_unit = 20
integer :: ioerr
character(len=len_read) :: line
!
integer :: res_no, chain_no
character(len=1) :: chain_id, chain_id_prev
character(len=5) :: res_no_char, res_no_char_prev
!
character(len=4) :: res_name, atom_name
integer :: i_res, i_atm

call allocate_protein_type(protein)

open(f_unit, file=trim(infile_pdb), iostat=ioerr)
if (ioerr < 0) then
    call terminate_with_error("Error: failed to open the input PDB file")
end if

chain_id_prev = '-'
res_no_char_prev = ''
res_no = 0
chain_no = 0

do
    read(f_unit, '(A120)', iostat=ioerr) line
    if (ioerr < 0) exit
    !
    if (line(1:4) /= 'ATOM') cycle
    !
    atom_name = line(13:16)
    i_atm = backbone_atom_index(atom_name)
    if (i_atm == -1) cycle ! read backbone atoms only
    !
    res_name = line(18:20)
    i_res = residue_index(res_name)
    if (i_res == -1) cycle  ! read standard residues only
    !
    res_no_char = line(23:27)
    chain_id = line(22:22)
    if ((res_no_char /= res_no_char_prev) .or. (chain_id /= chain_id_prev)) then
        ! new residue
        res_no = res_no + 1
        if (res_no > protein%n_res) then
            call reallocate_protein_type(protein, n_res=protein%n_res+max_res)
        end if
        res_no_char_prev = res_no_char
        !
        protein%residue(res_no)%chain_id = chain_id
        protein%residue(res_no)%res_index = i_res
        protein%residue(res_no)%res_name = res_name
        protein%residue(res_no)%res_no_char = res_no_char
        call allocate_residue_type(protein%residue(res_no), n_backbone_atom)
        !
        if (chain_id /= chain_id_prev) then
            chain_no = chain_no + 1
            if (chain_no > protein%n_chain) then
                call reallocate_protein_type(protein, n_chain=protein%n_chain+max_chain)
            end if
            chain_id_prev = chain_id
            !
            protein%chain_res(1, chain_no) = res_no
        end if
        protein%chain_res(2, chain_no) = res_no
    end if
    !
    protein%residue(res_no)%atom_name(i_atm) = atom_name
    protein%residue(res_no)%n_atom = protein%residue(res_no)%n_atom + 1
    read(line(31:54),'(3F8.3)') protein%residue(res_no)%R(1:3,i_atm)
end do

close(f_unit)

call reallocate_protein_type(protein, n_res=res_no, n_chain=chain_no)

do res_no = 1, protein%n_res
    protein%residue(res_no)%ter_type = ''
end do
do chain_no = 1, protein%n_chain
    protein%residue(protein%chain_res(1,chain_no))%ter_type = 'N'
    protein%residue(protein%chain_res(2,chain_no))%ter_type = 'C'
end do

call initialize_protein_geometry(protein)

!-------------------------------------------------------------------------------
end subroutine read_pdb
!-------------------------------------------------------------------------------
subroutine open_write_pdb(f_unit, outfile_pdb)
!-------------------------------------------------------------------------------
integer, intent(in) :: f_unit
character(len=len_fname), intent(in) :: outfile_pdb
integer :: ioerr

open(f_unit, file=trim(outfile_pdb), iostat=ioerr, status='replace')
if (ioerr < 0) then
    call terminate_with_error("Error: failed to open the output PDB file")
end if

end subroutine open_write_pdb
!-------------------------------------------------------------------------------
subroutine close_write_pdb(f_unit)
!-------------------------------------------------------------------------------
integer, intent(in) :: f_unit

write(f_unit, '(A)') "END"
if (f_unit /= 6) close(f_unit)

end subroutine close_write_pdb
!-------------------------------------------------------------------------------
subroutine write_pdb(f_unit, protein, model_no)
!-------------------------------------------------------------------------------
integer, intent(in) :: f_unit
type(protein_type), intent(in) :: protein
integer, intent(in), optional :: model_no
integer :: res_no, atom_no, ia

if (present(model_no)) then
    write(f_unit, '(A,1x,I5)') 'MODEL', model_no
end if

ia = 0
do res_no = 1, protein%n_res
    do atom_no = 1, protein%residue(res_no)%n_atom
        ia = ia + 1
        write(f_unit, 99) 'ATOM  ', ia, &
            protein%residue(res_no)%atom_name(atom_no), &
            protein%residue(res_no)%res_name, &
            protein%residue(res_no)%chain_id, &
            protein%residue(res_no)%res_no_char, &
            protein%residue(res_no)%R(1:3, atom_no), &
            1.0, 0.0
    end do
    !
    if (protein%residue(res_no)%ter_type == 'C') then
        write(f_unit, '(A)') "TER"
    end if
end do

if (present(model_no)) then
    write(f_unit, '(A)') 'ENDMDL'
end if

99 format (A6,I5,1x,A4,1x,A3,1x,A1,A5,3x,3F8.3,2x,F4.2,F6.2)

end subroutine write_pdb
!-------------------------------------------------------------------------------
function residue_index(res_name)
!-------------------------------------------------------------------------------
integer :: residue_index
character(len=4), intent(in) :: res_name
integer :: i_res

residue_index = -1
do i_res = 1, n_stdres
    if (trim(res_name) == trim(stdres(i_res))) then
        residue_index = i_res
        return
    end if
end do

end function residue_index
!-------------------------------------------------------------------------------
function backbone_atom_index(atom_name)
!-------------------------------------------------------------------------------
integer :: backbone_atom_index
character(len=4), intent(in) :: atom_name
integer :: i_atom

backbone_atom_index = -1
do i_atom = 1, n_backbone_atom
    if (trim(atom_name) == trim(backbone(i_atom))) then
        backbone_atom_index = i_atom
        return
    end if
end do

end function backbone_atom_index
!-------------------------------------------------------------------------------
END MODULE IN_OUT
!-------------------------------------------------------------------------------
