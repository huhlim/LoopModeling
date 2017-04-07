!-------------------------------------------------------------------------------
PROGRAM LoopModeling
!-------------------------------------------------------------------------------
use globals
use logger
use random, only: initialize_random
use in_out, only: read_pdb, open_write_pdb, close_write_pdb, write_pdb
use geometry, only: cartesian2internal, internal2cartesian, report_protein_geometry
use mathfunction, only: quaternion, rotation_matrix
use loop_modeling, only: close_loop, close_loop_complete

implicit none

integer :: n_argc
character(len=len_fname) :: cmd
character(len=len_fname) :: infile_pdb
character(len=len_fname) :: outfile_pdb

type(protein_type) :: protein, ref
type(protein_type), allocatable :: model(:)

integer :: i,n, n_model
integer, allocatable :: comb(:,:)
integer, parameter :: print_unit = 77

n_argc = iargc()
call getarg(0, cmd)
if (n_argc >= 1) then
    call getarg(1, infile_pdb)
else
    call usage(cmd)
end if

me = 0
call initialize_random()

call read_pdb(infile_pdb, protein)
call cartesian2internal(protein)
call internal2cartesian(protein)

ref = protein
call randomize_torsion(protein)

outfile_pdb = 'out.pdb'
call open_write_pdb(print_unit, outfile_pdb)

call close_loop_complete(protein, 47, 52, n_model, model, unperturbed=ref)
do i = 1, n_model
    call write_pdb(print_unit, model(i), i)
end do
!
!call close_loop(protein, 47, 52, unperturbed=ref)
!call write_pdb(print_unit, protein)

call close_write_pdb(print_unit)

!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
subroutine usage(cmd)
!-------------------------------------------------------------------------------
character(len=len_fname) :: cmd

write(log_msg, '(A,1x,A)') trim(cmd), '[Input PDB]'
call terminate_with_error(log_msg)

end subroutine usage
!-------------------------------------------------------------------------------
subroutine randomize_torsion(protein)
!-------------------------------------------------------------------------------
type(protein_type), intent(inout) :: protein
integer :: i_res

do i_res = 48, 51
    protein%residue(i_res)%t_ang(1) = &
        protein%residue(i_res)%t_ang(1) - 0.5d0
    protein%residue(i_res)%t_ang(3) = &
        protein%residue(i_res)%t_ang(3) - 0.5d0
end do

call internal2cartesian(protein, 47, 52)
call cartesian2internal(protein)

end subroutine randomize_torsion
!-------------------------------------------------------------------------------
END PROGRAM LoopModeling
!-------------------------------------------------------------------------------
