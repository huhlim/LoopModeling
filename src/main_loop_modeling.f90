!-------------------------------------------------------------------------------
PROGRAM LoopModeling
!-------------------------------------------------------------------------------
use globals
use logger
use random, only: initialize_random
use in_out, only: read_pdb, open_write_pdb, close_write_pdb, write_pdb
use geometry, only: cartesian2internal, internal2cartesian, calc_RMSD_CA
use mathfunction, only: quaternion, rotation_matrix
use loop_modeling, only: close_loop, close_loop_complete

use cluster, only: hierarchical_clustering

implicit none

integer :: n_argc
character(len=len_fname) :: cmd
character(len=len_fname) :: infile_pdb
character(len=len_fname) :: outfile_pdb

type(protein_type) :: protein, ref
type(protein_type) :: loop, ref_loop
type(protein_type), allocatable :: model(:)

integer :: i, j, n_model, n_sum, res_i, res_j
integer, allocatable :: comb(:,:)
integer, parameter :: print_unit = 77
real(dp), allocatable :: dmtx(:,:)
integer :: center(10)

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
call extract_loop_from_protein(ref, ref_loop, 47, 52, res_i, res_j)
call randomize_torsion(protein)
call extract_loop_from_protein(protein, loop, 47, 52, res_i, res_j)

outfile_pdb = 'out.pdb'
call open_write_pdb(print_unit, outfile_pdb)

call write_pdb(print_unit, ref, 0)

call close_loop_complete(loop, res_i, res_j, n_model, model, unperturbed=ref_loop)

allocate(dmtx(n_model, n_model))
dmtx = 0.0d0
do i = 1, n_model-1
    do j = i+1, n_model
        dmtx(j,i) = calc_RMSD_CA(model(i), model(j), res_i, res_j)
    end do
end do
dmtx = dmtx + transpose(dmtx)

call hierarchical_clustering(n_model, 10, dmtx, center)

do i = 1, n_model
    call copy_loop_to_protein(model(i), protein, 47, 52)
    call write_pdb(print_unit, protein, i)
end do
deallocate(model)

call close_write_pdb(print_unit)

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
