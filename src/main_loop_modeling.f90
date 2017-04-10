!-------------------------------------------------------------------------------
PROGRAM LoopModeling
!-------------------------------------------------------------------------------
use globals
use logger
use random, only: initialize_random
use in_out, only: read_pdb, open_write_pdb, close_write_pdb, write_pdb
use geometry, only: cartesian2internal, internal2cartesian, calc_rmsd_CA, calc_rmsd
use mathfunction, only: quaternion, rotation_matrix
use loop_modeling, only: close_loop, close_loop_complete

use cluster, only: hierarchical_clustering, hierarchical_clustering_cutoff

implicit none

integer :: n_argc
character(len=len_fname) :: cmd
character(len=len_fname) :: infile_pdb, native_pdb
character(len=len_fname) :: outfile_pdb

type(protein_type) :: protein, ref, native
type(protein_type) :: loop, ref_loop, native_loop
type(protein_type), allocatable :: model(:)

integer :: i, j, n_model, n_sum, res_i, res_j, n_cluster
integer, allocatable :: comb(:,:)
integer, parameter :: print_unit = 77
real(dp), allocatable :: dmtx(:,:)
integer, allocatable :: center(:)
integer, parameter :: res_start=47, res_end=52
integer, parameter :: res_start=44, res_end=55

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

native_pdb = 'ref.pdb'
call read_pdb(native_pdb, native)
call cartesian2internal(native)
call internal2cartesian(native)

ref = protein
call extract_loop_from_protein(ref, ref_loop, res_start, res_end, res_i, res_j)
call extract_loop_from_protein(native, native_loop, res_start, res_end, res_i, res_j)
call extract_loop_from_protein(protein, loop, res_start, res_end, res_i, res_j)

outfile_pdb = 'output.pdb'
call open_write_pdb(print_unit, outfile_pdb)

call close_loop_complete(loop, res_i, res_j, n_model, model, unperturbed=ref_loop)

allocate(dmtx(n_model, n_model))
dmtx = 0.0d0
do i = 1, n_model
    do j = 1, i-1
        dmtx(j,i) = calc_rmsd_CA(model(i), model(j), res_i, res_j)
    end do
end do

dmtx = dmtx + transpose(dmtx)

!call hierarchical_clustering(n_model, 10, dmtx, center)
call hierarchical_clustering_cutoff(n_model, 0.5d0, dmtx, n_cluster, center)
deallocate(dmtx)

do i = 1, n_cluster
    write(*,'(I4, 2x,F7.3)') i, calc_rmsd(native_loop, model(center(i)), res_i, res_j)
    call copy_loop_to_protein(model(center(i)), protein, res_start, res_end)
    call write_pdb(print_unit, protein, i)
end do
!do i = 1, n_model
!    write(*,'(I4, 2x,F7.3)') i, calc_rmsd(native_loop, model(i), res_i, res_j)
!    call copy_loop_to_protein(model(i), protein, res_start, res_end)
!    call write_pdb(print_unit, protein, i)
!end do
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
