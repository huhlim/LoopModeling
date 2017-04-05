!-------------------------------------------------------------------------------
PROGRAM LoopModeling
!-------------------------------------------------------------------------------
use globals
use logger
use random, only: initialize_random
use in_out, only: read_pdb, open_write_pdb, close_write_pdb, write_pdb
use geometry, only: cartesian2internal, internal2cartesian
use mathfunction, only: quaternion, rotation_matrix
use loop_modeling, only: close_loop, close_loop_complete

implicit none

integer :: n_argc
character(len=len_fname) :: cmd
character(len=len_fname) :: infile_pdb

type(protein_type) :: protein
type(protein_type), allocatable :: model(:)

integer :: i,n, n_model
integer, allocatable :: comb(:,:)

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

!call close_loop(protein, 46, 52)
call close_loop_complete(protein, 46, 52, n_model, model)
print*, n_model

do i = 1, n_model
    call write_pdb(print_screen, model(i), i)
end do

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
END PROGRAM LoopModeling
!-------------------------------------------------------------------------------
