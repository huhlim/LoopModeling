!-------------------------------------------------------------------------------
PROGRAM LoopModeling
!-------------------------------------------------------------------------------
use globals
use logger
use in_out, only: read_pdb, open_write_pdb, close_write_pdb, write_pdb
use geometry, only: cartesian2internal, internal2cartesian
use mathfunction, only: quaternion, rotation_matrix
use loop_modeling, only: close_loop

implicit none

integer :: n_argc
character(len=len_fname) :: cmd
character(len=len_fname) :: infile_pdb

type(protein_type) :: protein

integer :: i,n
integer, allocatable :: comb(:,:)

n_argc = iargc()
call getarg(0, cmd)
if (n_argc >= 1) then
    call getarg(1, infile_pdb)
else
    call usage(cmd)
end if

me = 0

call read_pdb(infile_pdb, protein)
call cartesian2internal(protein)
call internal2cartesian(protein)

call close_loop(protein, 48, 50)

call write_pdb(print_screen, protein)

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
