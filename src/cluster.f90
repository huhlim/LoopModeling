!-------------------------------------------------------------------------------
MODULE CLUSTER
!-------------------------------------------------------------------------------
use globals

implicit none
private

real(dp), parameter :: maximum_distance = 9999999999.9d0

public :: hierarchical_clustering
public :: hierarchical_clustering_cutoff

CONTAINS
!-------------------------------------------------------------------------------
subroutine hierarchical_clustering(n_conf0, n_cluster, dmatrix, center)
!-------------------------------------------------------------------------------
integer, intent(in) :: n_conf0, n_cluster
real(dp), intent(in) :: dmatrix(:,:)
integer, intent(out) :: center(n_cluster)
integer :: i, k, n_conf, clij(2)
real(dp) :: dij
integer, allocatable :: memb(:,:), clsize(:), i_conf(:)
real(dp), allocatable :: dmtx(:,:)

call remove_redundant(n_conf0, dmatrix, n_conf, dmtx, i_conf)
allocate(memb(n_conf, n_conf))
allocate(clsize(n_conf))
memb(:,:) = 0
clsize(:) = 1
do i = 1, n_conf
    memb(1,i) = i
end do

do i = n_conf, n_cluster+1, -1
    !call linkage_maximum(dmtx, n_conf, memb, clsize, clij, dij)
    call linkage_minimum(dmtx, n_conf, memb, clsize, clij, dij)
    !call linkage_average(dmtx, n_conf, memb, clsize, clij, dij)
    !call linkage_centroid(dmtx, n_conf, memb, clsize, clij, dij)
    ! for linkage_ward
    call update_ward_distance(dmtx, n_conf, clij(1), clij(1), clij(2), memb, clsize)
    !
    memb(clsize(clij(1))+1:clsize(clij(1))+clsize(clij(2)), clij(1)) = memb(1:clsize(clij(2)), clij(2))
    clsize(clij(1)) = clsize(clij(1)) + clsize(clij(2))
    clsize(clij(2)) = 0
    !
    ! for linkage_centroid
    !call get_centroid(dmtx, clsize(clij(1)), memb(:,clij(1)))
end do

k = 0
do i = 1, n_cluster
    if (clsize(i) == 0) cycle
    !
    k = k + 1
    call get_centroid(dmtx, clsize(i), memb(:,i))
    center(k) = i_conf(memb(1,i))
end do

deallocate(memb)
deallocate(clsize)

end subroutine hierarchical_clustering
!-------------------------------------------------------------------------------
subroutine hierarchical_clustering_cutoff(n_conf0, cutoff, dmatrix, n_cluster, center)
!-------------------------------------------------------------------------------
integer, intent(in) :: n_conf0
real(dp), intent(in) :: dmatrix(:,:), cutoff
integer, intent(out) :: n_cluster
integer, intent(out), allocatable :: center(:)
integer :: i, k, n_conf, clij(2)
real(dp) :: dij
integer, allocatable :: memb(:,:), clsize(:), i_conf(:)
real(dp), allocatable :: dmtx(:,:)

call remove_redundant(n_conf0, dmatrix, n_conf, dmtx, i_conf)
allocate(memb(n_conf, n_conf))
allocate(clsize(n_conf))
memb(:,:) = 0
clsize(:) = 1
do i = 1, n_conf
    memb(1,i) = i
end do

n_cluster = n_conf
do i = n_conf, 1, -1
    !call linkage_maximum(dmtx, n_conf, memb, clsize, clij, dij)
    call linkage_minimum(dmtx, n_conf, memb, clsize, clij, dij)
    !call linkage_average(dmtx, n_conf, memb, clsize, clij, dij)
    !call linkage_centroid(dmtx, n_conf, memb, clsize, clij, dij)
    write(*,'(3(I5,2x),F8.5)') n_cluster, clij, dij
    if (dij > cutoff) exit
    !
    ! for linkage_ward
    call update_ward_distance(dmtx, n_conf, clij(1), clij(1), clij(2), memb, clsize)
    !
    n_cluster = n_cluster - 1
    memb(clsize(clij(1))+1:clsize(clij(1))+clsize(clij(2)), clij(1)) = memb(1:clsize(clij(2)), clij(2))
    clsize(clij(1)) = clsize(clij(1)) + clsize(clij(2))
    clsize(clij(2)) = 0
    !
    ! for linkage_centroid
    !call get_centroid(dmtx, clsize(clij(1)), memb(:,clij(1)))
end do

allocate(center(n_cluster))

k = 0
do i = 1, n_conf
    if (clsize(i) == 0) cycle
    !
    k = k + 1
    call get_centroid(dmtx, clsize(i), memb(:,i))
    center(k) = i_conf(memb(1,i))
end do

deallocate(memb)
deallocate(clsize)

end subroutine hierarchical_clustering_cutoff
!-------------------------------------------------------------------------------
subroutine linkage_maximum(dmtx, n_conf, memb, clsize, ij_min, d_min)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dmtx(:,:)
integer, intent(in) :: n_conf, memb(:,:), clsize(:)
integer, intent(out) :: ij_min(2)
real(dp), intent(out) :: d_min
integer :: i, j, ii, jj
real(dp) :: d

ij_min = 0
d_min = maximum_distance
do i = 1, n_conf-1
    if (clsize(i) == 0) cycle
    do j = i+1, n_conf
        if (clsize(j) == 0) cycle
        !
        d = 0.0d0
        do ii = 1, clsize(i)
            do jj = 1, clsize(j)
                d = max(dmtx(memb(ii,i),memb(jj,j)), d)
            end do
        end do
        !
        if (d < d_min) then
            ij_min = (/i,j/)
            d_min = d
        end if
    end do
end do

end subroutine linkage_maximum
!-------------------------------------------------------------------------------
subroutine linkage_minimum(dmtx, n_conf, memb, clsize, ij_min, d_min)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dmtx(:,:)
integer, intent(in) :: n_conf, memb(:,:), clsize(:)
integer, intent(out) :: ij_min(2)
real(dp), intent(out) :: d_min
integer :: i, j, ii, jj
real(dp) :: d

ij_min = 0
d_min = maximum_distance
do i = 1, n_conf-1
    if (clsize(i) == 0) cycle
    do j = i+1, n_conf
        if (clsize(j) == 0) cycle
        !
        d = maximum_distance
        do ii = 1, clsize(i)
            do jj = 1, clsize(j)
                d = min(dmtx(memb(ii,i),memb(jj,j)), d)
            end do
        end do
        !
        if (d < d_min) then
            ij_min = (/i,j/)
            d_min = d
        end if
    end do
end do

end subroutine linkage_minimum
!-------------------------------------------------------------------------------
subroutine linkage_average(dmtx, n_conf, memb, clsize, ij_min, d_min)
!-------------------------------------------------------------------------------
! UPGMA
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dmtx(:,:)
integer, intent(in) :: n_conf, memb(:,:), clsize(:)
integer, intent(out) :: ij_min(2)
real(dp), intent(out) :: d_min
integer :: i, j, ii, jj
real(dp) :: d

ij_min = 0
d_min = maximum_distance
do i = 1, n_conf-1
    if (clsize(i) == 0) cycle
    do j = i+1, n_conf
        if (clsize(j) == 0) cycle
        !
        d = 0.0d0
        do ii = 1, clsize(i)
            do jj = 1, clsize(j)
                d = d + dmtx(memb(ii,i), memb(jj,j))
            end do
        end do
        d = d / (clsize(i)*clsize(j))
        !
        if (d < d_min) then
            ij_min = (/i,j/)
            d_min = d
        end if
    end do
end do

end subroutine linkage_average
!-------------------------------------------------------------------------------
subroutine linkage_centroid(dmtx, n_conf, memb, clsize, ij_min, d_min)
!-------------------------------------------------------------------------------
! UPGMC
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dmtx(:,:)
integer, intent(in) :: n_conf, memb(:,:), clsize(:)
integer, intent(out) :: ij_min(2)
real(dp), intent(out) :: d_min
integer :: i, j
real(dp) :: d

ij_min = 0
d_min = maximum_distance
do i = 1, n_conf-1
    if (clsize(i) == 0) cycle
    do j = i+1, n_conf
        if (clsize(j) == 0) cycle
        !
        d = dmtx(memb(1,i), memb(1,j))
        if (d < d_min) then
            ij_min = (/i,j/)
            d_min = d
        end if
    end do
end do

end subroutine linkage_centroid
!-------------------------------------------------------------------------------
subroutine update_ward_distance(dmtx, n_conf, u,s,t, memb, clsize)
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: dmtx(:,:)
integer, intent(in) :: n_conf, u,s,t
integer, intent(in) :: memb(:,:), clsize(:)
integer :: v

do v = 1, n_conf
    if (clsize(v) == 0) cycle
    if (v == s .or. v == t) cycle
    dmtx(v,u) = ward_distance(dmtx, s,t,v, memb, clsize)
    dmtx(u,v) = dmtx(v,u)
end do

end subroutine update_ward_distance
!-------------------------------------------------------------------------------
function ward_distance(dmtx, s,t,v, memb, clsize)
!-------------------------------------------------------------------------------
real(dp) :: ward_distance
real(dp), intent(in) :: dmtx(:,:)
integer, intent(in) :: s, t, v
integer, intent(in) :: memb(:,:), clsize(:)

ward_distance = 0.0d0
ward_distance = ward_distance + (clsize(s)+clsize(v)) * dmtx(s,v)**2
ward_distance = ward_distance + (clsize(t)+clsize(v)) * dmtx(t,v)**2
ward_distance = ward_distance + (clsize(s)+clsize(t)) * dmtx(s,t)**2
ward_distance = dsqrt(ward_distance / (clsize(s)+clsize(t)+clsize(v)))

end function ward_distance
!-------------------------------------------------------------------------------
subroutine get_centroid(dmtx, n_conf, memb)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: dmtx(:,:)
integer, intent(in) :: n_conf
integer, intent(inout) :: memb(:)
integer :: i, j, i_cntr, cntr
real(dp) :: dsum(n_conf)

do i = 1, n_conf
    do j = 1, n_conf
        dsum(i) = dsum(i) + dmtx(memb(j),i)
    end do
end do

i_cntr = minloc(dsum, dim=1)
if (i_cntr /= 1) then
    cntr = memb(i_cntr)
    memb(2:i_cntr) = memb(1:i_cntr-1)
    memb(1) = cntr
end if

end subroutine get_centroid
!-------------------------------------------------------------------------------
subroutine remove_redundant(n_conf0, dmtx0, n_conf, dmtx, i_conf)
!-------------------------------------------------------------------------------
integer, intent(in) :: n_conf0
real(dp), intent(in) :: dmtx0(:,:)
integer, intent(out) :: n_conf
real(dp), intent(out), allocatable :: dmtx(:,:)
integer, intent(out), allocatable :: i_conf(:)
integer, allocatable :: uniq(:)
real(dp), parameter :: redundant_cutoff = 0.001d0
integer :: i, j
logical :: is_uniq

allocate(uniq(n_conf0))
n_conf = 0

do i = 1, n_conf0
    is_uniq = .true.
    do j = 1, i-1
        if (dmtx0(j,i) < redundant_cutoff) then
            is_uniq = .false.
            exit
        end if
    end do
    !
    if (is_uniq) then
        n_conf = n_conf + 1
        uniq(n_conf) = i
    end if
end do

allocate(i_conf(n_conf))
i_conf(1:n_conf) = uniq(1:n_conf)
deallocate(uniq)

allocate(dmtx(n_conf,n_conf))
do i = 1, n_conf
    do j = 1, n_conf
        dmtx(j,i) = dmtx0(i_conf(j), i_conf(i))
    end do
end do

end subroutine remove_redundant
!-------------------------------------------------------------------------------
END MODULE CLUSTER
!-------------------------------------------------------------------------------
