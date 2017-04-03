!-------------------------------------------------------------------------------
MODULE MathFunction
!-------------------------------------------------------------------------------
use globals

implicit none
public

CONTAINS
!-------------------------------------------------------------------------------
function v_size(v)
!-------------------------------------------------------------------------------
real(dp) :: v_size
real(dp), intent(in) :: v(:)

v_size = dsqrt(dot_product(v,v))

end function v_size
!-------------------------------------------------------------------------------
function v_norm(v)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: v(:)
real(dp), dimension(size(v)) :: v_norm

v_norm = v / v_size(v)

end function v_norm
!-------------------------------------------------------------------------------
function cross(u,v)
!-------------------------------------------------------------------------------
real(dp) :: cross(3)
real(dp), intent(in) :: u(3), v(3)

cross(1) = u(2)*v(3) - u(3)*v(2)
cross(2) = u(3)*v(1) - u(1)*v(3)
cross(3) = u(1)*v(2) - u(2)*v(1)

end function cross
!-------------------------------------------------------------------------------
function quaternion(axis, angle)
!-------------------------------------------------------------------------------
real(dp) :: quaternion(4)
real(dp), intent(in) :: axis(3)
real(dp), intent(in) :: angle

quaternion(1) = dcos(angle/2.0d0)
quaternion(2:4) = v_norm(axis(:)) * dsin(angle/2.0d0)

end function quaternion
!-------------------------------------------------------------------------------
function q_product(q1, q2)
!-------------------------------------------------------------------------------
real(dp) :: q_product(4)
real(dp), intent(in) :: q1(4), q2(4)
real(dp) :: tmp(3)

q_product(1) = q1(1)*q2(1) - dot_product(q1(2:4), q2(2:4))
tmp = cross(q1(2:4), q2(2:4))
q_product(2:4) = q1(1)*q2(2:4) + q2(1)*q1(2:4) + tmp

end function q_product
!-------------------------------------------------------------------------------
function rotation_matrix(q)
!-------------------------------------------------------------------------------
real(dp) :: rotation_matrix(3,3)
real(dp), intent(in) :: q(4)

rotation_matrix(1,1) = q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2
rotation_matrix(1,2) = 2.0d0*(q(2)*q(3) - q(1)*q(4))
rotation_matrix(1,3) = 2.0d0*(q(2)*q(4) + q(1)*q(3))

rotation_matrix(2,1) = 2.0d0*(q(2)*q(3) + q(1)*q(4))
rotation_matrix(2,2) = q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2
rotation_matrix(2,3) = 2.0d0*(q(3)*q(4) - q(1)*q(2))

rotation_matrix(3,1) = 2.0d0*(q(2)*q(4) - q(1)*q(3))
rotation_matrix(3,2) = 2.0d0*(q(3)*q(4) + q(1)*q(2))
rotation_matrix(3,3) = q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2

end function rotation_matrix
!-------------------------------------------------------------------------------
function rotation_matrix2(axis, angle)
!-------------------------------------------------------------------------------
real(dp) :: rotation_matrix2(3,3)
real(dp), intent(in) :: axis(3)
real(dp), intent(in) :: angle

rotation_matrix2 = rotation_matrix(quaternion(axis, angle))

end function rotation_matrix2
!-------------------------------------------------------------------------------
function rotate(v, U)
!-------------------------------------------------------------------------------
real(dp) :: rotate(3)
real(dp), intent(in) :: v(3), U(3,3)

rotate(1) = dot_product(U(1,:), v)
rotate(2) = dot_product(U(2,:), v)
rotate(3) = dot_product(U(3,:), v)

end function rotate
!-------------------------------------------------------------------------------
function rotate2(v, axis, angle)
!-------------------------------------------------------------------------------
real(dp) :: rotate2(3)
real(dp), intent(in) :: v(3), axis(3), angle

rotate2 = rotate(v, rotation_matrix2(axis, angle))

end function rotate2
!-------------------------------------------------------------------------------
function bound_angle(ang)
!-------------------------------------------------------------------------------
real(dp) :: bound_angle
real(dp), intent(in) :: ang

if (ang < -pi) then
    bound_angle = ang + 2.0d0*pi
else if (ang > pi) then
    bound_angle = ang - 2.0d0*pi
else
    bound_angle = ang
end if

end function bound_angle
!-------------------------------------------------------------------------------
function bound_angle90(ang)
!-------------------------------------------------------------------------------
real(dp) :: bound_angle90
real(dp), intent(in) :: ang
real(dp), parameter :: bound = 0.5d0*pi

if (ang < -bound) then
    bound_angle90 = pi + ang
else if (ang > bound) then
    bound_angle90 = pi - ang
else
    bound_angle90 = ang
end if

end function bound_angle90
!-------------------------------------------------------------------------------
function atan2(y, x)
!-------------------------------------------------------------------------------
! Arctangent for x and its perpendicular y
!-------------------------------------------------------------------------------
real(dp) :: atan2
real(dp), intent(in) :: x,y
real(dp) :: phi

phi = atan(abs(y/x))

if (x > 0) then
    atan2 = phi
else if (x < 0) then
    atan2 = pi - phi
else
    atan2 = 0.5d0 * pi
end if

atan2 = sign(atan2, y)

end function atan2
!-------------------------------------------------------------------------------
recursive subroutine combinations(n, r, n_combinations, comb)
!-------------------------------------------------------------------------------
integer, intent(in) :: n, r
integer, intent(out) :: n_combinations
integer, intent(out), allocatable :: comb(:,:)
integer :: i, j, k, n_sub
integer, allocatable :: sub(:,:)

if (r == 0) return

! n_combinations = n!/r!/(n-r)!
n_combinations = 1
do i = 0, r-1
    n_combinations = n_combinations * (n-i)
end do
do i = 1, r
    n_combinations = n_combinations / i
end do

allocate(comb(r,n_combinations))

if (r > 1) then
    call combinations(n, r-1, n_sub, sub)
    k = 0
    do i = 1, n
        do j = 1, n_sub
            if (i >= sub(1,j)) cycle
            !
            k = k + 1
            comb(1,k) = i 
            comb(2:r,k) = sub(1:r-1,j)
        end do
    end do
    deallocate(sub)
else
    do i = 1, n
        comb(1,i) = i
    end do
end if

end subroutine combinations
!-------------------------------------------------------------------------------
END MODULE MathFunction
!-------------------------------------------------------------------------------
