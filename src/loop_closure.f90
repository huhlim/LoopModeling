!-------------------------------------------------------------------------------
MODULE LOOP_CLOSURE
!-------------------------------------------------------------------------------
use globals
use logger, only: terminate_with_error
use geometry, only: calc_angle, calc_torsion, rotation_matrix2
use mathfunction, only: atan2, cross

implicit none
private

integer, parameter :: max_soln = 16
integer, parameter :: deg_poly = 16

!-------------------------------------------------------------------------------
! DERIVED TYPES
!-------------------------------------------------------------------------------
type tripep_geometry_type
!-------------------------------------------------------------------------------
real(dp) :: b_len0(6), b_ang0(7), t_ang0(2) ! Input geometry
real(dp) :: r32(3), r12(3), r43(3)          ! vectors between input coordinates
real(dp) :: u32(3), u12(3), u43(3)          ! unit vectors of rXX
real(dp) :: len_na(3), len_ac(3), len_aa(3)

end type tripep_geometry_type
!-------------------------------------------------------------------------------
type tripep_angle_type
!-------------------------------------------------------------------------------
! Defining tripep_closure related angles
!-------------------------------------------------------------------------------
real(dp) :: angle
real(dp) :: sin_angle
real(dp) :: cos_angle
!-------------------------------------------------------------------------------
end type tripep_angle_type
!-------------------------------------------------------------------------------
type(tripep_geometry_type) :: geom
!
type(tripep_angle_type) :: delta_lc(0:3)
type(tripep_angle_type) :: xi_lc(3)
type(tripep_angle_type) :: eta_lc(3)
type(tripep_angle_type) :: alpha_lc(3)
type(tripep_angle_type) :: theta_lc(3)

!-------------------------------------------------------------------------------
! parameters for tripeptide loop (including bond lengths & angles)
real(dp) :: dsq_r32_min, dsq_r32_max

! used for polynomial coefficients
real(dp) :: C0(0:2,3), C1(0:2,3), C2(0:2,3)
real(dp) :: Q(0:16,0:4), RM(0:16,0:2)

real(dp) :: r_tri_rot0(3,3,3,max_soln), r_tri_rot(3,3,3,43,max_soln)
logical :: present_tmp = .false.

!-------------------------------------------------------------------------------
public :: solve_tripep_closure
public :: max_soln

CONTAINS
!-------------------------------------------------------------------------------
subroutine solve_tripep_closure(b_len, b_ang, t_ang, r_anchor, n_soln, soln)
!-------------------------------------------------------------------------------
! Solve the tripeptide closure problem
!  Input: r_anchor(3,4) -> coord. for 1st vertex N/CA, 3rd vertex CA/C
!  Output: n_soln -> number of solutions
!          soln(3,3,3,max_soln) -> solutions; xyz, N/CA/C, 1/2/3
!-------------------------------------------------------------------------------
real(dp), intent(in) :: b_len(6), b_ang(7), t_ang(2)
real(dp), intent(in)  :: r_anchor(3,4)
!
integer,  intent(out) :: n_soln
real(dp), intent(out) :: soln(3,3,3,max_soln)
!
real(dp) :: poly_coeff(0:deg_poly), roots(max_soln)

call initialize_tripep_closure(b_len, b_ang, t_ang)

call get_input_angles(r_anchor, n_soln)
if (n_soln == 0) return

call get_poly_coeff(poly_coeff)

call solve_sturm(deg_poly, n_soln, poly_coeff, roots)
if (n_soln == 0) return

call coord_from_poly_roots(n_soln, roots, r_anchor, soln)

end subroutine solve_tripep_closure
!-------------------------------------------------------------------------------
subroutine initialize_tripep_closure(b_len, b_ang, t_ang)
!-------------------------------------------------------------------------------
! Input angles for the given bond lengths and angles 
!-------------------------------------------------------------------------------
real(dp), intent(in) :: b_len(6), b_ang(7), t_ang(2)
!
real(dp) :: len1, len2, a_min, a_max
real(dp), dimension(3) :: axis, rr_a1, rr_c1, rr_n2, rr_a2, rr_n2a2_ref, rr_c1a1
real(dp), dimension(3) :: rr_a1a2, dr, bb_c1a1, bb_a1a2, bb_a2n2
real(dp) :: Us(3,3)
real(dp), parameter :: tol_secant = 1.0d-15
integer, parameter :: max_iter_sturm = 100, max_iter_secant = 20
integer :: i
  
call initialize_sturm(tol_secant, max_iter_sturm, max_iter_secant)

geom%b_len0(1:6) = b_len(1:6)
geom%b_ang0(1:7) = b_ang(1:7)
geom%t_ang0(1:2) = t_ang(1:2)

rr_c1(1:3) = 0.0d0
axis(1:3) = (/ 1.0d0, 0.0d0, 0.0d0 /)

do i = 0, 1
    rr_a1(1:3) = (/ cos(geom%b_ang0(3*i+2))*geom%b_len0(3*i+1), &
                    sin(geom%b_ang0(3*i+2))*geom%b_len0(3*i+1), 0.0d0 /)
    rr_n2(1:3) = (/ geom%b_len0(3*i+2), 0.0d0, 0.0d0 /)
    rr_c1a1(:) = rr_a1(:) - rr_c1(:)
    rr_n2a2_ref(1:3) = (/ -cos(geom%b_ang0(3*i+3))*geom%b_len0(3*i+3), &
                           sin(geom%b_ang0(3*i+3))*geom%b_len0(3*i+3), 0.0d0 /)
    Us = rotation_matrix2(axis, geom%t_ang0(i+1))
    rr_a2(:) =  matmul(Us, rr_n2a2_ref) + rr_n2(:)
    rr_a1a2(:) = rr_a2(:) - rr_a1(:)
    dr(:) = rr_a1a2(:)
    len2 = dot_product(dr, dr)
    len1 = sqrt(len2)
    ! len_aa
    geom%len_aa(i+2) = len1
    bb_c1a1(:) = rr_c1a1(:)/geom%b_len0(3*i+1)
    bb_a1a2(:) = rr_a1a2(:)/len1
    bb_a2n2(:) = (rr_n2(:) - rr_a2(:))/geom%b_len0(3*i+3)
    ! xi
    xi_lc(i+2)%angle = calc_angle(bb_a1a2, bb_a2n2)

    ! eta
    eta_lc(i+1)%angle = calc_angle(-bb_a1a2, -bb_c1a1)

    ! delta: pi -  dih of N(1)CA(1)CA(3)C(3)
    delta_lc(i+1)%angle = pi - calc_torsion(bb_c1a1, bb_a1a2, bb_a2n2)
end do

a_min = b_ang(4) - (xi_lc(2)%angle + eta_lc(2)%angle)
a_max = min(b_ang(4) + (xi_lc(2)%angle + eta_lc(2)%angle), pi)

! min/max of base length
dsq_r32_min = geom%len_aa(2)**2 + geom%len_aa(3)**2 - 2.0d0*geom%len_aa(2)*geom%len_aa(3)*cos(a_min)
dsq_r32_max = geom%len_aa(2)**2 + geom%len_aa(3)**2 - 2.0d0*geom%len_aa(2)*geom%len_aa(3)*cos(a_max)

end subroutine initialize_tripep_closure
!-------------------------------------------------------------------------------
subroutine get_input_angles(r_anchor, n_soln)
!-------------------------------------------------------------------------------
real(dp), intent(in)  :: r_anchor(3,4)
integer,  intent(out) :: n_soln
!
real(dp) :: dsq_r32
integer :: i

n_soln = max_soln

! Virtual bond
geom%r32(:) = r_anchor(:,3) - r_anchor(:,2)
dsq_r32 = dot_product(geom%r32, geom%r32)
geom%len_aa(1) = sqrt(dsq_r32)

if (dsq_r32 < dsq_r32_min .or. dsq_r32 > dsq_r32_max) then
    n_soln = 0
    return
end if

! Bond lengths
geom%r12(:)    = r_anchor(:,1) - r_anchor(:,2)
geom%len_na(1) = sqrt(dot_product(geom%r12, geom%r12))
geom%len_na(2) = geom%b_len0(3)
geom%len_na(3) = geom%b_len0(6)
!
geom%r43(:)    = r_anchor(:,4) - r_anchor(:,3)
geom%len_ac(1) = geom%b_len0(1)
geom%len_ac(2) = geom%b_len0(4)
geom%len_ac(3) = sqrt(dot_product(geom%r43, geom%r43))

! Unit vectors
geom%u12(:) = geom%r12(:)/geom%len_na(1)
geom%u43(:) = geom%r43(:)/geom%len_ac(3)
geom%u32(:) = geom%r32(:)/geom%len_aa(1)

! delta_lc(3): dih of N(1)CA(1)CA(3)C(3)
delta_lc(3)%angle = calc_torsion(-geom%u12, geom%u32, geom%u43)
delta_lc(0)%angle = delta_lc(3)%angle

! xi_lc(1)
xi_lc(1)%angle = calc_angle(geom%u32, geom%u12)
 
! eta_lc(3)
eta_lc(3)%angle = calc_angle(-geom%u32, geom%u43)

do i = 1, 3
    delta_lc(i)%cos_angle = cos(delta_lc(i)%angle)
    delta_lc(i)%sin_angle = sin(delta_lc(i)%angle)
    xi_lc(i)%cos_angle = cos(xi_lc(i)%angle)
    xi_lc(i)%sin_angle = sin(xi_lc(i)%angle)
    eta_lc(i)%cos_angle = cos(eta_lc(i)%angle)
    eta_lc(i)%sin_angle = sin(eta_lc(i)%angle)
end do
delta_lc(0)%cos_angle = delta_lc(3)%cos_angle
delta_lc(0)%sin_angle = delta_lc(3)%sin_angle

! theta (N, CA, C) bond angle
theta_lc(1)%angle = geom%b_ang0(1)
theta_lc(2)%angle = geom%b_ang0(4)
theta_lc(3)%angle = geom%b_ang0(7)
do i = 1, 3
    theta_lc(i)%cos_angle = cos(theta_lc(i)%angle)
end do

! alpha_lc
alpha_lc(1)%cos_angle = -(geom%len_aa(1)**2 + geom%len_aa(2)**2 - geom%len_aa(3)**2)/(2.0d0*geom%len_aa(1)*geom%len_aa(2))
alpha_lc(1)%angle     = acos(alpha_lc(1)%cos_angle)
alpha_lc(1)%sin_angle =  sin(alpha_lc(1)%angle)

alpha_lc(2)%cos_angle = (geom%len_aa(2)**2 + geom%len_aa(3)**2 - geom%len_aa(1)**2)/(2.0d0*geom%len_aa(2)*geom%len_aa(3))
alpha_lc(2)%angle     = acos(alpha_lc(2)%cos_angle)
alpha_lc(2)%sin_angle =  sin(alpha_lc(2)%angle)

alpha_lc(3)%angle     = pi - alpha_lc(1)%angle + alpha_lc(2)%angle
alpha_lc(3)%cos_angle = cos(alpha_lc(3)%angle)
alpha_lc(3)%sin_angle = sin(alpha_lc(3)%angle)

! check for existence of soln
do i = 1, 3
    call test_two_cone_existence_soln(theta_lc(i)%angle, xi_lc(i)%angle, eta_lc(i)%angle, alpha_lc(i)%angle, n_soln)
    if (n_soln == 0) return
end do

end subroutine get_input_angles
!-------------------------------------------------------------------------------
subroutine test_two_cone_existence_soln(tt, kx, et, ap, n_soln)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: tt, kx, et, ap
integer, intent(out) :: n_soln
real(dp) :: at, ex, abs_at, ap1, kx1, et1
real(dp) :: cos_tx1, cos_tx2, cos_te1, cos_te2, cos_ea1, cos_ea2, cos_xa1, cos_xa2
logical :: s1, s2, t1, t2, complicated = .false.

n_soln = max_soln
 
ap1 = ap
kx1 = kx
et1 = et
  
at = ap1 - tt
ex = kx1 + et1
abs_at = abs(at)

! case of no soln
if (abs_at > ex) then
    n_soln = 0
    return
end if

if (complicated) then
    ! find type of intersection
    cos_tx1 = cos(tt+kx1)
    cos_tx2 = cos(tt-kx1)
    cos_te1 = cos(tt+et1)
    cos_te2 = cos(tt-et1)
    cos_ea1 = cos(et1+ap1)
    cos_ea2 = cos(et1-ap1)
    cos_xa1 = cos(kx1+ap1)
    cos_xa2 = cos(kx1-ap1)
    s1 = .false.; s2 = .false.; t1 = .false.; t2 = .false.
    if ((cos_te1-cos_xa2)*(cos_te1-cos_xa1) <= 0.0d0) s1 = .true.
    if ((cos_te2-cos_xa2)*(cos_te2-cos_xa1) <= 0.0d0) s2 = .true.
    if ((cos_tx1-cos_ea2)*(cos_tx1-cos_ea1) <= 0.0d0) t1 = .true.
    if ((cos_tx2-cos_ea2)*(cos_tx2-cos_ea1) <= 0.0d0) t2 = .true.
     
end if

end subroutine test_two_cone_existence_soln
!-------------------------------------------------------------------------------
subroutine get_poly_coeff(poly_coeff)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: poly_coeff(0:deg_poly)
integer :: i
real(dp) :: Ax(0:4), Axx(2,3)
real(dp) :: Bx(3,0:8)
real(dp), dimension(0:4,0:4) :: u11, u12, u13, u31, u32, u33
real(dp), dimension(0:4,0:4) :: q_tmp
real(dp), dimension(0:4,0:4, 6) :: um
integer :: p1(2), p2, p3(2), p4, p_um(2,6), p_final, p_fx(27), p_Q(2)
real(dp), dimension(0:16, 27) :: fx

! A0, B0
do i = 1, 3
    Ax(0) =  alpha_lc(i)%cos_angle*xi_lc(i)%cos_angle*eta_lc(i)%cos_angle - theta_lc(i)%cos_angle
    Ax(1) = -alpha_lc(i)%sin_angle*xi_lc(i)%cos_angle*eta_lc(i)%sin_angle
    Ax(2) =  alpha_lc(i)%sin_angle*xi_lc(i)%sin_angle*eta_lc(i)%cos_angle
    Ax(3) =                        xi_lc(i)%sin_angle*eta_lc(i)%sin_angle
    Ax(4) =  alpha_lc(i)%cos_angle*Ax(3)
    !
    Axx(1,1:3) = Ax(2:4)*delta_lc(i-1)%cos_angle
    Axx(2,1:3) = Ax(2:4)*delta_lc(i-1)%sin_angle
    !
    Bx(i,0) = Ax(0) + Axx(2,1) + Axx(1,2)
    Bx(i,1) = 2.0d0 *(Ax(1)    + Axx(2,3))
    Bx(i,2) = 2.0d0 *(Axx(2,2) - Axx(1,1))
    Bx(i,3) =-4.0d0 * Axx(1,3)
    Bx(i,4) = Ax(0) + Axx(2,1) - Axx(1,2)
    Bx(i,5) = Ax(0) - Axx(2,1) - Axx(1,2)
    Bx(i,6) =-2.0d0 *(Axx(1,1) + Axx(2,2))
    Bx(i,7) = 2.0d0 *(Ax(1)    - Axx(2,3))
    Bx(i,8) = Ax(0) - Axx(2,1) + Axx(1,2)
end do

! C0i
C0(0:2,1) = (/ Bx(1,0), Bx(1,2), Bx(1,5) /)
C1(0:2,1) = (/ Bx(1,1), Bx(1,3), Bx(1,7) /)
C2(0:2,1) = (/ Bx(1,4), Bx(1,6), Bx(1,8) /)
do i = 2, 3
    C0(0:2,i) = (/ Bx(i,0), Bx(i,1), Bx(i,4) /)
    C1(0:2,i) = (/ Bx(i,2), Bx(i,3), Bx(i,6) /)
    C2(0:2,i) = (/ Bx(i,5), Bx(i,7), Bx(i,8) /)
end do

! first determinant
do i = 0, 2
    u11(i,0) = C0(i,1)
    u12(i,0) = C1(i,1)
    u13(i,0) = C2(i,1)
    u31(0,i) = C0(i,2)
    u32(0,i) = C1(i,2)
    u33(0,i) = C2(i,2)
end do

p1(1:2) = (/ 2, 0 /)
p3(1:2) = (/ 0, 2 /)

call poly_mul_sub2(u32, u32, u31, u33, p3, p3, p3, p3, um(:,:,1), p_um(:,1))
call poly_mul_sub2(u12, u32, u11, u33, p1, p3, p1, p3, um(:,:,2), p_um(:,2))
call poly_mul_sub2(u12, u33, u13, u32, p1, p3, p1, p3, um(:,:,3), p_um(:,3))
call poly_mul_sub2(u11, u33, u31, u13, p1, p3, p3, p1, um(:,:,4), p_um(:,4))
!
call poly_mul_sub2(u13, um(:,:,1), u33, um(:,:,2), p1, p_um(:,1), p3, p_um(:,2), um(:,:,5), p_um(:,5))
call poly_mul_sub2(u13, um(:,:,4), u12, um(:,:,3), p1, p_um(:,4), p1, p_um(:,3), um(:,:,6), p_um(:,6))
!
call poly_mul_sub2(u11, um(:,:,5), u31, um(:,:,6), p1, p_um(:,5), p3, p_um(:,6), q_tmp, p_Q)
!
Q(0:4,0:4) = q_tmp(0:4,0:4)

! second determinant
RM(:,:) = 0.0d0
RM(0:2,0) = C0(0:2,3)
RM(0:2,1) = C1(0:2,3)
RM(0:2,2) = C2(0:2,3)
p2 = 2
p4 = 4

call poly_mul_sub1(RM(:,1), RM(:,1), RM(:,0), RM(:,2), p2, p2,      p2, p2,      fx(:,1), p_fx(1))
call poly_mul1    (RM(:,1), RM(:,2), p2, p2,      fx(:,2), p_fx(2))
call poly_mul_sub1(RM(:,1), fx(:,1), RM(:,0), fx(:,2), p2, p_fx(1), p2, p_fx(2), fx(:,3), p_fx(3))
call poly_mul1    (RM(:,2), fx(:,1), p2, p_fx(1), fx(:,4), p_fx(4))
call poly_mul_sub1(RM(:,1), fx(:,3), RM(:,0), fx(:,4), p2, p_fx(3), p2, p_fx(4), fx(:,5), p_fx(5))
!
call poly_mul_sub1(Q(:,1), RM(:,1), Q(:,0),  RM(:,2), p4, p2,      p4, p2,      fx(:,6), p_fx(6))
call poly_mul_sub1(Q(:,2), fx(:,1), RM(:,2), fx(:,6), p4, p_fx(1), p2, p_fx(6), fx(:,7), p_fx(7))
call poly_mul_sub1(Q(:,3), fx(:,3), RM(:,2), fx(:,7), p4, p_fx(3), p2, p_fx(7), fx(:,8), p_fx(8))
call poly_mul_sub1(Q(:,4), fx(:,5), RM(:,2), fx(:,8), p4, p_fx(5), p2, p_fx(8), fx(:,9), p_fx(9))
!
call poly_mul_sub1(Q(:,3), RM(:,1), Q(:,4),  RM(:,0),  p4, p2,      p4, p2,       fx(:,10), p_fx(10))
call poly_mul_sub1(Q(:,2), fx(:,1), RM(:,0), fx(:,10), p4, p_fx(1), p2, p_fx(10), fx(:,11), p_fx(11))
call poly_mul_sub1(Q(:,1), fx(:,3), RM(:,0), fx(:,11), p4, p_fx(3), p2, p_fx(11), fx(:,12), p_fx(12))
!
call poly_mul_sub1(Q(:,2), RM(:,1),  Q(:,1),  RM(:,2),  p4, p2,       p4, p2,       fx(:,13), p_fx(13))
call poly_mul_sub1(Q(:,3), fx(:,1),  RM(:,2), fx(:,13), p4, p_fx(1),  p2, p_fx(13), fx(:,14), p_fx(14))
call poly_mul_sub1(Q(:,3), RM(:,1),  Q(:,2),  RM(:,2),  p4, p2,       p4, p2,       fx(:,15), p_fx(15))
call poly_mul_sub1(Q(:,4), fx(:,1),  RM(:,2), fx(:,15), p4, p_fx(1),  p2, p_fx(15), fx(:,16), p_fx(16))
call poly_mul_sub1(Q(:,1), fx(:,14), Q(:,0),  fx(:,16), p4, p_fx(14), p4, p_fx(16), fx(:,17), p_fx(17))
!
call poly_mul_sub1(Q(:,2), RM(:,2),  Q(:,3), RM(:,1),  p4, p2,       p4, p2,       fx(:,18), p_fx(18))
call poly_mul_sub1(Q(:,1), RM(:,2),  Q(:,3), RM(:,0),  p4, p2,       p4, p2,       fx(:,19), p_fx(19))
call poly_mul_sub1(Q(:,3), fx(:,19), Q(:,2), fx(:,18), p4, p_fx(19), p4, p_fx(18), fx(:,20), p_fx(20))
call poly_mul_sub1(Q(:,1), RM(:,1),  Q(:,2), RM(:,0),  p4, p2,       p4, p2,       fx(:,21), p_fx(21))
!
call poly_mul1(Q(:,4),   fx(:,21), p4,       p_fx(21), fx(:,22), p_fx(22))
call poly_sub1(fx(:,20), fx(:,22), p_fx(20), p_fx(22), fx(:,23), p_fx(23))
call poly_mul1(RM(:,0),  fx(:,23), p2,       p_fx(23), fx(:,24), p_fx(24))
call poly_sub1(fx(:,17), fx(:,24), p_fx(17), p_fx(24), fx(:,25), p_fx(25))
!
call poly_mul_sub1(Q(:,4), fx(:,12), RM(:,2), fx(:,25), p4, p_fx(12), p2, p_fx(25), fx(:,26),   p_fx(26))
call poly_mul_sub1(Q(:,0), fx(:,9),  RM(:,0), fx(:,26), p4, p_fx(9),  p2, p_fx(26), poly_coeff, p_final)

if (p_final /= deg_poly) then
    call terminate_with_error('Error: degree of polynomial is not 16!')
end if

if (poly_coeff(16) < 0.0d0) then
    poly_coeff(0:16) = -poly_coeff(0:16)
end if

end subroutine get_poly_coeff
!-------------------------------------------------------------------------------
subroutine coord_from_poly_roots(n_soln, roots, r_anchor, soln)
!-------------------------------------------------------------------------------
integer, intent(in) :: n_soln
real(dp), intent(in) :: roots(n_soln), r_anchor(3,4)
real(dp), intent(out) :: soln(3,3,3,max_soln)

real(dp) :: ex(3), ey(3), ez(3), b_a1a2(3), b_a3a2(3), r_tmp(3)
real(dp) :: p_s(3,3), s1(3,3), s2(3,3), p_t(3,3), t1(3,3), t2(3,3)
real(dp) :: p_s_c(3,3), s1_s(3,3), s2_s(3,3), p_t_c(3,3), t1_s(3,3), t2_s(3,3)
real(dp) :: angle, sig1_init, half_tan(3)
real(dp) :: cos_tau(0:3), sin_tau(0:3), cos_sig(3), sin_sig(3), ht, tmp, sig1
real(dp) :: r_s(3), r_t(3), r0(3), r_n(3,3), r_a(3,3), r_c(3,3), Us(3,3)
integer :: i_soln, i, j
real(dp) :: tau1, tau2, tau3, tau1_curr, tau2_curr, tau3_curr
integer :: i_frame, ii

if (n_soln == 0) return

! Define body frame (ex, ey, ez)
ex(:) = geom%u32(:)
ez = cross(geom%r12, ex)
ez(:) = ez(:)/sqrt(dot_product(ez,ez))
ey = cross(ez, ex)
! virtual bond vectors in the reference plane
b_a1a2(:) = -alpha_lc(1)%cos_angle*ex(:) + alpha_lc(1)%sin_angle*ey(:)
b_a3a2(:) =  alpha_lc(3)%cos_angle*ex(:) + alpha_lc(3)%sin_angle*ey(:)
!! Define cone coordinates for each angle joint.
! (p_s,s1,s2) and (p_t,t1,t2):  Right Orthonormal systems
! residue 1
p_s(:,1) = -ex(:)
s1(:,1)  = ez(:)  ! (p_s)X(p_t)/||(p_s)X(p_t)||
s2(:,1)  = ey(:)  ! p_s X s1
p_t(:,1) = b_a1a2(:)
t1(:,1)  = ez(:)  ! s1
t2(:,1)  = alpha_lc(1)%sin_angle*ex(:) + alpha_lc(1)%cos_angle*ey(:) ! p_t X t1
! residue 2
p_s(:,2) = -b_a1a2(:)
s1(:,2)  = -ez(:)
s2(:,2)  = t2(:,1)  ! sina1*ex(:) + cosa1*ey(:)
p_t(:,2) = -b_a3a2(:)
t1(:,2)  = -ez(:)
t2(:,2)  = alpha_lc(3)%sin_angle*ex(:) - alpha_lc(3)%cos_angle*ey(:)
! residue 3
p_s(:,3) = b_a3a2(:)
s2(:,3)  = t2(:,2)   ! sina3*ex(:) + cosa3*ey(:)
s1(:,3)  = ez(:)
p_t(:,3) = ex(:)
t1(:,3) =  ez(:)
t2(:,3) = -ey(:)
! scale vectors
do i = 1, 3
    p_s_c(:,i) = p_s(:,i)*xi_lc(i)%cos_angle
    s1_s(:,i)  =  s1(:,i)*xi_lc(i)%sin_angle
    s2_s(:,i)  =  s2(:,i)*xi_lc(i)%sin_angle
    p_t_c(:,i) = p_t(:,i)*eta_lc(i)%cos_angle
    t1_s(:,i)  =  t1(:,i)*eta_lc(i)%sin_angle
    t2_s(:,i)  =  t2(:,i)*eta_lc(i)%sin_angle
end do

! initial sig(1)
r_tmp(:) = (geom%r12(:)/geom%len_na(1) - p_s_c(:,1))/xi_lc(1)%sin_angle
angle = calc_angle(-s1(:,1), r_tmp)
sig1_init = sign(angle, dot_product(r_tmp(:),s2(:,1)))

! CA
r_a(:,1) = r_anchor(:,2)
r_a(:,2) = r_anchor(:,2) + geom%len_aa(2)*b_a1a2(:)
r_a(:,3) = r_anchor(:,3)
r0(:) = r_anchor(:,2)
  
if (present_tmp) then
    do i_soln = 1, n_soln
        do i = 1, 3
            cos_tau(i) = 0.0d0
            sin_tau(i) = 1.0d0
        end do
        do i = 1, 3
            j = i - 1
            cos_sig(i) = delta_lc(j)%cos_angle*cos_tau(j) + delta_lc(j)%sin_angle*sin_tau(j)
            sin_sig(i) = delta_lc(j)%sin_angle*cos_tau(j) - delta_lc(j)%cos_angle*sin_tau(j)
        end do
        do i = 1, 3
            r_s(:) = p_s_c(:,i) + cos_sig(i)*s1_s(:,i) + sin_sig(i)*s2_s(:,i)
            r_t(:) = p_t_c(:,i) + cos_tau(i)*t1_s(:,i) + sin_tau(i)*t2_s(:,i)
            r_n(:,i) = r_s(:)*geom%len_na(i) + r_a(:,i)
            r_c(:,i) = r_t(:)*geom%len_ac(i) + r_a(:,i)
            r_tri_rot0(:,1,i,i_soln) = r_n(:,i)
            r_tri_rot0(:,2,i,i_soln) = r_a(:,i)
            r_tri_rot0(:,3,i,i_soln) = r_c(:,i)
        end do
    end do
     
    i_frame = 0
    do i_soln = 1, n_soln
        do i = 1, 3
            cos_tau(i) = 0.0d0
            sin_tau(i) = 1.0d0
        end do
        half_tan(3) = roots(i_soln)
        half_tan(2) = calc_t2(half_tan(3))
        half_tan(1) = calc_t1(half_tan(3), half_tan(2))
        do i = 1, 1
            ht = half_tan(i)
            tmp = 1.0d0 + ht*ht
            cos_tau(i) = (1.0d0 - ht*ht)/tmp
            sin_tau(i) = 2.0d0*ht/tmp
        end do
        tau1_curr = atan2(sin_tau(1), cos_tau(1))
        cos_tau(0) = cos_tau(3)
        sin_tau(0) = sin_tau(3)

        do ii = 0, 10
            tau1 = pi*0.5d0 + (tau1_curr - pi*0.5d0)*dble(ii)/10.0d0
            cos_tau(1) = cos(tau1)
            sin_tau(1) = sin(tau1)
            do i = 1, 3
                j = i - 1
                cos_sig(i) = delta_lc(j)%cos_angle*cos_tau(j) + delta_lc(j)%sin_angle*sin_tau(j)
                sin_sig(i) = delta_lc(j)%sin_angle*cos_tau(j) - delta_lc(j)%cos_angle*sin_tau(j)
            end do
            i_frame = i_frame + 1
            do i = 1, 3
                r_s(:) = p_s_c(:,i) + cos_sig(i)*s1_s(:,i) + sin_sig(i)*s2_s(:,i)
                r_t(:) = p_t_c(:,i) + cos_tau(i)*t1_s(:,i) + sin_tau(i)*t2_s(:,i)
                r_n(:,i) = r_s(:)*geom%len_na(i) + r_a(:,i)
                r_c(:,i) = r_t(:)*geom%len_ac(i) + r_a(:,i)
                r_tri_rot(:,1,i,i_frame,i_soln) = r_n(:,i)
                r_tri_rot(:,2,i,i_frame,i_soln) = r_a(:,i)
                r_tri_rot(:,3,i,i_frame,i_soln) = r_c(:,i)
            end do
        end do

        do ii = 1, 5
            i_frame = i_frame + 1
            do i = 1, 3
                r_tri_rot(:,1,i,i_frame,i_soln) = r_n(:,i)
                r_tri_rot(:,2,i,i_frame,i_soln) = r_a(:,i)
                r_tri_rot(:,3,i,i_frame,i_soln) = r_c(:,i)
            end do
        end do

        do i = 1, 3
            cos_tau(i) = 0.0d0
            sin_tau(i) = 1.0d0
        end do
        half_tan(3) = roots(i_soln)
        half_tan(2) = calc_t2(half_tan(3))
        half_tan(1) = calc_t1(half_tan(3), half_tan(2))
        do i = 1, 2
            ht = half_tan(i)
            tmp = 1.0d0 + ht*ht
            cos_tau(i) = (1.0d0 - ht*ht)/tmp
            sin_tau(i) = 2.0d0*ht/tmp
        end do
        cos_tau(0) = cos_tau(3)
        sin_tau(0) = sin_tau(3)
        tau2_curr = atan2(sin_tau(2), cos_tau(2))

        do ii = 0, 10
            tau2 = pi*0.5d0 + (tau2_curr - pi*0.5d0)*dble(ii)/10.0d0
            cos_tau(2) = cos(tau2)
            sin_tau(2) = sin(tau2)

            do i = 1, 3
                j = i - 1
                cos_sig(i) = delta_lc(j)%cos_angle*cos_tau(j) + delta_lc(j)%sin_angle*sin_tau(j)
                sin_sig(i) = delta_lc(j)%sin_angle*cos_tau(j) - delta_lc(j)%cos_angle*sin_tau(j)
            end do
            i_frame = i_frame + 1
            do i = 1, 3
                r_s(:) = p_s_c(:,i) + cos_sig(i)*s1_s(:,i) + sin_sig(i)*s2_s(:,i)
                r_t(:) = p_t_c(:,i) + cos_tau(i)*t1_s(:,i) + sin_tau(i)*t2_s(:,i)
                r_n(:,i) = r_s(:)*geom%len_na(i) + r_a(:,i)
                r_c(:,i) = r_t(:)*geom%len_ac(i) + r_a(:,i)
                r_tri_rot(:,1,i,i_frame,i_soln) = r_n(:,i)
                r_tri_rot(:,2,i,i_frame,i_soln) = r_a(:,i)
                r_tri_rot(:,3,i,i_frame,i_soln) = r_c(:,i)
            end do
        end do

        do ii = 1, 10
            i_frame = i_frame + 1
            do i = 1, 3
                r_tri_rot(:,1,i,i_frame,i_soln) = r_n(:,i)
                r_tri_rot(:,2,i,i_frame,i_soln) = r_a(:,i)
                r_tri_rot(:,3,i,i_frame,i_soln) = r_c(:,i)
            end do
        end do

        half_tan(3) = roots(i_soln)
        half_tan(2) = calc_t2(half_tan(3))
        half_tan(1) = calc_t1(half_tan(3), half_tan(2))
        do i = 1, 3
            ht = half_tan(i)
            tmp = 1.0d0 + ht*ht
            cos_tau(i) = (1.0d0 - ht*ht)/tmp
            sin_tau(i) = 2.0d0*ht/tmp
        end do
        cos_tau(0) = cos_tau(3)
        sin_tau(0) = sin_tau(3)
        do i = 1, 3
            j = i - 1
            cos_sig(i) = delta_lc(j)%cos_angle*cos_tau(j) + delta_lc(j)%sin_angle*sin_tau(j)
            sin_sig(i) = delta_lc(j)%sin_angle*cos_tau(j) - delta_lc(j)%cos_angle*sin_tau(j)
        end do
        do i = 1, 3
            r_s(:) = p_s_c(:,i) + cos_sig(i)*s1_s(:,i) + sin_sig(i)*s2_s(:,i)
            r_t(:) = p_t_c(:,i) + cos_tau(i)*t1_s(:,i) + sin_tau(i)*t2_s(:,i)
            r_n(:,i) = r_s(:)*geom%len_na(i) + r_a(:,i)
            r_c(:,i) = r_t(:)*geom%len_ac(i) + r_a(:,i)
        end do
        
        ! rotate back atoms by -(sig(1) - sig1_init) around -ex
        sig1 = atan2(sin_sig(1), cos_sig(1))
        tau3_curr = -(sig1 - sig1_init)

        do ii = 0, 20
            tau3 = tau3_curr*dble(ii)/10.0d0
            Us = rotation_matrix2(-ex, tau3)
        
            soln(:,1,1,i_soln) = r_anchor(:,1)
            soln(:,2,1,i_soln) = r_anchor(:,2)
            soln(:,3,1,i_soln) = matmul(Us, r_c(:,1) - r0(:)) + r0(:)
            soln(:,1,2,i_soln) = matmul(Us, r_n(:,2) - r0(:)) + r0(:)
            soln(:,2,2,i_soln) = matmul(Us, r_a(:,2) - r0(:)) + r0(:)
            soln(:,3,2,i_soln) = matmul(Us, r_c(:,2) - r0(:)) + r0(:)
            soln(:,1,3,i_soln) = matmul(Us, r_n(:,3) - r0(:)) + r0(:)
            soln(:,2,3,i_soln) = r_anchor(:,3)
            soln(:,3,3,i_soln) = r_anchor(:,4)

            i_frame = i_frame + 1
            do i = 1, 3
                r_tri_rot(:,1,i,i_frame,i_soln) = soln(:,1,i,i_soln)
                r_tri_rot(:,2,i,i_frame,i_soln) = soln(:,2,i,i_soln)
                r_tri_rot(:,3,i,i_frame,i_soln) = soln(:,3,i,i_soln)
            end do
        end do

    end do

end if

do i_soln = 1, n_soln
    half_tan(3) = roots(i_soln)
    half_tan(2) = calc_t2(half_tan(3))
    half_tan(1) = calc_t1(half_tan(3), half_tan(2))
    do i = 1, 3
        ht = half_tan(i)
        tmp = 1.0d0 + ht*ht
        cos_tau(i) = (1.0d0 - ht*ht)/tmp
        sin_tau(i) = 2.0d0*ht/tmp
    end do
    cos_tau(0) = cos_tau(3)
    sin_tau(0) = sin_tau(3)
    do i = 1, 3
        j = i - 1
        cos_sig(i) = delta_lc(j)%cos_angle*cos_tau(j) + delta_lc(j)%sin_angle*sin_tau(j)
        sin_sig(i) = delta_lc(j)%sin_angle*cos_tau(j) - delta_lc(j)%cos_angle*sin_tau(j)
    end do
    do i = 1, 3
        r_s(:) = p_s_c(:,i) + cos_sig(i)*s1_s(:,i) + sin_sig(i)*s2_s(:,i)
        r_t(:) = p_t_c(:,i) + cos_tau(i)*t1_s(:,i) + sin_tau(i)*t2_s(:,i)
        r_n(:,i) = r_s(:)*geom%len_na(i) + r_a(:,i)
        r_c(:,i) = r_t(:)*geom%len_ac(i) + r_a(:,i)
    end do

      
    ! rotate back atoms by -(sig(1) - sig1_init) around -ex
    sig1 = atan2(sin_sig(1), cos_sig(1))
    Us = rotation_matrix2(-ex, sig1_init-sig1)
     
    soln(:,1,1,i_soln) = r_anchor(:,1)
    soln(:,2,1,i_soln) = r_anchor(:,2)
    soln(:,3,1,i_soln) = matmul(Us, r_c(:,1) - r0(:)) + r0(:)
    soln(:,1,2,i_soln) = matmul(Us, r_n(:,2) - r0(:)) + r0(:)
    soln(:,2,2,i_soln) = matmul(Us, r_a(:,2) - r0(:)) + r0(:)
    soln(:,3,2,i_soln) = matmul(Us, r_c(:,2) - r0(:)) + r0(:)
    soln(:,1,3,i_soln) = matmul(Us, r_n(:,3) - r0(:)) + r0(:)
    soln(:,2,3,i_soln) = r_anchor(:,3)
    soln(:,3,3,i_soln) = r_anchor(:,4)
end do

end subroutine coord_from_poly_roots
!-------------------------------------------------------------------------------
subroutine poly_mul_sub2(u1, u2, u3, u4, p1, p2, p3, p4, u5, p5)
!-------------------------------------------------------------------------------
real(dp), dimension(0:4,0:4), intent(in) :: u1, u2, u3, u4
integer, dimension(2), intent(in) :: p1, p2, p3, p4
real(dp), dimension(0:4,0:4), intent(out) :: u5
integer, dimension(2), intent(out) :: p5
real(dp), dimension(0:4,0:4) :: d1, d2
integer, dimension(2) :: pd1, pd2

call poly_mul2(u1, u2, p1, p2, d1, pd1)
call poly_mul2(u3, u4, p3, p4, d2, pd2)
call poly_sub2(d1, d2, pd1, pd2, u5, p5)

end subroutine poly_mul_sub2
!-------------------------------------------------------------------------------
subroutine poly_mul2(u1, u2, p1, p2, u3, p3)
!-------------------------------------------------------------------------------
real(dp), dimension(0:4,0:4), intent(in) :: u1, u2
integer, dimension(2), intent(in) :: p1, p2
real(dp), dimension(0:4,0:4), intent(out) :: u3
integer, intent(out) :: p3(2)
integer :: i1, j1, i2, j2, i3, j3, p11, p12, p21, p22
real(dp) :: u1ij

p3(:) = p1(:) + p2(:)
u3(:,:) = 0.0d0

p11 = p1(1)
p12 = p1(2)
p21 = p2(1)
p22 = p2(2)

do i1 = 0, p12
    do j1 = 0, p11
        u1ij = u1(j1,i1)
        do i2 = 0, p22
            i3 = i1 + i2
            do j2 = 0, p21
                j3 = j1 + j2
                u3(j3,i3) = u3(j3,i3) + u1ij*u2(j2,i2)
            end do
        end do
    end do
end do

end subroutine poly_mul2
!-------------------------------------------------------------------------------
subroutine poly_sub2(u1, u2, p1, p2, u3, p3)
!-------------------------------------------------------------------------------
real(dp), dimension(0:4,0:4), intent(in) :: u1, u2
integer, intent(in) :: p1(2), p2(2)
real(dp), dimension(0:4,0:4), intent(out) :: u3
integer, intent(out) :: p3(2)
integer :: i, j, p11, p12, p21, p22
logical :: i1_ok, i2_ok

p11 = p1(1)
p12 = p1(2)
p21 = p2(1)
p22 = p2(2)
p3(1) = max(p11,p21)
p3(2) = max(p12,p22)

do i = 0, p3(2)
    i1_ok = (i > p12)
    i2_ok = (i > p22)
    do j = 0, p3(1)
        if (i2_ok .or. (j > p21)) then
            u3(j,i) = u1(j,i)
        else if (i1_ok .or. (j > p11)) then
            u3(j,i) = -u2(j,i)
        else
            u3(j,i) = u1(j,i) - u2(j,i)
        end if
    end do
end do

end subroutine poly_sub2
!-------------------------------------------------------------------------------
subroutine poly_mul_sub1(u1, u2, u3, u4, p1, p2, p3, p4, u5, p5)
!-------------------------------------------------------------------------------
real(dp), dimension(0:16), intent(in) :: u1, u2, u3, u4
integer, intent(in) :: p1, p2, p3, p4
real(dp), dimension(0:16), intent(out) :: u5
integer, intent(out) :: p5
real(dp), dimension(0:16) :: d1, d2
integer :: pd1, pd2

call poly_mul1(u1, u2, p1, p2, d1, pd1)
call poly_mul1(u3, u4, p3, p4, d2, pd2)
call poly_sub1(d1, d2, pd1, pd2, u5, p5)

end subroutine poly_mul_sub1
!-------------------------------------------------------------------------------
subroutine poly_mul1(u1, u2, p1, p2, u3, p3)
!-------------------------------------------------------------------------------
real(dp), dimension(0:16), intent(in) :: u1, u2
integer, intent(in) :: p1, p2
real(dp), dimension(0:16), intent(out) :: u3
integer, intent(out) :: p3
integer :: i1, i2, i3
real(dp) :: u1i

p3 = p1 + p2
u3(:) = 0.0d0

do i1 = 0, p1
    u1i = u1(i1)
    do i2 = 0, p2
        i3 = i1 + i2
        u3(i3) = u3(i3) + u1i*u2(i2)
    end do
end do

end subroutine poly_mul1
!-------------------------------------------------------------------------------
subroutine poly_sub1(u1, u2, p1, p2, u3, p3)
!-------------------------------------------------------------------------------
real(dp), dimension(0:16), intent(in) :: u1, u2
integer, intent(in) :: p1, p2
real(dp), dimension(0:16), intent(out) :: u3
integer, intent(out) :: p3
integer :: i

p3 = max(p1, p2)

do i = 0, p3
    if (i > p2) then
        u3(i) = u1(i)
    else if (i > p1) then
        u3(i) = -u2(i)
    else
        u3(i) = u1(i) - u2(i)
    end if
end do

end subroutine poly_sub1
!-------------------------------------------------------------------------------
function calc_t2(t0)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: t0
real(dp) :: calc_t2
real(dp) :: B0, B1, B2, A0, A1, A2, A3, A4, B2_2, B2_3
real(dp) :: K0, K1, K2, K3, t0_2, t0_3, t0_4

t0_2 = t0*t0
t0_3 = t0_2*t0
t0_4 = t0_3*t0

A0 = Q(0,0) + Q(1,0)*t0 + Q(2,0)*t0_2 + Q(3,0)*t0_3 + Q(4,0)*t0_4
A1 = Q(0,1) + Q(1,1)*t0 + Q(2,1)*t0_2 + Q(3,1)*t0_3 + Q(4,1)*t0_4
A2 = Q(0,2) + Q(1,2)*t0 + Q(2,2)*t0_2 + Q(3,2)*t0_3 + Q(4,2)*t0_4
A3 = Q(0,3) + Q(1,3)*t0 + Q(2,3)*t0_2 + Q(3,3)*t0_3 + Q(4,3)*t0_4
A4 = Q(0,4) + Q(1,4)*t0 + Q(2,4)*t0_2 + Q(3,4)*t0_3 + Q(4,4)*t0_4

B0 = RM(0,0) + RM(1,0)*t0 + RM(2,0)*t0_2
B1 = RM(0,1) + RM(1,1)*t0 + RM(2,1)*t0_2
B2 = RM(0,2) + RM(1,2)*t0 + RM(2,2)*t0_2

B2_2 = B2*B2
B2_3 = B2_2*B2

K0 = A2*B2 - A4*B0
K1 = A3*B2 - A4*B1
K2 = A1*B2_2 - K1*B0
K3 = K0*B2 - K1*B1
  
calc_t2 = (K3*B0 - A0*B2_3)/(K2*B2 - K3*B1)

end function calc_t2
!-------------------------------------------------------------------------------
function calc_t1(t0, t2)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: t0, t2
real(dp) :: calc_t1
real(dp) :: U11, U12, U13, U31, U32, U33
real(dp) :: t0_2, t2_2

t0_2 = t0*t0
t2_2 = t2*t2

U11 = C0(0,1) + C0(1,1)*t0 + C0(2,1)*t0_2
U12 = C1(0,1) + C1(1,1)*t0 + C1(2,1)*t0_2
U13 = C2(0,1) + C2(1,1)*t0 + C2(2,1)*t0_2
U31 = C0(0,2) + C0(1,2)*t2 + C0(2,2)*t2_2
U32 = C1(0,2) + C1(1,2)*t2 + C1(2,2)*t2_2
U33 = C2(0,2) + C2(1,2)*t2 + C2(2,2)*t2_2

calc_t1 = (U31*U13-U11*U33)/(U12*U33-U13*U32)

end function calc_t1
!-------------------------------------------------------------------------------
END MODULE LOOP_CLOSURE
!-------------------------------------------------------------------------------
