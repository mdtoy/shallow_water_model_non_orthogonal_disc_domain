module prognostics

!-----------------------------------------------------------------
!
!   Contains arrays related to the prognostic variables, and
!       those diagnostic variables derived from the prognostic
!       variables.  Variables are initialized to zero here.
!       Surface geopotential also declared here.
!
!-----------------------------------------------------------------


use kinds
use model_parameters
use physical_parameters


implicit none
save


!
! Declare prognostic variables
!
real (kind = dbl_kind), dimension(im,jm,ntprog) ::          &
                     sqrt_G_h_star   ! sqrt_G*thickness of fluid (m) sqrt_G*(h-hs)
real (kind = dbl_kind), dimension(im+1,jm,ntprog) ::        &
                     u1_cov          ! x1-covariant velocity components (m/s)
real (kind = dbl_kind), dimension(im,jm+1,ntprog) ::        &
                     u2_cov          ! x2-covariant velocity components (m/s)

!
! Declare tendencies of prognostic variables
!
real (kind = dbl_kind), dimension(im+1,jm,nttend) ::          &
                     u1_cov_f         ! d/dt (u1_cov)   (m/s^2)
real (kind = dbl_kind), dimension(im,jm+1,nttend) ::          &
                     u2_cov_f         ! d/dt (u2_cov)   (m/s^2)
real (kind = dbl_kind), dimension(im,jm,nttend) ::            &
                     sqrt_G_h_star_f      ! d/dt (sqrt_G*h_star)   (m/s)

!
! Declare diagnostic variables
!
real (kind = dbl_kind), dimension(im,jm) ::                 &
                     ke_horiz,  &     ! contribution to kinetic energy
                                      ! from horizontal velocity (J/kg)
                     h,         &     ! height of free surface (m)
                     h_star           ! thickness of fluid (m)

real (kind = dbl_kind), dimension(im+1,jm) ::               &
                     u_init_u1, &     ! initial Cartesian u and v (m/s)
                     v_init_u1, &     ! at u1 points
                     u1_cont,   &     ! contravariant x1 velocity component
                     u                ! Cartesian u velocity (m/s)

real (kind = dbl_kind), dimension(im,jm+1) ::               &
                     u_init_u2, &     ! initial Cartesian u and v (m/s)
                     v_init_u2, &     ! at u2 points
                     u2_cont,   &     ! contravariant x2 velocity component
                     v                ! Cartesian v velocity (m/s)

real (kind = dbl_kind), dimension(im+1,jm+1) ::             &
                     zeta,      &     ! relative vorticity (s^-1)
                     pv               ! vertical component of potential 
                                      ! vorticity (m^2 eta/kg/s)


!
! Declare topographic height
!
real (kind = dbl_kind), dimension (im,jm) ::                &
                     hs             ! topographic height (m)






contains


!======================================================================
! BEGINNING OF INIT_PROGNOSTICS
!======================================================================

subroutine init_prognostics

implicit none

! initialize prognostic arrays
u1_cov = c0
u2_cov = c0
sqrt_G_h_star = c0

u1_cov_f = c0
u2_cov_f = c0
sqrt_G_h_star_f = c0


end subroutine init_prognostics

!======================================================================
! END OF INIT_PROGNOSTICS
!======================================================================


end module prognostics
