module momentum_tendency

!-----------------------------------------------------------------------
! PURPOSE: Calculates the tendencies of u and v.
!-----------------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters


implicit none
save



contains



!======================================================================
! BEGINNING OF GET_U1COVF_U2COVF
!======================================================================

subroutine get_u1covf_u2covf ( u1_cov, u2_cov, u1_cont, u2_cont,     &
                 F_u1, F_u2, h, h_star, sqrt_G_h_star,               &
                 pv, zeta, ke_horiz, u1_cov_f, u2_cov_f )

!---------------------------------------------------------------------------
! PURPOSE:
!   Computes the time tendency of the x and y components of velocity,
!   i.e. u and v respectively.
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT IN
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im+1,jm), intent(in) ::            &
          u1_cov,   &           ! covariant velocity components (m/s)
          u1_cont               ! contravariant velocity components (m/s)
real (kind=dbl_kind), dimension(im,jm+1), intent(in) ::            &
          u2_cov,   &           ! covariant velocity components (m/s)
          u2_cont               ! contravariant velocity components (m/s)

real (kind=dbl_kind), dimension(im+1,jm), intent(in) ::            &
          F_u1           ! Mass flux in u1 direction
real (kind=dbl_kind), dimension(im,jm+1), intent(in) ::            &
          F_u2           ! Mass flux in u2 direction
real (kind=dbl_kind), dimension(im,jm), intent(in) ::              &
          h,           & ! height of free surface (m)
          h_star,      & ! thickness of fluid (m) (h-hs)
          sqrt_G_h_star  ! sqrt_G times h_star (m)

!---------------------------------------------------------------------------
! INTENT OUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im+1,jm+1), intent(out) ::         &
          pv,      &     ! vert. component of pot. vorticity (m-1 s-1)
          zeta           ! relative vorticity (s^-1)
real (kind=dbl_kind), dimension(im,jm), intent(out) ::             &
          ke_horiz       ! contribution to kinetic energy
                         ! from horizontal velocity (J/kg)

real (kind=dbl_kind), dimension(im+1,jm), intent(out) ::           &
          u1_cov_f       ! tendency of u1_cov (m/s^2)
real (kind=dbl_kind), dimension(im,jm+1), intent(out) ::           &
          u2_cov_f       ! tendency of u2_cov (m/s^2)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: i,j

real (kind=dbl_kind), parameter ::                                   &
          inv24 = c1/24.00000_dbl_kind

real (kind=dbl_kind), dimension(im+1,jm+1) ::                        &
          abs_vort,        &                  ! absolute vorticity (s^-1)
          h_q                                 ! SW height interpolated to q-pts (m)
real (kind=dbl_kind), dimension(1:im,jm) ::                          &
          alfa, delt                   ! linear combinations of pv
real (kind=dbl_kind), dimension(2:im+1,jm) ::                        &
          beta, gamm                   ! linear combinations of pv
real (kind=dbl_kind), dimension(im,jm) ::                            &
          epsln, fi                    ! linear combinations of pv
real (kind=dbl_kind) :: temp1, temp2



!---------------------------------------------------------------------------
! Calculate the (PV k) cross (m times velocity) term.       (term 1 of 5)
! Potential enstrophy and energy conserving scheme of Arakawa and 
! Lamb (1981) is used.
!---------------------------------------------------------------------------

! calculate potential vorticity

! Interior points
do j = 2,jm
   do i = 2,im
      zeta(i,j) = inv_sqrt_G_q(i,j) *                                &
                  ( invdx2 * ( u1_cov(i,j-1) - u1_cov(i,j) )  +      &
                    invdx1 * ( u2_cov(i,j) - u2_cov(i-1,j) ) )
   end do
end do
! Left edge
i=1
do j = 2,jm
   zeta(i,j) = inv_sqrt_G_q(i,j) *                                   &
      ( invdx2 * p5*( u1_cov(1,j-1) - u1_cov(1,j) +                  &
                      u1_cov(im+1,j-1) - u1_cov(im+1,j) ) +          &
        invdx1 * ( u2_cov(1,j) - u2_cov(im,j) ) )
   ! Right edge
   zeta(im+1,j) = zeta(i,j)
end do
! Bottom edge
j=1
do i = 2,im
   zeta(i,j) = inv_sqrt_G_q(i,j) *                                   &
         ( invdx2 * ( u1_cov(i,jm) - u1_cov(i,1) )  +                &
           invdx1 * p5*( u2_cov(i,1) - u2_cov(i-1,1) +               &
                         u2_cov(i,jm+1) - u2_cov(i-1,jm+1) )  )
   ! Top edge
   zeta(i,jm+1) = zeta(i,j)
end do
! Corners
zeta(1,1) = inv_sqrt_G_q(1,1) *                                      &
     ( invdx2 * p5*( u1_cov(1,jm) - u1_cov(1,1) +                    &
                     u1_cov(im+1,jm) - u1_cov(im+1,1) )  +           &
       invdx1 * p5*( u2_cov(1,1) - u2_cov(im,1) +                    &
                     u2_cov(1,jm+1) - u2_cov(im,jm+1) ) )
zeta(1,jm+1) = zeta(1,1)
zeta(im+1,1) = zeta(1,1)
zeta(im+1,jm+1) = zeta(1,1)

abs_vort(:,:) = f_cor(:,:) + zeta(:,:)

! Calculate mass interpolated to q-points
! Interior points
do j = 2,jm
   do i = 2,im
      h_q(i,j) = inv_sqrt_G_q(i,j) * p25 * ( sqrt_G_h_star(i,j) +    &
                  sqrt_G_h_star(i-1,j) + sqrt_G_h_star(i-1,j-1) +    &
                  sqrt_G_h_star(i,j-1) )
   end do
end do
! Left edge
i=1
do j = 2,jm
   h_q(i,j) = inv_sqrt_G_q(i,j) * p25 * ( sqrt_G_h_star(i,j) +       &
                 sqrt_G_h_star(im,j) + sqrt_G_h_star(im,j-1) +       &
                 sqrt_G_h_star(i,j-1) )
   ! Right edge
   h_q(im+1,j) = h_q(i,j)
end do
! Bottom edge
j=1
do i = 2,im
   h_q(i,j) = inv_sqrt_G_q(i,j) * p25 * ( sqrt_G_h_star(i,j) +       &
                sqrt_G_h_star(i-1,j) + sqrt_G_h_star(i-1,jm) +       &
                sqrt_G_h_star(i,jm) )
   ! Top edge
   h_q(i,jm+1) = h_q(i,j)
end do
! Corners
h_q(1,1) = inv_sqrt_G_q(1,1) * p25 * ( sqrt_G_h_star(1,1) +       &
                  sqrt_G_h_star(im,1) + sqrt_G_h_star(im,jm) +    &
                  sqrt_G_h_star(1,jm) )
h_q(1,jm+1) = h_q(1,1)
h_q(jm+1,1) = h_q(1,1)
h_q(jm+1,im+1) = h_q(1,1)


pv(:,:) = abs_vort(:,:) / h_q(:,:)


! calculate linear combinations of pv ( eqn. (3.34) of AL (1981) )
! Interior points
do j = 1,jm
   do i = 2,im
      alfa(i,j)  = inv24 *                                           &
            ( c2*pv(i+1,j+1) + pv(i,j+1) + c2*pv(i,j) + pv(i+1,j) )
      beta(i,j)  = inv24 *                                           &
            ( pv(i,j+1) + c2*pv(i-1,j+1) + pv(i-1,j) + c2*pv(i,j) )
      gamm(i,j)  = inv24 *                                           &
            ( c2*pv(i,j+1) + pv(i-1,j+1) + c2*pv(i-1,j) + pv(i,j) )
      delt(i,j)  = inv24 *                                           &
            ( pv(i+1,j+1) + c2*pv(i,j+1) + pv(i,j) + c2*pv(i+1,j) )
   end do
   ! Left edge
   i=1
   alfa(i,j)  = inv24 *                                              &
          ( c2*pv(i+1,j+1) + pv(i,j+1) + c2*pv(i,j) + pv(i+1,j) )
   delt(i,j)  = inv24 *                                              &
          ( pv(i+1,j+1) + c2*pv(i,j+1) + pv(i,j) + c2*pv(i+1,j) )
   ! Right edge
   i=im+1
   beta(i,j)  = inv24 *                                              &
          ( pv(i,j+1) + c2*pv(i-1,j+1) + pv(i-1,j) + c2*pv(i,j) )
   gamm(i,j)  = inv24 *                                              &
          ( c2*pv(i,j+1) + pv(i-1,j+1) + c2*pv(i-1,j) + pv(i,j) )
end do
do j = 1,jm
   do i = 1,im
      epsln(i,j) = inv24 *                                           &
           ( pv(i+1,j+1) + pv(i,j+1) - pv(i,j) - pv(i+1,j) )
      fi(i,j)    = inv24 *                                           &
           (-pv(i+1,j+1) + pv(i,j+1) + pv(i,j) - pv(i+1,j) ) 
   end do
end do


!
! calculate u and v tendencies -- term 1 of 5 ( see eqns. (3.5) and (3.6)
!                                 of AL (1981) )

! Interior points
do j = 1,jm
   do i = 2,im
      u1_cov_f(i,j) = alfa(i,j)*F_u2(i,j+1) + beta(i,j)*F_u2(i-1,j+1) +    &
                      gamm(i,j)*F_u2(i-1,j) + delt(i,j)*F_u2(i,j) -        &
                      epsln(i,j)*F_u1(i+1,j) + epsln(i-1,j)*F_u1(i-1,j)
   end do
   ! Left edge
   i=1
   temp1 = c2*alfa(i,j)*F_u2(i,j+1) + c2*delt(i,j)*F_u2(i,j) -       &
                   epsln(i,j)*F_u1(i+1,j) + epsln(im,j)*F_u1(im,j)
   ! Right edge
   i=im+1
   temp2 = c2*beta(i,j)*F_u2(i-1,j+1) + c2*gamm(i,j)*F_u2(i-1,j) -   &
                   epsln(1,j)*F_u1(2,j) + epsln(im,j)*F_u1(im,j)
   ! Change in plan -- use local coordinate at edge
   u1_cov_f(1,j) = p5*(temp1+temp2)
   u1_cov_f(im+1,j) = u1_cov_f(1,j)
end do

! Interior points
do i = 1,im
   do j = 2,jm
      u2_cov_f(i,j) = -gamm(i+1,j)*F_u1(i+1,j) - delt(i,j)*F_u1(i,j) -  &
           alfa(i,j-1)*F_u1(i,j-1) - beta(i+1,j-1)*F_u1(i+1,j-1) -      &
           fi(i,j)*F_u2(i,j+1) + fi(i,j-1)*F_u2(i,j-1)
   end do
   ! Bottom edge
   j=1
   temp1 = -c2*gamm(i+1,j)*F_u1(i+1,j) - c2*delt(i,j)*F_u1(i,j) -    &
                    fi(i,j)*F_u2(i,j+1) + fi(i,jm)*F_u2(i,jm)
   
   ! Top edge
   j=jm+1
   temp2 = - c2*alfa(i,j-1)*F_u1(i,j-1) -                            &
                     c2*beta(i+1,j-1)*F_u1(i+1,j-1) -                &
                     fi(i,1)*F_u2(i,2) + fi(i,jm)*F_u2(i,jm)
   ! Change in plan -- use local coordinate at edge
   u2_cov_f(i,1) = p5*(temp1+temp2)
   u2_cov_f(i,jm+1) = u2_cov_f(i,1)
end do


!---------------------------------------------------------------------------
! Add contribution of the horiz. gradient of kinetic energy.   (term 2 of 5)
! And also divergence damping!!!
!---------------------------------------------------------------------------


! Calculate contribution to kinetic energy
! from the horizontal velocity  ( see eqn (3.41) of AL (1981) )
! Note:  expect SICK to result from use of this K.E.

! Interior points
do j = 1,jm
   do i = 1,im
      ke_horiz(i,j) = p25 * ( u1_cov(i,j)*u1_cont(i,j) +             &
                              u1_cov(i+1,j)*u1_cont(i+1,j) ) +       &
                      p25 * ( u2_cov(i,j)*u2_cont(i,j) +             &
                              u2_cov(i,j+1)*u2_cont(i,j+1) )
   end do
end do



! Add contribution of horiz. gradient of K.E.
! Interior points
do j = 1,jm
   do i = 2,im
      u1_cov_f(i,j) = u1_cov_f(i,j) - invdx1*(ke_horiz(i,j)-ke_horiz(i-1,j))
   end do
   ! Left edge
   i=1
   u1_cov_f(i,j) = u1_cov_f(i,j) - invdx1*(ke_horiz(i,j)-ke_horiz(im,j))
   ! Right edge
   i=im+1
   u1_cov_f(i,j) = u1_cov_f(i,j) - invdx1*(ke_horiz(1,j)-ke_horiz(im,j))
end do
! Interior points
do i = 1,im
   do j = 2,jm
      u2_cov_f(i,j) = u2_cov_f(i,j) - invdx2*(ke_horiz(i,j)-ke_horiz(i,j-1))
   end do
   ! Bottom edge
   j=1
   u2_cov_f(i,j) = u2_cov_f(i,j) - invdx2*(ke_horiz(i,j)-ke_horiz(i,jm))
   ! Top edge
   j=jm+1
   u2_cov_f(i,j) = u2_cov_f(i,j) - invdx2*(ke_horiz(i,1)-ke_horiz(i,jm))
end do


!---------------------------------------------------------------------------
! Add contributions due to vert. advection of horiz. momentum and
! subgrid-scale turbulent momentum flux.            (terms 3 & 4 of 5)
! XXXXXXXXXXXXXXX NOT APPLICABLE TO SHALLOW WATER EQUATIONS XXXXXXXXXXXX
!---------------------------------------------------------------------------


!---------------------------------------------------------------------------
! Add contribution of the horizontal pressure gradient force.  (term 5 of 5)
!---------------------------------------------------------------------------

! Interior points
do j = 1,jm
   do i = 2,im
      u1_cov_f(i,j) = u1_cov_f(i,j) - invdx1*grav*(h(i,j)-h(i-1,j))
   end do
   ! Left edge
   i=1
   u1_cov_f(i,j) = u1_cov_f(i,j) - invdx1*grav*(h(i,j)-h(im,j))
   ! Right edge
   i=im+1
   u1_cov_f(i,j) = u1_cov_f(i,j) - invdx1*grav*(h(1,j)-h(im,j))
end do
! Interior points
do i = 1,im
   do j = 2,jm
      u2_cov_f(i,j) = u2_cov_f(i,j) - invdx2*grav*(h(i,j)-h(i,j-1))
   end do
   ! Bottom edge
   j=1
   u2_cov_f(i,j) = u2_cov_f(i,j) - invdx2*grav*(h(i,j)-h(i,jm))
   ! Top edge
   j=jm+1
   u2_cov_f(i,j) = u2_cov_f(i,j) - invdx2*grav*(h(i,1)-h(i,jm))
end do

end subroutine get_u1covf_u2covf

!======================================================================
! END OF GET_UF_VF
!======================================================================



end module momentum_tendency
