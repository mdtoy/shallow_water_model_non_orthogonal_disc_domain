module initialize

!-----------------------------------------------------------------
!   This module sets the initial conditions of the model and
!   and sets the values of the eta levels.
!-----------------------------------------------------------------

use kinds
use physical_parameters
use model_parameters
use prognostics
use step

implicit none
save


contains


!===================================================================
! BEGINNING OF INITIALIZE_MODEL
!===================================================================

subroutine initialize_model(step_count,tau_initial)

implicit none

!-------------------------------------------------------------------
! INTENT IN 
!-------------------------------------------------------------------
integer (kind = int_kind), intent(in) :: step_count

!-------------------------------------------------------------------
! INTENT INOUT
!-------------------------------------------------------------------
real (kind = dbl_kind), intent(inout) :: tau_initial

!-------------------------------------------------------------------
! LOCAL
!-------------------------------------------------------------------
integer :: i, j

logical, parameter :: init_from_ic_prog_file = .false.
logical, parameter :: init_from_bc_surf_file = .false.
logical, parameter :: water_mtn = .false.
logical, parameter :: conical_mtn = .true.

real (kind=dbl_kind), parameter ::    &
     r0 = (pi/9._dbl_kind)*a_rad,   &  ! radius of mountain (cone-shaped TC5 mountain)
     r0_G = (p5/log(c2)**0.5)*r0,    &  ! radius of Gaussian mountain (matches half-width of TC5 mtn)
     x_ctr = -p5*pi*a_rad,    &         ! x-coord center of mountain
     y_ctr = (pi/6._dbl_kind)*a_rad,   &   ! y-coord center of mountain
     mtn_ht = 2000._dbl_kind,  &    !  2000._dbl_kind     ! height of topo. mountain (m)
     water_mtn_ht = 2000._dbl_kind     ! height of water mountain (m)
real (kind=dbl_kind) :: r_test     ! test to see if within perturbation

real (kind=dbl_kind), dimension(im,    jm) :: Jtemp_h, Ktemp_h, Ltemp_h, Mtemp_h
real (kind=dbl_kind), dimension(im+1,  jm) :: Jtemp_u, Ktemp_u, Ltemp_u, Mtemp_u
real (kind=dbl_kind), dimension(im,  jm+1) :: Jtemp_v, Ktemp_v, Ltemp_v, Mtemp_v
real (kind=dbl_kind), dimension(im+1,jm+1) :: Jtemp_q, Ktemp_q, Ltemp_q, Mtemp_q

! Components of the transformation Jacobian along the edge q-points only
real (kind=dbl_kind), dimension(1:im+1) ::                                           &
   dx1_dx_q_top_edge, dx_dx2_q_bottom_edge, dx1_dy_q_top_edge, dy_dx2_q_bottom_edge, &
   dx1_dx_q_bottom_edge, dx_dx2_q_top_edge, dx1_dy_q_bottom_edge, dy_dx2_q_top_edge
real (kind=dbl_kind), dimension(1:jm+1) ::                                           &
   dx2_dx_q_right_edge, dx_dx1_q_left_edge, dx2_dy_q_right_edge, dy_dx1_q_left_edge, &
   dx2_dx_q_left_edge, dx_dx1_q_right_edge, dx2_dy_q_left_edge, dy_dx1_q_right_edge

! Cross-edge metric tensors
real (kind=dbl_kind), dimension (1:im+1) ::    &
         a1_dot_a2_top_edge, a1_dot_a2_bottom_edge
real (kind=dbl_kind), dimension (1:jm+1) ::    &
         a2_dot_a1_left_edge, a2_dot_a1_right_edge



call init_prognostics



! Initialize prognostic and tendency time step indices
n4 = mod(step_count+1,2) + 1
n3 = mod(step_count,2) + 1
n3_f = mod(step_count+2,3) + 1
n2_f = mod(step_count+1,3) + 1
n1_f = mod(step_count,3) + 1




! Initialize coordinate x1-x2 locations
do i = 1,im+1
   x1_u1(i) = (i-1)*dx1 - p5*im*dx1 + alf*a_rad  ! origin at center of domain
end do
do i = 1,im
   x1_h(i) = (i-1)*dx1 + p5*dx1 - p5*im*dx1 + alf*a_rad ! origin at center of domain
end do
do j = 1,jm+1
   x2_u2(j) = (j-1)*dx2 - p5*jm*dx2 + bet*a_rad  ! origin at center of domain
end do
do j = 1,jm
   x2_h(j) = (j-1)*dx2 + p5*dx2 - p5*jm*dx2 + bet*a_rad ! origin at center of domain
end do



! Calculate Cartesian-coordinate (x-y) locations
do j = 1,jm
   do i = 1,im
      call coordinate_solve ( x1_h(i),  x2_h(j),  x_h(i,j), y_h(i,j) )  ! mass points
   end do
end do
do j = 1,jm
   do i = 1,im+1
      call coordinate_solve ( x1_u1(i), x2_h(j),  x_u(i,j), y_u(i,j) )  ! u1 points
   end do
end do
do j = 1,jm+1
   do i = 1,im
      call coordinate_solve ( x1_h(i),  x2_u2(j), x_v(i,j), y_v(i,j) )  ! u2 points
   end do
end do
do j = 1,jm+1
   do i = 1,im+1
      call coordinate_solve ( x1_u1(i), x2_u2(j), x_q(i,j), y_q(i,j) )  ! vorticity points
   end do
end do



! Now that we have Cartesian-coordinate locations,
! we can initialize the metric terms
Jtemp_h(:,:) = c1 + C_fac*Cos(p5*(y_h(:,:)*inv_a_rad-bet))*          &
                          Cos(x_h(:,:)*inv_a_rad-alf)
Ktemp_h(:,:) = -p5*C_fac*Sin(x_h(:,:)*inv_a_rad-alf)*                &
                         Sin(p5*(y_h(:,:)*inv_a_rad-bet))
Ltemp_h(:,:) = c1 + D_fac*Cos(p5*(x_h(:,:)*inv_a_rad-alf))*          &
                          Cos(y_h(:,:)*inv_a_rad-bet)
Mtemp_h(:,:) = -p5*D_fac*Sin(y_h(:,:)*inv_a_rad-bet)*                &
                         Sin(p5*(x_h(:,:)*inv_a_rad-alf))
inv_sqrt_G_h(:,:) = Jtemp_h(:,:)*Ltemp_h(:,:) - Ktemp_h(:,:)*Mtemp_h(:,:)
sqrt_G_h(:,:) = c1 / inv_sqrt_G_h(:,:)

Jtemp_q(:,:) = c1 + C_fac*Cos(p5*(y_q(:,:)*inv_a_rad-bet))*          &
                          Cos(x_q(:,:)*inv_a_rad-alf)
Ktemp_q(:,:) = -p5*C_fac*Sin(x_q(:,:)*inv_a_rad-alf)*                &
                         Sin(p5*(y_q(:,:)*inv_a_rad-bet))
Ltemp_q(:,:) = c1 + D_fac*Cos(p5*(x_q(:,:)*inv_a_rad-alf))*          &
                          Cos(y_q(:,:)*inv_a_rad-bet)
Mtemp_q(:,:) = -p5*D_fac*Sin(y_q(:,:)*inv_a_rad-bet)*                &
                         Sin(p5*(x_q(:,:)*inv_a_rad-alf))
inv_sqrt_G_q(:,:) = Jtemp_q(:,:)*Ltemp_q(:,:) - Ktemp_q(:,:)*Mtemp_q(:,:)

Jtemp_u(:,:) = c1 + C_fac*Cos(p5*(y_u(:,:)*inv_a_rad-bet))*          &
                          Cos(x_u(:,:)*inv_a_rad-alf)
Ktemp_u(:,:) = -p5*C_fac*Sin(x_u(:,:)*inv_a_rad-alf)*                &
                         Sin(p5*(y_u(:,:)*inv_a_rad-bet))
Ltemp_u(:,:) = c1 + D_fac*Cos(p5*(x_u(:,:)*inv_a_rad-alf))*          &
                          Cos(y_u(:,:)*inv_a_rad-bet)
Mtemp_u(:,:) = -p5*D_fac*Sin(y_u(:,:)*inv_a_rad-bet)*                &
                         Sin(p5*(x_u(:,:)*inv_a_rad-alf))
G_cont_11_u1(:,:) = Jtemp_u(:,:)**2 + Ktemp_u(:,:)**2
G_cont_12_u1(:,:) = Jtemp_u(:,:)*Mtemp_u(:,:) + Ktemp_u(:,:)*Ltemp_u(:,:)

Jtemp_v(:,:) = c1 + C_fac*Cos(p5*(y_v(:,:)*inv_a_rad-bet))*          &
                          Cos(x_v(:,:)*inv_a_rad-alf)
Ktemp_v(:,:) = -p5*C_fac*Sin(x_v(:,:)*inv_a_rad-alf)*                &
                         Sin(p5*(y_v(:,:)*inv_a_rad-bet))
Ltemp_v(:,:) = c1 + D_fac*Cos(p5*(x_v(:,:)*inv_a_rad-alf))*          &
                          Cos(y_v(:,:)*inv_a_rad-bet)
Mtemp_v(:,:) = -p5*D_fac*Sin(y_v(:,:)*inv_a_rad-bet)*                &
                         Sin(p5*(x_v(:,:)*inv_a_rad-alf))
G_cont_22_u2(:,:) = Mtemp_v(:,:)**2 + Ltemp_v(:,:)**2
G_cont_12_u2(:,:) = Jtemp_v(:,:)*Mtemp_v(:,:) + Ktemp_v(:,:)*Ltemp_v(:,:)

! Initialize components of the transformation Jacobian
dx1_dx_u1(:,:) = Jtemp_u(:,:)
dx1_dy_u2(:,:) = Ktemp_v(:,:)
dx2_dx_u1(:,:) = Mtemp_u(:,:)
dx2_dy_u2(:,:) = Ltemp_v(:,:)
dx_dx1_u1(:,:) = Ltemp_u(:,:) / ( Ltemp_u(:,:)*Jtemp_u(:,:) -        &
                                      Ktemp_u(:,:)*Mtemp_u(:,:) )
dx_dx2_u2(:,:) = - Ktemp_v(:,:) / ( Ltemp_v(:,:)*Jtemp_v(:,:) -      &
                                      Ktemp_v(:,:)*Mtemp_v(:,:) )
dy_dx1_u1(:,:) = - Mtemp_u(:,:) / ( Ltemp_u(:,:)*Jtemp_u(:,:) -      &
                                      Ktemp_u(:,:)*Mtemp_u(:,:) )
dy_dx2_u2(:,:) = Jtemp_v(:,:) / ( Ltemp_v(:,:)*Jtemp_v(:,:) -        &
                                      Ktemp_v(:,:)*Mtemp_v(:,:) )

! Edge values
dx1_dx_q_top_edge(:) = Jtemp_q(1:im+1,jm+1)
dx1_dx_q_bottom_edge(:) = Jtemp_q(1:im+1,1)
dx1_dy_q_top_edge(:) = Ktemp_q(1:im+1,jm+1)
dx1_dy_q_bottom_edge(:) = Ktemp_q(1:im+1,1)
dx_dx2_q_top_edge(:) = - Ktemp_q(1:im+1,jm+1) / ( Ltemp_q(1:im+1,jm+1)*   &
     Jtemp_q(1:im+1,jm+1) - Ktemp_q(1:im+1,jm+1)*Mtemp_q(1:im+1,jm+1) )
dx_dx2_q_bottom_edge(:) = - Ktemp_q(1:im+1,1) / ( Ltemp_q(1:im+1,1)*      &
     Jtemp_q(1:im+1,1) - Ktemp_q(1:im+1,1)*Mtemp_q(1:im+1,1) )
dy_dx2_q_top_edge(:) = Jtemp_q(1:im+1,jm+1) / ( Ltemp_q(1:im+1,jm+1)*     &
     Jtemp_q(1:im+1,jm+1) - Ktemp_q(1:im+1,jm+1)*Mtemp_q(1:im+1,jm+1) )
dy_dx2_q_bottom_edge(:) = Jtemp_q(1:im+1,1) / ( Ltemp_q(1:im+1,1)*        &
     Jtemp_q(1:im+1,1) - Ktemp_q(1:im+1,1)*Mtemp_q(1:im+1,1) )
dx2_dx_q_right_edge(:) = Mtemp_q(im+1,1:jm+1)
dx_dx1_q_left_edge(:) = Ltemp_q(1,1:jm+1) / ( Ltemp_q(1,1:jm+1)*          &
     Jtemp_q(1,1:jm+1) - Ktemp_q(1,1:jm+1)*Mtemp_q(1,1:jm+1) )
dx2_dy_q_right_edge(:) = Ltemp_q(im+1,1:jm+1)
dy_dx1_q_left_edge(:) = - Mtemp_q(1,1:jm+1) / ( Ltemp_q(1,1:jm+1)*        &
     Jtemp_q(1,1:jm+1) - Ktemp_q(1,1:jm+1)*Mtemp_q(1,1:jm+1) )
dx2_dx_q_left_edge(:) = Mtemp_q(1,1:jm+1)
dx_dx1_q_right_edge(:) = Ltemp_q(im+1,1:jm+1) / ( Ltemp_q(im+1,1:jm+1)*   &
     Jtemp_q(im+1,1:jm+1) - Ktemp_q(im+1,1:jm+1)*Mtemp_q(im+1,1:jm+1) )
dx2_dy_q_left_edge(:) = Ltemp_q(1,1:jm+1)
dy_dx1_q_right_edge(:) = - Mtemp_q(im+1,1:jm+1) / ( Ltemp_q(im+1,1:jm+1)* &
     Jtemp_q(im+1,1:jm+1) - Ktemp_q(im+1,1:jm+1)*Mtemp_q(im+1,1:jm+1) )


! Cross-edge metric tensors
a1_dot_a2_top_edge(:) = dx1_dx_q_bottom_edge(:)*dx_dx2_q_top_edge(:) +    &
                        dx1_dy_q_bottom_edge(:)*dy_dx2_q_top_edge(:)
a1_dot_a2_bottom_edge(:) = dx1_dx_q_top_edge(:)*dx_dx2_q_bottom_edge(:) + &
                           dx1_dy_q_top_edge(:)*dy_dx2_q_bottom_edge(:)
a2_dot_a1_left_edge(:) = dx2_dx_q_right_edge(:)*dx_dx1_q_left_edge(:) +   &
                         dx2_dy_q_right_edge(:)*dy_dx1_q_left_edge(:)
a2_dot_a1_right_edge(:) = dx2_dx_q_left_edge(:)*dx_dx1_q_right_edge(:) +  &
                          dx2_dy_q_left_edge(:)*dy_dx1_q_right_edge(:)


! Calculate modified G_cont_11 and G_cont_22 metric tensors near edges
G_cont_11_u1_tilde_top(:) = G_cont_11_u1(:,jm) + p25*G_cont_12_u1(:,jm)*    &
                                                    a1_dot_a2_top_edge(:)
G_cont_11_u1_tilde_bottom(:) = G_cont_11_u1(:,1) + p25*G_cont_12_u1(:,1)*   &
                                                    a1_dot_a2_bottom_edge(:)
G_cont_22_u2_tilde_left(:) = G_cont_22_u2(1,:) + p25*G_cont_12_u2(1,:)*     &
                                                    a2_dot_a1_left_edge(:)
G_cont_22_u2_tilde_right(:) = G_cont_22_u2(im,:) + p25*G_cont_12_u2(im,:)*  &
                                                    a2_dot_a1_right_edge(:)



! Initialize G_cont_12 weightings for averaging contravariant/covariant velocities
! Arithmetic mean now used
! Interior u1-points
do j = 1,jm
   do i = 2,im
      mean_G_cont_12_u1_nw(i,j) = p5*(G_cont_12_u1(i,j)+G_cont_12_u2(i-1,j+1))
      mean_G_cont_12_u1_ne(i,j) = p5*(G_cont_12_u1(i,j)+G_cont_12_u2(i,j+1))
      mean_G_cont_12_u1_sw(i,j) = p5*(G_cont_12_u1(i,j)+G_cont_12_u2(i-1,j))
      mean_G_cont_12_u1_se(i,j) = p5*(G_cont_12_u1(i,j)+G_cont_12_u2(i,j))
   end do
   ! Left edge
   i=1
   mean_G_cont_12_u1_nw(i,j) = p5*(G_cont_12_u1(im+1,j)+G_cont_12_u2(im,j+1))
   mean_G_cont_12_u1_ne(i,j) = p5*(G_cont_12_u1(i,j)+G_cont_12_u2(i,j+1))
   mean_G_cont_12_u1_sw(i,j) = p5*(G_cont_12_u1(im+1,j)+G_cont_12_u2(im,j))
   mean_G_cont_12_u1_se(i,j) = p5*(G_cont_12_u1(i,j)+G_cont_12_u2(i,j))
   ! Right edge
   mean_G_cont_12_u1_nw(im+1,j) = mean_G_cont_12_u1_nw(i,j)
   mean_G_cont_12_u1_ne(im+1,j) = mean_G_cont_12_u1_ne(i,j)
   mean_G_cont_12_u1_sw(im+1,j) = mean_G_cont_12_u1_sw(i,j)
   mean_G_cont_12_u1_se(im+1,j) = mean_G_cont_12_u1_se(i,j)
end do

! Interior u2-points
do i = 1,im
   do j = 2,jm
      mean_G_cont_12_u2_nw(i,j) = p5*(G_cont_12_u2(i,j)+G_cont_12_u1(i,j))
      mean_G_cont_12_u2_ne(i,j) = p5*(G_cont_12_u2(i,j)+G_cont_12_u1(i+1,j))
      mean_G_cont_12_u2_sw(i,j) = p5*(G_cont_12_u2(i,j)+G_cont_12_u1(i,j-1))
      mean_G_cont_12_u2_se(i,j) = p5*(G_cont_12_u2(i,j)+G_cont_12_u1(i+1,j-1))
   end do
   ! Bottom edge
   j=1
   mean_G_cont_12_u2_nw(i,j) = p5*(G_cont_12_u2(i,j)+G_cont_12_u1(i,j))
   mean_G_cont_12_u2_ne(i,j) = p5*(G_cont_12_u2(i,j)+G_cont_12_u1(i+1,j))
   mean_G_cont_12_u2_sw(i,j) = p5*(G_cont_12_u2(i,jm+1)+G_cont_12_u1(i,jm))
   mean_G_cont_12_u2_se(i,j) = p5*(G_cont_12_u2(i,jm+1)+G_cont_12_u1(i+1,jm))
   ! Top edge
   mean_G_cont_12_u2_nw(i,jm+1) = mean_G_cont_12_u2_nw(i,j)
   mean_G_cont_12_u2_ne(i,jm+1) = mean_G_cont_12_u2_ne(i,j)
   mean_G_cont_12_u2_sw(i,jm+1) = mean_G_cont_12_u2_sw(i,j)
   mean_G_cont_12_u2_se(i,jm+1) = mean_G_cont_12_u2_se(i,j)
end do



if ( init_from_ic_prog_file ) then

   print *
   print *, "Initializing prognostic variables from ic_prog file"
   print *

   open (unit = 20, file = "./data/ic_prog", action = "read", form = "unformatted")


   ! Read in initial value of time (tau)
   read (20) tau_initial


   ! Read initial conditions prognostic variables
   ! do j = 1,jm
   !    read (20) (u(i,j,n4),  i = 1,im)
   !    read (20) (v(i,j,n4),  i = 1,im)
   !    read (20) (h_star(i,j,n4),  i = 1,im)
   ! end do

   close (20)

else    !  (.not.init_from_ic_prog_file)

   print *
   print *, "Initializing prognostic variables within model"
   print *

! Mimicking Williamson et al. (1992) test case 5 -- Zonal Flow over an Isolated Mtn.

   tau_initial = c0


   if ( water_mtn ) then   ! case with no topography, but perturbed fluid depth

      ! Initialize Cartesian velocity components
      u_init_u1(:,:) = c0 ! 30._dbl_kind*                                      &
                          !  Sin(inv_a_rad*5*(y_u(:,:)+30220._dbl_kind))*  &
                          !  Cos(inv_a_rad*2*(x_u(:,:)-10210._dbl_kind))  ! c0
      u_init_u2(:,:) = c0 ! 30._dbl_kind*                                      &
                          !  Sin(inv_a_rad*5*(y_v(:,:)+30220._dbl_kind))*  &
                          !  Cos(inv_a_rad*2*(x_v(:,:)-10210._dbl_kind))  ! c0
      v_init_u1(:,:) = c0 ! 10._dbl_kind*                                      &
                          !  Sin(inv_a_rad*7*(x_u(:,:)-55002._dbl_kind))*  &
                          !  Sin(inv_a_rad*y_u(:,:))  ! c0
      v_init_u2(:,:) = c0 ! 10._dbl_kind*                                      &
                          !  Sin(inv_a_rad*7*(x_v(:,:)-55002._dbl_kind))*  &
                          !  Sin(inv_a_rad*y_v(:,:))  ! c0

      ! Initialize prognostic covariant velocity components
      u1_cov(:,:,n4) = dx_dx1_u1(:,:)*u_init_u1(:,:) +               &
                       dy_dx1_u1(:,:)*v_init_u1(:,:)
      u2_cov(:,:,n4) = dx_dx2_u2(:,:)*u_init_u2(:,:) +               &
                       dy_dx2_u2(:,:)*v_init_u2(:,:)

      ! Temporary lines
      u1_cov(im+1,:,n4) = p5*(u1_cov(im+1,:,n4)+u1_cov(1,:,n4))
      u1_cov(1,:,n4) = u1_cov(im+1,:,n4)
      u2_cov(:,jm+1,n4) = p5*(u2_cov(:,jm+1,n4)+u2_cov(:,1,n4))
      u2_cov(:,1,n4) = u2_cov(:,jm+1,n4)

      sqrt_G_h_star(:,:,n4) = 8420._dbl_kind*sqrt_G_h(:,:)
      do j = 1, jm
         do i = 1, im
            r_test = ( (x_h(i,j)-x_ctr)**2 + (y_h(i,j)-y_ctr)**2 ) ** 0.5
            if ( r_test.le.r0 ) then
               sqrt_G_h_star(i,j,n4) = sqrt_G_h_star(i,j,n4) +             &
                              water_mtn_ht*(c1-r_test/r0)*sqrt_G_h(i,j)
            end if
         end do
      end do

   else   !  (.not.water_mtn), i.e., topographical mountain instead

      ! Initialize Cartesian velocity components
      u_init_u1(:,:) = 20._dbl_kind*Cos(inv_a_rad*y_u(:,:))
      u_init_u2(:,:) = 20._dbl_kind*Cos(inv_a_rad*y_v(:,:))
      v_init_u1(:,:) = c0
      v_init_u2(:,:) = c0

      ! Initialize prognostic covariant velocity components
      u1_cov(:,:,n4) = dx_dx1_u1(:,:)*u_init_u1(:,:) +               &
                       dy_dx1_u1(:,:)*v_init_u1(:,:)
      u2_cov(:,:,n4) = dx_dx2_u2(:,:)*u_init_u2(:,:) +               &
                       dy_dx2_u2(:,:)*v_init_u2(:,:)

      ! Temporary lines
      u1_cov(im+1,:,n4) = p5*(u1_cov(im+1,:,n4)+u1_cov(1,:,n4))
      u1_cov(1,:,n4) = u1_cov(im+1,:,n4)
      u2_cov(:,jm+1,n4) = p5*(u2_cov(:,jm+1,n4)+u2_cov(:,1,n4))
      u2_cov(:,1,n4) = u2_cov(:,jm+1,n4)

      h_star(:,:) = 5960._dbl_kind    ! initially h_0

      ! Put mountain in place

      if ( conical_mtn ) then

         do j = 1, jm
            do i = 1, im
      
               r_test = ( (x_h(i,j)-x_ctr)**2 + (y_h(i,j)-y_ctr)**2 ) ** 0.5
               if ( r_test.le.r0 ) then
                  hs(i,j) = mtn_ht*(c1-r_test/r0)
                  h_star(i,j) = h_star(i,j) - hs(i,j)
               else
                  hs(i,j) = c0
               end if
         
               ! adjust height for geostrophic balance
               h_star(i,j) = h_star(i,j) -                                    &
                  invgrav*Omega*a_rad*20._dbl_kind*Sin(inv_a_rad*y_h(i,j))**2

            end do
         end do

      else   ! Gaussian mountain case

         do j = 1, jm
            do i = 1, im

               r_test = ( (x_h(i,j)-x_ctr)**2 + (y_h(i,j)-y_ctr)**2 ) ** 0.5
               if ( r_test.le.(0.375_dbl_kind*pi*a_rad) ) then
                  hs(i,j) = mtn_ht*exp( - (r_test/r0_G)**2 )
                  h_star(i,j) = h_star(i,j) - hs(i,j)
               else
                  hs(i,j) = c0
               end if
         
               ! adjust height for geostrophic balance
               h_star(i,j) = h_star(i,j) -                                    &
                  invgrav*Omega*a_rad*20._dbl_kind*Sin(inv_a_rad*y_h(i,j))**2

            end do
         end do

      end if    ! if ( conical_mtn )

      sqrt_G_h_star(:,:,n4) = h_star(:,:)*sqrt_G_h(:,:)

   end if    ! if ( water_mtn )


end if    ! if ( init_from_ic_prog_file )



! Initialize Coriolis parameter
! f_cor(:,:) = 3.24e-04_dbl_kind     ! f-plane
f_cor(:,:) = c2*Omega*Sin(inv_a_rad*y_q(:,:))


end subroutine initialize_model



subroutine coordinate_solve (x1_in, x2_in, x_out, y_out)

! Subroutine to iteratively solve for Cartesian-coordinate values
! of model grid points (x1_in,x2_in)

real (kind=dbl_kind), intent(in) ::        &
          x1_in, x2_in     ! non-orthogonal coordinate grid-point values in
          
real (kind=dbl_kind), intent(out) ::       &
          x_out, y_out     ! Cartesian coordinate grid-point values out

real (kind=dbl_kind) ::                    &
          error_x1, error_x2, inv_Jac_11, inv_Jac_12, inv_Jac_21, inv_Jac_22
          
real (kind=dbl_kind) ::                    &
          temp, temp1, temp2, temp3, temp4, x1_test, x2_test

real (kind=dbl_kind), parameter ::         &
          conv = 1.E-08_dbl_kind   ! convergence criterion

integer :: iter

integer, parameter :: max_iter = 100


iter = 0
error_x1 = 1.E+10_dbl_kind
error_x2 = 1.E+10_dbl_kind

! First guess for x_out and y_out
x_out = x1_in
y_out = x2_in

do while ( (abs(error_x1).gt.conv).or.(abs(error_x2).gt.conv) )

   temp1 = c1 + C_fac*Cos(p5*(y_out*inv_a_rad-bet))*Cos(x_out*inv_a_rad-alf)
   temp2 = -p5*C_fac*Sin(x_out*inv_a_rad-alf)*Sin(p5*(y_out*inv_a_rad-bet))
   temp3 = -p5*D_fac*Sin(y_out*inv_a_rad-bet)*Sin(p5*(x_out*inv_a_rad-alf))
   temp4 = c1 + D_fac*Cos(p5*(x_out*inv_a_rad-alf))*Cos(y_out*inv_a_rad-bet)
   temp = c1 / ( temp1*temp4 - temp2*temp3 )
   inv_Jac_11 = temp * ( c1 + D_fac*Cos(p5*(x_out*inv_a_rad-alf))*   &
                                    Cos(y_out*inv_a_rad-bet) )
   inv_Jac_12 = temp * p5*C_fac*Sin(x_out*inv_a_rad-alf)*            &
                                Sin(p5*(y_out*inv_a_rad-bet))
   inv_Jac_21 = temp * p5*D_fac*Sin(y_out*inv_a_rad-bet)*            &
                                Sin(p5*(x_out*inv_a_rad-alf))
   inv_Jac_22 = temp * ( c1 + C_fac*Cos(p5*(y_out*inv_a_rad-bet))*   &
                              Cos(x_out*inv_a_rad-alf) )

   x1_test = x_out + a_rad*C_fac*Cos(p5*(y_out*inv_a_rad-bet))*      &
                                 Sin(x_out*inv_a_rad - alf)
   x2_test = y_out + a_rad*D_fac*Cos(p5*(x_out*inv_a_rad-alf))*      &
                                 Sin(y_out*inv_a_rad - bet)

   error_x1 = x1_test - x1_in
   error_x2 = x2_test - x2_in

   x_out = x_out - inv_Jac_11*error_x1 - inv_Jac_12*error_x2
   y_out = y_out - inv_Jac_21*error_x1 - inv_Jac_22*error_x2

   iter = iter + 1
   
   if ( iter.gt.max_iter ) then
      print *
      print *, "Max iterations reached in coordinate_solve"
      print *, "iter = ", iter
      print *, "   Program stopping"
      print *
      stop
   end if

end do

end subroutine coordinate_solve


!===================================================================
! END OF INITIALIZE_MODEL
!===================================================================


end module initialize
