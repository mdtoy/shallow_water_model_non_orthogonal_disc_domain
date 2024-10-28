program calc_error_norms

!--------------------------------------------------------------------
! This program calculates L1, L2 and Linf norms of input variable.
! Data is read from netCDF files.
!--------------------------------------------------------------------



implicit none

include 'netcdf.inc'




integer, parameter :: dbl_kind = selected_real_kind(13)



integer :: i
integer :: im,jm

integer :: ncid_in,err
integer :: varid
integer :: dimid


real (kind = dbl_kind), allocatable ::                               &
           h_star(:), h_star_conv(:)

real (kind = dbl_kind) ::                                            &
           L1, L2, Linf, error, L1_denom, L2_denom, Linf_denom



! read in array limits
err = NF_OPEN("h_stars_for_error_norm_calc.TC_5_Gauss_mtn.dx_50km.nc",NF_NOWRITE,ncid_in)
err = NF_INQ_DIMID(ncid_in,'x_h',dimid)
err = NF_INQ_DIMLEN(ncid_in,dimid,im)
err = NF_INQ_DIMID(ncid_in,'y_h',dimid)
err = NF_INQ_DIMLEN(ncid_in,dimid,jm)




print *
print *, "   im =", im, "   jm =", jm
print *



! read in variables
allocate (h_star(im*jm))
allocate (h_star_conv(im*jm))
err = NF_INQ_VARID(ncid_in,'h_star_50km',varid)
err = NF_GET_VAR_DOUBLE(ncid_in,varid,h_star)
err = NF_INQ_VARID(ncid_in,'h_star_conv_interp',varid)
err = NF_GET_VAR_DOUBLE(ncid_in,varid,h_star_conv)





! Calculate error norms (Williamson et al. (1992))
L1   = 0._dbl_kind
L2   = 0._dbl_kind
Linf = 0._dbl_kind
L1_denom  = 0._dbl_kind
L2_denom  = 0._dbl_kind
Linf_denom = 0._dbl_kind
do i =  1, im*jm
   ! Flag masked values
   if ( h_star_conv(i) .gt. 1.E+30_dbl_kind ) then
      print *, "Masked value at i =", i, h_star_conv(i)
      cycle
   end if
   error = h_star(i) - h_star_conv(i)
   L1 = L1 + abs(error)
   L2 = L2 + error**2
   Linf = max(Linf,abs(error))
   L1_denom = L1_denom + abs(h_star_conv(i))
   L2_denom = L2_denom + h_star_conv(i)**2
   Linf_denom = max(Linf_denom,abs(h_star_conv(i)))
end do
L1 = L1 / L1_denom
L2 = L2**0.5 / L2_denom**0.5
Linf = Linf / Linf_denom


print *
print *, "L1 =", L1
print *, "L2 =", L2
print *, "Linf =", Linf
print *


! Temporary check
! print*, "Linf_denom =", Linf_denom
! do i = 1, im*jm
!    error = h_star(i) - h_star_conv(i)
!    if ( (abs(error)/Linf_denom).ge.7.E-04_dbl_kind ) then
!       print *, "i =", i, "     i/800 =", i/800, "    remainder =",   &
!            nint(800*(i/800._dbl_kind - i/800)),  &
!             "     abs(error) =",  &
!            abs(error), "     abs(error)/Linf_denom =", abs(error)/Linf_denom
!    end if
! end do



err = NF_CLOSE(ncid_in)


end program calc_error_norms
