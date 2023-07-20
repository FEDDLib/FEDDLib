      module c_routines_DeformationDiffusionNeoHook
        use, intrinsic :: iso_c_binding
        use f_routines_DeformationDiffusionNeoHook
        implicit none
      contains
        subroutine skr_DDNH(v,d,ul,ul0,xl,s,p,ht,hp,deltat) 
     &  bind(c,name='skr_DDNH')
            implicit NONE
            real(c_double), intent(inout) :: v(2238)
            real(c_double), intent(inout) :: d(6) 
            real(c_double), intent(inout) :: ul(3,20)
            real(c_double), intent(inout) :: ul0(3,20)
            real(c_double), intent(inout) :: xl(3,20)
            real(c_double), intent(inout) :: s(60,60)
            real(c_double), intent(inout) :: p(60)
            real(c_double), intent(inout) :: ht(*)
            real(c_double), intent(inout) :: hp(*)
            real(c_double), intent(inout) :: deltat(1)
            call SKR10(v,d,ul,ul0,xl,s,p,ht,hp,deltat)
         end subroutine skr_DDNH
      end module c_routines_DeformationDiffusionNeoHook   
