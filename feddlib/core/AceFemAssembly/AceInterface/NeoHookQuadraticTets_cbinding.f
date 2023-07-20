      module c_routines
        use, intrinsic :: iso_c_binding
        use f_routines
        implicit none
      contains
        subroutine skr1(v,d,ul,ul0,xl,s,p,ht,hp) bind(c,name='skr1')
            implicit NONE
            real(c_double), intent(inout) :: v(1002)
            real(c_double), intent(inout) :: d(2)
            real(c_double), intent(inout) :: ul(3,10)
            real(c_double), intent(inout) :: ul0(3,10)
            real(c_double), intent(inout) :: xl(3,10)
            real(c_double), intent(inout) :: s(30,30)
            real(c_double), intent(inout) :: p(30)
            real(c_double), intent(inout) :: ht(*)
            real(c_double), intent(inout) :: hp(*)
            call SKR10(v,d,ul,ul0,xl,s,p,ht,hp)
        end subroutine skr1

        subroutine spp1(v,d,ul,ul0,xl,s,p,ht,hp,sg,sg0,sxd,gpost,npost) 
     &  bind(c,name='spp1')
            implicit NONE
            real(c_double), intent(inout) :: v(1002)
            real(c_double), intent(inout) :: d(2)
            real(c_double), intent(inout) :: ul(3,10)
            real(c_double), intent(inout) :: ul0(3,10)
            real(c_double), intent(inout) :: xl(3,10)
            real(c_double), intent(inout) :: s(30,30)
            real(c_double), intent(inout) :: p(30)
            real(c_double), intent(inout) :: ht(*)
            real(c_double), intent(inout) :: hp(*)
            real(c_double), intent(inout) :: sg(*)
            real(c_double), intent(inout) :: sg0(*)
            real(c_double), intent(inout) :: sxd(30)
            real(c_double), intent(inout) :: gpost(64,21)
            real(c_double), intent(inout) :: npost(10,6)
            call SPP10(v,d,ul,ul0,xl,s,p,ht,hp,sg,sg0,sxd,gpost,npost)
        end subroutine spp1
      end module c_routines
