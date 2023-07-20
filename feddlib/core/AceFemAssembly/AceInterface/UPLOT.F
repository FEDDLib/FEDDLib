c$Id: uplot.f,v 1.1 1997/07/29 19:15:14 rlt Exp $
      subroutine uplot(ct)

      implicit none
	include 'sms.h'
      include 'cdat1.h'
      include 'ndata.h'
      include 'part0.h'
      include 'umac1.h'
      include 'prlod.h'
      include 'pview.h'
      include 'pdata3.h'

      real*8  ct(3)
      integer k1,k2,k3,ts,nn,i,k4,j,nix,nk,nxd
      logical ll,palloc

      write(*,*) ' User plot '
c     [uplo,k1,k2] - k1 = sensitivity number (1-nsenpa)
c     [uplo,k1,k2] - k2 = dimension number (1-ndf)
c                    k3 = # lines (fill if <,= 0)
8     k2  = ct(3)
      nk   = abs(k2)
      nn   = max(1,nk)

      i   = max(1,abs(int(ct(1))))
      i   = min(i,nsenalo)

      j=max(1,min(ndf,int(ct(2))))

      ts=up(upsens)+(i-1)*nneq
      if(k2.le.0) then
        nn = -1
      endif
      k4=0
      call rprint(hr(ts+j-1),hr(np(43)),numnp,ndm,ndf,k4)

      call pdefm(hr(np(43)),hr(np(40)),cs,hr(np(45)),ndm,ndf,numnp,
     &           hr(np(53)))
      
	if(hide) then
	   nix=54
	   nxd=7
	else
	   nix=33
	   nxd=nen1
	endif  
      k4=1      
	call pltcon(hr(np(53)),mr(np(32)),mr(np(nix)),mr(np(62)),
     -               hr(ts+j-1),
     &              nie,3,ndf,nxd,k4,nn,i,2,.false.)

      end
