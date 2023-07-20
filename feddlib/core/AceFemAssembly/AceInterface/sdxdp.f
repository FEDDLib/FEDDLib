! beam - dimensions (0,L)*(-1,1)
! height discretised by the quadratic functions
! 3 shape parameters are needed
! 1 material
      subroutine sensspec1(ma,tsenp,senstype,sensti)
      integer ma,tsenp,senstype,sensti,types(1,3)
      integer params(1,3)
      data types/2,2,2/
      data params/1,2,3/
      senstype=types(ma,tsenp)
      sensti=params(ma,tsenp)
      return
      end
      subroutine sdxdp1(ts,xl,sxd,p,ns,nd,nn)
      implicit logical (b)
      real*8 xl(nd,nn),sxd(nd,nn),p(ns),x,y,vx,vy
      real*8 l,yr
      integer ts,nd,nn,ns
      l=10.d0
      do i=1,nn
        x=xl(1,i)
        y=xl(2,i)
        yr=y/(1.d0+p(1)*x/l+p(2)*(1-x/l)+p(3)*x*(1-x/l))
        vx=0.d0
        if(ts.eq.1) then
          vy=yr*x/l
        elseif(ts.eq.2) then
          vy=yr*(1-x/l)
        elseif(ts.eq.3) then
          vy=yr*x*(1-x/l)
        else
          write(*,*)'dxdp - unknown parameter',ts
        endif
        sxd(1,i)=vx
        sxd(2,i)=vy
      enddo
      end



      subroutine sdxdp2(ts,xl,sxd,p,ns,nd,nn)
!2
! beam - dimensions (0,L)*(0,H)
! length discretised by the linear functions 1 parameter
      real*8 xl(nd,nn),sxd(nd,nn),p(ns),x,y,vx,vy
      real*8 l
      integer ts,nd,nn,i,ns
	l=10.d0
      do i=1,nn
        x=xl(1,i)
        y=xl(2,i)
        sxd(1,i)=x/l
        sxd(2,i)=0.d0
      enddo
      end 

      subroutine sdxdp3(ts,xl,sxd,p,ns,nd,nn)
!3
! beam - dimensions (0,L)*(0,H)
! length discretised by the exponential functions (x/l)**n
      real*8 xl(nd,nn),sxd(nd,nn),p(ns),x,y,vx,vy
      real*8 l
      integer ts,nd,nn,i,ns
	l=10.d0
      do i=1,nn
        x=xl(1,i)
        y=xl(2,i)
        sxd(1,i)=(x/l)**20
        sxd(2,i)=0.d0
      enddo
      end 


! beam - dimensions (0,L)*(0,H)
! length discretised by the linear functions 1 parameter
! 2 materials of equal length
! 3th material param in section 2
! 2th material param in section 1
      subroutine sensspec4(ma,tsenp,senstype,sensti)
      integer ma,tsenp,senstype,sensti
      integer params(3,2),types(3,2)
      data types/1,3,2, 3,1,2/
      data params/3,0,1, 0,1,1/
      senstype=types(tsenp,ma)
      sensti=params(tsenp,ma)
      return
      end
      subroutine sdxdp4(ts,xl,sxd,p,ns,nd,nn)
      real*8 xl(nd,nn),sxd(nd,nn),p(ns),x,y,vx,vy
      real*8 l
      integer ts,nd,nn,i,ns
      l=10.d0
      do i=1,nn
        x=xl(1,i)
        y=xl(2,i)
        sxd(1,i)=x/l
        sxd(2,i)=0.d0
      enddo
      end 

!generic test g2d4ni
      subroutine sensspec(ma,tsenp,senstype,sensti)
      integer ma,tsenp,senstype,sensti
      integer params(3,1),types(3,1)
      data types/1,1,2/
      data params/1,2,1/
      senstype=types(tsenp,ma)
      sensti=params(tsenp,ma)
      return
      end
      subroutine sdxdp(ts,xl,sxd,p,ns,nd,nn)
      real*8 xl(nd,nn),sxd(nd,nn),p(ns),x,y,vx,vy
      real*8 l
      integer ts,nd,nn,i,ns
      do i=1,nn
       do j=1,nd
        sxd(j,i)=0.d0
        sxd(j,i)=0.d0
       enddo
      enddo
      sxd(1,2)=1.0
      end 


