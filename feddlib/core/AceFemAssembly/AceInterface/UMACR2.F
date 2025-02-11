c$Id: umacr2.f,v 1.1 1997/07/29 19:15:14 rlt Exp $
      subroutine umacr2(lct,ctl,prt)
 
c      * * F E A P * * A Finite Element Analysis Program
 
c....  Copyright (c) 1984-1997: Robert L. Taylor
 
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose:  User interface for adding solution command language
c                instructions.
 
c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters
c         prt       - Flag, output if true
 
c      Outputs:
c         N.B.  Users are responsible for command actions.  See
c               programmers manual for example.
c-----[--+---------+---------+---------+---------+---------+---------+-]
 
      implicit   none
 
      logical    pcomp,prt
      character  lct*15
      real*8     ctl(*)
 
      include   'iofile.h'
      include   'umac1.h'
      include   'tdata.h'
      include   'ldata.h'
      include   'auto2.h'

      integer  n0, n1, n2, n3,lloadnbpo,lloadprs,llodout
      real *8 llodmax,llodtol
      common/llod/n0,n1,n2,n3,lloadnbpo,lloadprs,llodmax,llodtol
      integer loopterm
      common/cloopterm/loopterm
	real*8 a,b,c
      save

      if(pcomp(uct,'mac2',4)) then
        uct='llod'
	  return
      endif
 
	if(lct(1:3).eq.'set') then
          if(ctl(1).eq.0.d0) then
              llodmax=1d0
          else
              llodmax=ctl(1)
          endif
          if(ctl(2).eq.0.d0) then
              llodtol=.01d0
          else
              llodtol=ctl(2)
          endif               
          if(int(ctl(3)).eq.0) then
              lloadnbpo=1
          else
              lloadnbpo=int(ctl(3))
          endif
          lloadprs=0           
          loopterm=0
          write(iow,*)'Set limit load tol=',llodtol,
     -         " bifurcation=",lloadnbpo
	else
          if(n1.lt.lloadnbpo) then
             if(ttim+dt.gt.llodmax)then
                dt=(llodmax-ttim)/2.d0
                lloadprs=1
             else
                if(lloadprs.ne.0) dt=dt/2.d0
             endif
           else
             lloadprs=1
             call autbac(0.d0)
	       dt=dt/2.d0
	     endif
           write(iow,*)' CURRENT_LOAD= ',ttim+dt
           write(*,*)' CURRENT_LOAD= ',ttim+dt
	     if(abs(dt/(ttim+dt)).le.llodtol) then
             write(iow,*)' LIMIT_LOAD= ',ttim+dt
             write(*,*)' LIMIT_LOAD= ',ttim+dt
             loopterm=1
           endif
      endif
      end