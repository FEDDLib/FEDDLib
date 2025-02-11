c$Id: umacr1.f,v 1.1 1997/07/29 19:15:14 rlt Exp $
      subroutine umacr1(lct,ctl,prnt)
 
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
 
      include 'sms.h'
      include 'ndata.h'
      include 'part0.h'
      include 'umac1.h'
      include 'prlod.h'
      include 'hdatam.h'

      logical    pcomp,prnt
      character  lct*15
      real*8     ctl
      integer i,j,k,ts
  
      save
 
c     Set command word
 
      if(pcomp(uct,'mac1',4)) then
          uct='tecp'
          return
      endif
 
 
	if(lct(1:3).eq.'set') then
          write(iow,*)'T e c p l o t   set'
          if(ior.lt.0 .and. prnt)then
              write(*,*)'T e c p l o t  set'
	    endif
          tecout=25
	    open(tecout,file="tecplot.dat")
	    write(tecout,*)"TITLE=""FEAP output"""
          tecnode=26
	    open(tecnode,access='direct',status='replace',
     *		file='tecnode.tmp',recl=2*8*ntecdata)
      else

	    do i=1,numnp
		    write(tecnode,rec=i)(0.d0,j=1,2*ntecdata)
	    enddo

          hflgu=.false.
          call formfe(hr(np(40)),hr(np(26)),hr(np(26)),hr(np(26)),
     - 		  .false.,.false.,.false.,.false.,22,1,numel,1)

          write(tecout,'(a,e10.3,a,i4,a,i4,a)')
     -      "ZONE T=""",ttim,""" N=",numnp,", E=",numel,
     -     	", F=FEPOINT, ET=QUADRILATERAL"

   	    do k=1,numnp
              do tsenp=1,nsenalo
                  ts=up(upsens)+(tsenp-1)*nneq
                  do i=0,ndf-1
                      tecx(ntecdata+ndf*(tsenp-1)+i)=
     -                    hr(np(ts)+ndf*(k-1)+i)
                  enddo
              enddo

		    read(tecnode,rec=k)(tecp(i),tecx(i),i=1,ntecdata)
		    write(tecout,'(100e12.5)')(tecx(i)/tecp(i),i=1,ntecdata),
     -            (tecx(i),i=ntecdata+1,ntecdata+nsenalo)
	    enddo


	    do k=1,numel
		    write(tecout,'(20i4)')(mr(np(33)+nen*(k-1)+i-1),i=1,nen)
	    enddo
          write(iow,*)'Tecplot end'
          if(ior.lt.0 .and. prnt)then
              write(*,*)'Tecplot end'
	    endif

      endif

      end
