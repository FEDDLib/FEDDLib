c$Id: umacr0.f,v 1.1 1997/07/29 19:15:14 rlt Exp $
      subroutine umacr0(lct,ctl,prnt)
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
	include 'eqsym.h'
      logical    pcomp,prnt,setvar,ualloc,DEBUG
      parameter (DEBUG=.true.)      
      character  lct*15
      real*8     ctl(3),rnorm
      integer    i,j,k,l,ts,ni,ns(100)
      save
c      write(*,'(a5,a5,3f10.5)')' sens ',lct,ctl
c     Set command word
      if(pcomp(uct,'mac0',4)) then
        uct='sens'
	  return
      endif
      if(lct(1:3).eq.'set') then
          upsens=2
          write(iow,*)'S e n s i t i v i t y   set space'
          if(ior.lt.0 .and. prnt)then
              write(*,*)'S e n s i t i v i t y   set space'
	    endif
          if(int(ctl(1)).eq.0) then
              setvar=ualloc(upsens,'SENS',nneq*nsenpa,ipr)
              setvar=ualloc(upsens+1,'SENS0',nneq*nsenpa,ipr)
              nsenalo=nsenpa
          else
              setvar=ualloc(upsens,'SENS',nneq*int(ctl(1)),ipr)
              setvar=ualloc(upsens+1,'SENS0',nneq*int(ctl(1)),ipr)
              nsenalo=int(ctl(1))
          endif
          senstype=1
          sensti=1
	elseif(lct(1:4).eq.'type') then
          if(int(ctl(1)).eq.0) then
              senstype=1
          else
              senstype=int(ctl(1))
          endif
	elseif(lct(1:4).eq.'para') then
          if(int(ctl(1)).eq.0) then
              sensti=1
          else
              sensti=int(ctl(1))
          endif
      elseif(lct(1:4).eq.'solv') then
          if(int(ctl(1)).eq.0) then
            ni=nsenalo
            do i=1,ni
                ns(i)=i
            enddo
          elseif(int(ctl(2)).eq.0) then
            ni=1
            ns(1)=int(ctl(1))
          else
            ni=int(ctl(2))-int(ctl(1))+1
            do i=0,ni-1
                ns(i+1)=int(ctl(1))+i
            enddo
          endif   
          do i=1,ni
              tsenp=ns(i)              
              write(iow,*)'S e n s i t i v i t y   external : ',tsenp
              if(ior.lt.0 .and. prnt)then
                  write(*,*)'S e n s i t i v i t y   external : ',tsenp
	        endif
              ts=up(upsens)+(tsenp-1)*nneq
              call pmove(hr(ts),hr(up(upsens+1)+(tsenp-1)*nneq),nneq)
	        call pzero(hr(np(26)),neq)
              hflgu=.true.
              call formfe(hr(np(40)),hr(np(26)),hr(np(26)),hr(np(26)),
     -      		  .false.,.true.,.false.,.false.,20,1,numel,1)
              call dasol (hr(nal),hr(nau),hr(na),
     -              hr(np(26)),mr(np(20+npart)),neqs,neq,rnorm)
              call pmovec(mr(np(31)),hr(np(26)),hr(ts),nneq)
              write(iow,*)'S e n s i t i v i t y   Internal : ',tsenp
              if(ior.lt.0 .and. prnt)then
                  write(*,*)'S e n s i t i v i t y   Internal : ',tsenp
	        endif
              hflgu=.true.
              call formfe(hr(np(40)),hr(np(26)),hr(np(26)),hr(np(26)),
     -      		  .false.,.false.,.false.,.false.,21,1,numel,1)
          enddo
      elseif(lct(1:4).eq.'inte') then

      elseif(lct(1:4).eq.'disp') then
          if(int(ctl(2)).eq.0) then
            ni=min(100,numnp)
            do i=1,ni
                ns(i)=i
            enddo
          elseif(int(ctl(3)).eq.0) then
            ni=1
            ns(1)=int(ctl(2))
          else
            ni=int(ctl(3))-int(ctl(2))+1
            do i=0,ni-1
                ns(i)=int(ctl(1))+i
            enddo
          endif   
          if(int(ctl(1)).eq.0) then
             do tsenp=1,nsenalo
              ts=up(upsens)+(tsenp-1)*nneq
              write(iow,*)'S e n s i t i v i t y   Display : ',tsenp
              if(ior.lt.0 .and. prnt)then
                  write(*,*)'S e n s i t i v i t y   Display : ',tsenp
	        endif
              call prtlis(hr(np(43)),hr(ts),ttim,prop,ndm,ndf,
     -              ni,ns,1,.false.)
             enddo
          else
             tsenp=int(ctl(1))
             ts=up(upsens)+(tsenp-1)*nneq
              write(iow,*)'S e n s i t i v i t y   Display : ',tsenp
              if(ior.lt.0 .and. prnt)then
                  write(*,*)'S e n s i t i v i t y   Display : ',tsenp
	        endif
             call prtlis(hr(np(43)),hr(ts),ttim,prop,ndm,ndf,
     -              ni,ns,1,.false.)
         endif
      else
          write(*,*) 'Unknown sensitivity command :',lct
	  endif
      end
      subroutine SMSMoveGlobal(ix,v,s,nen,ndf)
      implicit none
      integer nen,ix(nen),i,j,k,l,m,ndf
      real*8  v(*),s(*)
      l=1
      do i=1,nen
           j = abs(ix(i))
           if(j.gt.0) then
              m=ndf*(j-1)
              do k=1,ndf
                 s(l)=v(m+k)
                 l=l+1
              enddo
           endif
      enddo
      end