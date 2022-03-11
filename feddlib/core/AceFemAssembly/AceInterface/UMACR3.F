c$Id: umacr3.f,v 1.1 1997/07/29 19:15:14 rlt Exp $
      subroutine umacr3(lct,ctl,prt)
 
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

c....  2000 J. Korelc   check error status automatically set by SMS elements
c....
c....  chkc	  report arror status to screen and output file and clear all error flags
c....  chkc,clea  clear all error flags and write report to output file
c....  chkc,,tag  
c....  chkc,clea,tag   tag is an arbitrary number included in a report that can be used to
c....                  locate the error
c
c.... error status:
c....   0 - no special events were detected during the session
c....   1 - warnings were detected during the session, (evaluation is still performed in
c....        a regular way, time step cutting is recommended but not necessary)
c....   2 - fatal errors were detected during the session (time step cutting is necessary)
c....   3 - fatal error (terminate the process)
c

      implicit   none
 
      logical    pcomp,prt
      character  lct*15
      real*8     ctl(*)

      include   'umac1.h'
      include   'sms.h'
	character  mes(5),st(4)
	integer imes(5),i,j
      data mes/"Missing user subroutine.",
     &    "Divergence in sub-iterative process.",
     &     "Non-physical element state.",
     &     "Non-physical element shape.",
     &     "Non-physical material point state."/
      data imes/ID_MissingSubroutine,ID_SubDivergence,ID_ElementState
     &    ,ID_ElementShape,ID_MaterialState/
      data st/"No errors.","Warning.","Warning -step cut."
     &    ,"Fatal error."/

      save
 
c     Set command word
 
      if(pcomp(uct,'mac3',4)) then
        uct='chkc'
        return
      endif
      j=idata(ID_ErrorStatus)
      Do i=1,5
        if(idata(imes(i)).ne.0) then
          write(iow,100) ctl(1),st(j+1),mes(i),idata(imes(i))
          if(lct(1:4).ne.'clea')
     &      write(*,100) ctl(1),st(j+1),mes(i),idata(imes(i))
          idata(imes(i))=0;
        endif
      End Do
      idata(ID_ErrorStatus)=ERROR_NoError
100   format(" REPORT ERROR (Tag=",f8.2,") ",a,a," No. of events= ",i4)
      end
