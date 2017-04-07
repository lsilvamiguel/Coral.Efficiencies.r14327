c        Subset of functions for shower separation in ECAL2.  A.Lednev 21.10.2008.
c     cell, dcelldx, dcelldy, d2celldx2, d2celldy2, d2celldxdy, IsItLG, errcell2.

c        A.Lednev        2008.
c       Nucl.Instr.& Meth. A 366 (1995), pp.292-297
c        Function to calculate energy deposition in cell.
c        x,y are the shower coordinates (in mm) from the cell center.
c        LG = 1 for Lead glass, LG.ne.1 for shashlyk.
      real function cell(x,y,LG)
      implicit none
      integer LG
      real x, y
#include "parameters.inc"
#include "showerprof.inc"
      integer i
      real cum, p, sum, u, v, x1, x2, y1, y2
      cum(u,v,p)=atan(u/p)+atan(v/p)+atan2(u*v,p*sqrt(u*u+v*v+p*p))
      if(x**2+y**2.gt.40000.) then
        cell=0.
        return
      endif
      x1=x+dcell/2.
      y1=y+dcell/2.
      x2=x1-dcell
      y2=y1-dcell
      sum=0.
      if(LG.eq.1) then
        do i=1,2                ! Lead glass
          sum=sum+a(i)*(cum(x1,y1,b(i))+cum(x2,y2,b(i))
     -                 -cum(x1,y2,b(i))-cum(x2,y1,b(i)))
        enddo
      else
        do i=3,4                ! Shashlyk
          sum=sum+a(i)*(cum(x1,y1,b(i))+cum(x2,y2,b(i))
     -                 -cum(x1,y2,b(i))-cum(x2,y1,b(i)))
        enddo
      endif
      cell=sum/pi2
      return
      end
c
      real function dcelldx(x,y,LG)
      implicit none
      integer LG
      real x, y
#include "parameters.inc"
#include "showerprof.inc"
      integer i
      real cum, p, sum, u, v, x1, x2, y1, y2
      cum(u,v,p)=p*(1.+v/sqrt(p*p+u*u+v*v))/(p*p+u*u)
      if(x**2+y**2.gt.40000.) then
        dcelldx=0.
        return
      endif
      x1=x+dcell/2.
      y1=y+dcell/2.
      x2=x1-dcell
      y2=y1-dcell
      sum=0.
      if(LG.eq.1) then
        do i=1,2
          sum=sum+a(i)*(cum(x1,y1,b(i))+cum(x2,y2,b(i))
     -                 -cum(x1,y2,b(i))-cum(x2,y1,b(i)))
        enddo
      else
        do i=3,4
          sum=sum+a(i)*(cum(x1,y1,b(i))+cum(x2,y2,b(i))
     -                 -cum(x1,y2,b(i))-cum(x2,y1,b(i)))
        enddo
      endif
      dcelldx=sum/pi2
      return
      end
c
      real function dcelldy(x,y,LG)
      implicit none
      integer LG
      real x, y
#include "parameters.inc"
#include "showerprof.inc"
      integer i
      real cum, p, sum, u, v, x1, x2, y1, y2
      cum(u,v,p)=p*(1.+u/sqrt(p*p+u*u+v*v))/(p*p+v*v)
      if(x**2+y**2.gt.40000.) then
        dcelldy=0.
        return
      endif
      x1=x+dcell/2.
      y1=y+dcell/2.
      x2=x1-dcell
      y2=y1-dcell
      sum=0.
      if(LG.eq.1) then
        do i=1,2
          sum=sum+a(i)*(cum(x1,y1,b(i))+cum(x2,y2,b(i))
     -                 -cum(x1,y2,b(i))-cum(x2,y1,b(i)))
        enddo
      else
        do i=3,4
          sum=sum+a(i)*(cum(x1,y1,b(i))+cum(x2,y2,b(i))
     -                 -cum(x1,y2,b(i))-cum(x2,y1,b(i)))
        enddo
      endif
      dcelldy=sum/pi2
      return
      end
c
      real function d2celldx2(x,y,LG)
      implicit none
      integer LG
      real x, y
#include "parameters.inc"
#include "showerprof.inc"
      integer i
      real cum, p, sum, u, v, x1, x2, y1, y2
      cum(u,v,p)=-p*u/(p*p+u*u)**2*(2.+v*(2.+(p*p+u*u)/(p*p+u*u+v*v))
     /          /sqrt(p*p+u*u+v*v))
      if(x**2+y**2.gt.40000.) then
        d2celldx2=0.
        return
      endif
      x1=x+dcell/2.
      y1=y+dcell/2.
      x2=x1-dcell
      y2=y1-dcell
      sum=0.
      if(LG.eq.1) then
        do i=1,2
          sum=sum+a(i)*(cum(x1,y1,b(i))+cum(x2,y2,b(i))
     -                 -cum(x1,y2,b(i))-cum(x2,y1,b(i)))
        enddo
      else
        do i=3,4
          sum=sum+a(i)*(cum(x1,y1,b(i))+cum(x2,y2,b(i))
     -                 -cum(x1,y2,b(i))-cum(x2,y1,b(i)))
        enddo
      endif
      d2celldx2=sum/pi2
      return
      end
c
      real function d2celldy2(x,y,LG)
      implicit none
      integer LG
      real x, y
#include "parameters.inc"
#include "showerprof.inc"
      integer i
      real cum, p, sum, u, v, x1, x2, y1, y2
      cum(u,v,p)=-p*v/(p*p+v*v)**2*(2.+u*(2.+(p*p+v*v)/(p*p+u*u+v*v))
     /          /sqrt(p*p+u*u+v*v))
      if(x**2+y**2.gt.40000.) then
        d2celldy2=0.
        return
      endif
      x1=x+dcell/2.
      y1=y+dcell/2.
      x2=x1-dcell
      y2=y1-dcell
      sum=0.
      if(LG.eq.1) then
        do i=1,2
          sum=sum+a(i)*(cum(x1,y1,b(i))+cum(x2,y2,b(i))
     -                 -cum(x1,y2,b(i))-cum(x2,y1,b(i)))
        enddo
      else
        do i=3,4
          sum=sum+a(i)*(cum(x1,y1,b(i))+cum(x2,y2,b(i))
     -                 -cum(x1,y2,b(i))-cum(x2,y1,b(i)))
        enddo
      endif
      d2celldy2=sum/pi2
      return
      end
c
      real function d2celldxdy(x,y,LG)
      implicit none
      integer LG
      real x, y
#include "parameters.inc"
#include "showerprof.inc"
      integer i
      real cum, p, sum, u, v, x1, x2, y1, y2
      cum(u,v,p)=p/sqrt(p*p+u*u+v*v)**3
      if(x**2+y**2.gt.40000.) then
        d2celldxdy=0.
        return
      endif
      x1=x+dcell/2.
      y1=y+dcell/2.
      x2=x1-dcell
      y2=y1-dcell
      sum=0.
      if(LG.eq.1) then
        do i=1,2
          sum=sum+a(i)*(cum(x1,y1,b(i))+cum(x2,y2,b(i))
     -                 -cum(x1,y2,b(i))-cum(x2,y1,b(i)))
        enddo
      else
        do i=3,4
          sum=sum+a(i)*(cum(x1,y1,b(i))+cum(x2,y2,b(i))
     -                 -cum(x1,y2,b(i))-cum(x2,y1,b(i)))
        enddo
      endif
      d2celldxdy=sum/pi2
      return
      end
c
      integer function IsItLG(iaddr)
      implicit none
      integer iaddr,ix,iy
      IsItLG=0
      ix=iaddr/48
      iy=iaddr-ix*48
      if(ix.ge.20.and.ix.le.43.and.iy.ge.12.and.iy.le.35) return
      if(ix.ge.16.and.ix.le.47.and.iy.ge.16.and.iy.le.31) return
      if(ix.ge. 8.and.ix.le.55.and.iy.ge.19.and.iy.le.28) return
      IsItLG=1
      return
      end
c
      real function errcell2(e)
      implicit none
      real e
      real err_a, err_b, err_c
      data err_a,err_b,err_c/0.15,0.015,0.05/
      errcell2=err_a**2*e+(err_b*e)**2+err_c**2
      return
      end
