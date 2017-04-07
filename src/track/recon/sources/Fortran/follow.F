      subroutine follow(pb,lpst,lpnd,xadd,yadd,zadd,dya1,dza1)
**************************************************************************************
*                                                                                    *
*---- recompute positions yadd,zadd and slopes dya1,dza1  for part of track in SM1...*
*                                                                                    *
**************************************************************************************
#include "PARAM.h"
#include "PLANES.h"

      dimension xadd(*),yadd(*),zadd(*),dya1(*),dza1(*)
      DOUBLE PRECISION  xx(3),h_field(3)
*------------------
***       print '(" follow..:",2i4)',lpst,lpnd

          do kv = lpst,lpnd
          xx(3) = xadd(kv)
          xx(1) = yadd(kv)
          xx(2) = zadd(kv)
         call mag_field(xx,h_field)
         fac = 0.3*sqrt(1.+dya1(kv)**2 + dza1(kv)**2)
         pys = fac*(h_field(3)*dza1(kv) + h_field(1) *
     1     dya1(kv) * dza1(kv) - h_field(2) *
     2     (1. + dya1(kv)**2))
         pzs = fac*(-h_field(3) * dya1(kv) - h_field(2) *
     1     dya1(kv) * dza1(kv) + h_field(1) *
     2     (1. + dza1(kv)**2))
         ys = pys * pb
         zs = pzs * pb
         dx = xadd(kv+1) - xadd(kv)
         dya1(kv+1) = dya1(kv) + ys * dx
         dza1(kv+1) = dza1(kv) + zs * dx
         yadd(kv+1) = yadd(kv) + dya1(kv) * dx + 0.5 * ys * dx**2
         zadd(kv+1) = zadd(kv) + dza1(kv) * dx + 0.5 * zs * dx**2
 
           hx = h_field(1)
           hy = h_field(2)
           hz = h_field(3)

***       print '("follow",i3,6f10.5)',kv,yadd(kv),zadd(kv),xadd(kv),
***  1      hx,hy,hz         
          enddo
      end
