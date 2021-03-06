c
c copyright notice
c                      Copyright 2013 Sandia Corporation             
c                                                                    
c     Under the terms of Contract DE-AC04-94AL85000 with Sandia      
c     Corporation, the U.S. Government retains certain rights in this
c     software.                                                      
c                                                                    
c     Redistribution and use in source and binary forms, with or     
c     without modification, are permitted provided that the following
c     conditions are met:                                            
c                                                                    
c      1. Redistributions of source code must retain the above       
c         copyright notice, this list of conditions and the following
c         disclaimer.                                                
c                                                                    
c      2. Redistributions in binary form must reproduce the above    
c         copyright notice, this list of conditions and the following
c         disclaimer in the documentation and/or other materials     
c         provided with the distribution.                            
c                                                                    
c      3. Neither the name of the Corporation nor the names of the   
c         contributors may be used to endorse or promote products    
c         derived from this software without specific prior written  
c         permission.                                                
c                                                                    
c     THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
c     EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,  
c     THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A    
c     PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA    
c     CORPORATION OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT,      
c     INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL     
c     DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF         
c     SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;   
c     OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  
c     LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      
c     (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF  
c     THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
c     SUCH DAMAGE.                                                   
c     
c
c************************************************************
c Subroutines in this file:
c  1. fmj
c  2. fmk 
c
c************************************************************
c
c
c****************************************************************
      subroutine fmj(jmax,kmax,x,y,q,rtxy,dj,j,jinc,kbeg,kend,
     &               x0,y0,fxp,fyp,fxs,fys,fmom)
c****************************************************************
c  Computes force and moment components due to pressure and
c  computes moment about the point xmom, ymom.
c
c  Assumes no-slip wall at j, fluid side of the wall at j+jinc.
c
c  Periodic case is handled automatically if kbeg,kend = 1,kmax
c----------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      __REAL x(jmax,kmax), y(jmax,kmax), q(jmax,kmax,3),
     &  rtxy(jmax,kmax,2,3), dj(jmax,kmax)
c
c      __REAL t11(jkmax), t22(jkmax), t12(jkmax)
      __REAL, allocatable:: t11(:),t22(:),t12(:)
 
      allocate (t11(kend),t22(kend),t12(kend))

c-----
c  Integrate pressure forces
c-----
      sgn = float( jinc )
      do 10 k=kbeg,kend-1
c
         pave = 0.5*( q(j,k+1,1) + q(j,k,1) )
c
         da1 = 0.5*( rtxy(j,k+1,1,2) + rtxy(j,k,1,2) )
         da2 = 0.5*( rtxy(j,k+1,1,3) + rtxy(j,k,1,3) )
c
         xave = 0.5*( x(j,k+1) + x(j,k) ) - x0
         yave = 0.5*( y(j,k+1) + y(j,k) ) - y0
c
c         fx = sgn*sw*pave*da1
c         fy = sgn*sw*pave*da2
         fx = sgn*pave*da1
         fy = sgn*pave*da2
c
         fxp = fxp - fx
         fyp = fyp - fy
         fmom = fmom - ( xave*fy - yave*fx )
10    continue
c-----
c  Inviscid
c-----
      if(ivis .eq. 0) then
         fxs = 0.0
         fys = 0.0
         return
      endif
c-----
c  Compute shear stress tensor
c-----
      do 20 k=kbeg,kend
         uxi = 0.5*sgn*( -3.*q(j,k,2) + 4.*q(j+jinc,k,2)
     &                                   - q(j+jinc+jinc,k,2) )
         vxi = 0.5*sgn*( -3.*q(j,k,3) + 4.*q(j+jinc,k,3)
     &                                   - q(j+jinc+jinc,k,3) )
c
         t11(k) = 2.*vnu*rtxy(j,k,1,2)*dj(j,k)*uxi
         t22(k) = 2.*vnu*rtxy(j,k,1,3)*dj(j,k)*vxi
         t12(k) = vnu*( rtxy(j,k,1,3)*uxi + rtxy(j,k,1,2)*vxi )*dj(j,k)
20    continue
c-----
c  Integrate skin-friction forces
c-----
      do 30 k=kbeg,kend-1
c
         t11ave = 0.5*( t11(k+1) + t11(k) )
         t22ave = 0.5*( t22(k+1) + t22(k) )
         t12ave = 0.5*( t12(k+1) + t12(k) )
         t21ave = t12ave
c
         da1 = 0.5*( rtxy(j,k+1,1,2) + rtxy(j,k,1,2) )
         da2 = 0.5*( rtxy(j,k+1,1,3) + rtxy(j,k,1,3) )
c
         xave = 0.5*( x(j,k+1) + x(j,k) ) - x0
         yave = 0.5*( y(j,k+1) + y(j,k) ) - y0
c
         fx = sgn*( t11ave*da1 + t12ave*da2 )
         fy = sgn*( t21ave*da1 + t22ave*da2 )
c
         fxs = fxs + fx
         fys = fys + fy
         fmom = fmom + xave*fy - yave*fx
30    continue

      deallocate (t11,t22,t12)

c-----
c  End of fmj
c-----
      return
      end
c
c
c****************************************************************
      subroutine fmk(jmax,kmax,x,y,q,rtxy,dj,k,kinc,jbeg,jend,
     &               x0,y0,fxp,fyp,fxs,fys,fmom)
c****************************************************************
c  Computes force and moment components due to pressure and
c  computes moment about the point xmom, ymom.
c
c  Assumes no-slip wall at k, fluid side of the wall at k+kinc.
c----------------------------------------------------------------
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      __REAL x(jmax,kmax), y(jmax,kmax), q(jmax,kmax,3),
     &  rtxy(jmax,kmax,2,3), dj(jmax,kmax)
c
c      __REAL t11(jkmax), t22(jkmax), t12(jkmax)
      __REAL, allocatable:: t11(:),t22(:),t12(:)

      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

      allocate (t11(jend),t22(jend),t12(jend))

c-----
c  Integrate pressure forces
c-----

      sgn = float( kinc )
c zero pressure forces
      fxp = 0.
      fyp = 0.
c zero viscous forces
      fxs = 0.
      fys = 0.

      if (DEBUG) then
       print*,''
       print*,' >> In fmk: nodeid,jbeg,jend = ',nodeid,jbeg,jend
       print*,''
      endif

      do 10 j=jbeg,jend-1
c
       pave = 0.5*( q(j+1,k,1) + q(j,k,1) )
c
       da1 = 0.5*( rtxy(j+1,k,2,2) + rtxy(j,k,2,2) )
       da2 = 0.5*( rtxy(j+1,k,2,3) + rtxy(j,k,2,3) )
c
       xave = 0.5*( x(j+1,k) + x(j,k) ) - x0
       yave = 0.5*( y(j+1,k) + y(j,k) ) - y0
c
       fx = sgn*pave*da1

       fy = sgn*pave*da2
c
       fxp = fxp - fx
       fyp = fyp - fy
       fmom = fmom - ( xave*fy - yave*fx )
10    continue

      if (DEBUG) then
       print*, 'In fmk: j,k,nodeid,q(j+1,k,1) = ',j,k,nodeid,q(jbeg,1,1)
       print*, 'In fmk: nodeid,rtxy(j+1,k,2,2) = ',
     &  nodeid,rtxy(jbeg+1,1,2,2)
       print*, 'In fmk: nodeid,x(j+1,k) = ',nodeid,x(jbeg+1,1)
       print*, 'In fmk: nodeid,y(j+1,k) = ',nodeid,y(jbeg+1,1)
       print*, 'In fmk: nodeid,fxp,fyp = ',nodeid,fxp,fyp
       call flush(6)
c       stop ' stop: in fm.f/fmk'
      endif

c-----
c  Inviscid
c-----
      if(ivis .eq. 0) then
         fxs = 0.0
         fys = 0.0
         return
      endif
c-----
c  Compute shear stress tensor
c-----
      do 20 j=jbeg,jend
         ueta = 0.5*sgn*( -3.*q(j,k,2) + 4.*q(j,k+kinc,2)
     &                                   - q(j,k+kinc+kinc,2) )
         veta = 0.5*sgn*( -3.*q(j,k,3) + 4.*q(j,k+kinc,3)
     &                                   - q(j,k+kinc+kinc,3) )
c
         t11(j) = 2.*vnu*rtxy(j,k,2,2)*dj(j,k)*ueta
         t22(j) = 2.*vnu*rtxy(j,k,2,3)*dj(j,k)*veta
         t12(j) = vnu*( rtxy(j,k,2,3)*ueta+rtxy(j,k,2,2)*veta )*dj(j,k)
20    continue
c-----
c  Integrate skin-friction forces
c-----
      do 30 j=jbeg,jend-1
c
         t11ave = 0.5*( t11(j+1) + t11(j) )
         t22ave = 0.5*( t22(j+1) + t22(j) )
         t12ave = 0.5*( t12(j+1) + t12(j) )
         t21ave = t12ave
c
         da1 = 0.5*( rtxy(j+1,k,2,2) + rtxy(j,k,2,2) )
         da2 = 0.5*( rtxy(j+1,k,2,3) + rtxy(j,k,2,3) )
c
         xave = 0.5*( x(j+1,k) + x(j,k) ) - x0
         yave = 0.5*( y(j+1,k) + y(j,k) ) - y0
c
         fx = sgn*( t11ave*da1 + t12ave*da2 )
         fy = sgn*( t21ave*da1 + t22ave*da2 )
c
         fxs = fxs + fx
         fys = fys + fy
         fmom = fmom + xave*fy - yave*fx
30    continue

      deallocate (t11,t22,t12)
c-----
c  End of fmk
c-----
      return
      end
