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
c  1. lrelax
c  2. ljsweep
c  3. lksweep
c
c************************************************************
c
c
c************************************************************************
      subroutine lrelax(jmax,kmax,q,dq,s,
     &                  btc,bjm,bjp,bkm,bkp)
c************************************************************************
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      __INTEGER stat(MPI_STATUS_SIZE)

      __REAL q(jmax,kmax,3),dq(jmax,kmax,3),s(jmax,kmax,3),
     &  btc(jmax,kmax,3,3),bjm(jmax,kmax,3,3),bjp(jmax,kmax,3,3),
     &  bkm(jmax,kmax,3,3), bkp(jmax,kmax,3,3)
      __INTEGER jmax, kmax
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

c-----
c  Zero out dq array 
c-----
      do 9 i=1,3
      do 9 k=1,kmax
      do 9 j=1,jmax
         dq(j,k,i) = 0.0
9     continue

      if (DEBUG.and.nodeid.eq.0) then
       print*,''
       print*, 'DWB10.5: BEFORE call to btlu3j'
       print*, 'DWB10.5: jmax,kmax,dq(2,1,1) = ',jmax,kmax,dq(2,1,1)
       print*, 'DWB10.5: btc(2,2,2,2) =',btc(2,2,2,2)
c       print*, 'DWB10.5: bjm(2,2,2,2), btc(2,2,2,2), bjp(2,2,2,2) = ',
c     &  bjm(2,2,2,2),btc(2,2,2,2),bjp(2,2,2,2)
       print*, 'DWB10.5: btc(2,2,1,1) = ',btc(2,2,1,1)
       print*,''
c       stop 'stop: in miniSMAC lrelax.f/lrelax 10'
      endif 

c-----
c  J-sweep: solve for lines of varying j and constant k.
c  First perform lu decomposition of each j-line for all grid
c  points in zones where j-sweeps will be done.
c-----
      njsmax = 0
c njsp = njswp = 2: set in initia.f
         njsmax = max( njsmax, njsp )
         if(njsp .gt. 0) then
           call btlu3j(jmax,kmax,bjm,btc,bjp,1,jmax,1,kmax)
         endif
      if (DEBUG.and.nodeid.eq.0) then
       print*, 'DWB10.6: AFTER call to btlu3j'
       print*, 'DWB10.6: jmax,kmax,dq(2,1,1) = ',jmax,kmax,dq(2,1,1)
       print*, 'DWB10.6: bjm(2,2,2,2), btc(2,2,2,2), bjp(2,2,2,2) = ',
     &  bjm(2,2,2,2),btc(2,2,2,2),bjp(2,2,2,2)
      endif 

c-----
c  Top of j-sweep loop
c-----
      if (DEBUG.and.nodeid.eq.0) then
        print*, ' DWB11.5: nodeid,njsmax,njsp = ',
     &    nodeid,njsmax,njsp
        print*, ' DWB11.5: nodeid,jmax,kmax,dq(2,1,1) = ',
     &    jmax,kmax,dq(2,1,1)
      endif

c default value of njsmax = 2 sweeps; njsp=njswp=2

      do 30 nj=1,njsmax
           if (DEBUG.and.nodeid.eq.0) then
            print*, 'DWB11: N5 - before ljsweep, nj,dq(2,1,1) = ',
     &            nj,dq(2,1,1)
            print*, 'DWB11: N5 - before ljsweep, nj,dq(10,1,2) = ',
     &            nj,dq(10,1,2)
           endif
            if(njsp .ge. nj) then
             call ljsweep(jmax,kmax,1,jmax,1,kmax,nodeid+1,nj,
     &        dq,s,btc,bjm,bjp,bkm,bkp
     &        )
              if (DEBUG.and.nodeid.eq.0) then
               print*, 'DWB9: after ljsweep, njsp,nj,dq(2,1,1) = ',
     &          njsp,nj,dq(2,1,1)
                print*, 'DWB9: after ljsweep, njsp,nj,dq(10,1,2) = ',
     &          njsp,nj,dq(10,1,2)
c               stop 'stop: after ljsweep in lrelax.f'
              endif

c-----
c  Pass boundary info between grids after each sweep
c    Note: default nbcimp = 1 in initia.f
c-----
               neqs = 3
               if( (nj/nbcimp)*nbcimp .eq. nj) then
               call bcimpds(jmax,kmax,neqs,q,dq,s,nodeid+1)
               if (DEBUG.and.nodeid.eq.0) then
                print*, 'DWB10: after bcimpds, nodeid,njsp,nj,dq(2,1,1), 
     &nt = ',
     &              nodeid,njsp,nj,dq(2,1,1),nt
                print*, 'DWB10: after bcimpds,nodeid,njsp,nj,dq(10,1,2),
     &nt = ',
     &              nodeid,njsp,nj,dq(10,1,2),nt
               endif
               endif
            endif
30    continue

      if (DEBUG.and.nt.ge.1) then 
c       stop 'stop: after bcimpds in lrelax.f'
      endif

c-----
c  Undo decomposition for j-sweeps
c-----
c njsp = 2 --> default defined in initia.f; defines number of 
c              sweeps in j or k direction
         if(njsp .gt. 0) then
          call btund3j(jmax,kmax,bjm,btc,bjp,
     &     1,jmax,1,kmax)
         endif

c-----
c  K-sweep: solve for lines of varying k and constant j.
c  First perform lu decomposition of each k-line for all grid
c  points in zones where k-sweeps will be done.
c-----
      nksmax = 0
c NOTE: nksp = 2, set in initia.f, identical to nkswp
         nksmax = max( nksmax, nksp )
         if(nksp .gt. 0) then
               call btlu3k(jmax,kmax,bkm,
     &            btc,bkp,1,jmax,1,kmax)
         endif

c-----
c  Top of k-sweep loop
c-----
      do 70 nk=1,nksmax
            if(nksp .ge. nk) then
                  call lksweep(jmax,kmax,1,jmax,1,kmax,
     &            nodeid+1,nk,dq,s,btc,bjm,bjp,bkm,bkp
     &            )
c-----
c  Pass boundary info between grids
c    Note: default nbcimp = 1 in initia.f
c-----
               neqs = 3
               if( (nk/nbcimp)*nbcimp .eq. nk)
     &         call bcimpds(jmax,kmax,neqs,q,dq,s,nodeid+1)

            endif
70    continue

c-----
c  Undo decomposition for k-sweeps
c-----
         if(nksp .gt. 0)
     &      call btund3k(jmax,kmax,bkm,btc,bkp,1,jmax,1,kmax)

c-----
c  End of lrelax
c-----
      return
      end
c
c
c****************************************************************
      subroutine ljsweep(jmax,kmax,jbeg,jend,kbeg,kend,nz,nj,dq,s,
     &                   btc,bjm,bjp,bkm,bkp)
c****************************************************************
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      __REAL btc(jmax,kmax,3,3), bjm(jmax,kmax,3,3),
     &     bjp(jmax,kmax,3,3), bkm(jmax,kmax,3,3),
     &     bkp(jmax,kmax,3,3), dq(jmax,kmax,3),
     &       s(jmax,kmax,3)
      dimension a(jkmax,3,3),b(jkmax,3,3),c(jkmax,3,3),rhs(jkmax,3),
     &            d(jkmax,3,3),e(jkmax,3,3)
      common/btri/a,b,c,rhs,
     &            d,e
      double precision, parameter :: half=0.5
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.

c-----
c  Form forcing function by moving off-line implicit terms*delta q
c  to rhs and then solve system one line at a time.
c-----
      if( (nj/2)*2 .eq. nj) then
c sweep in -k direction
         kb = kend
         ke = kbeg
         ki = -1
         if (DEBUG) then
          if (nz-1.eq.8) then
            print*,''
            print*, ' DWB99: In lrelax.f/ljsweep: node,nj,kb,ke,ki = ',
     &        nz-1,nj,kb,ke,ki
            print*,''
          endif
         endif
      else
c sweep in +k direction
         kb = kbeg
         ke = kend
         ki = 1
         if (DEBUG) then
          if (nz-1.eq.8) then
            print*,''
            print*, ' DWB98: In lrelax.f/ljsweep: node,nj,kb,ke,ki = ',
     &        nz-1,nj,kb,ke,ki
            print*,''
          endif
         endif
      endif

c sweep over k
      do 20 kk=kb,ke,ki
         k = nint( mod( float(kk+kmax), float(kmax)+.1 ) )
         km = max(k-1, 1)
         kp = min(k+1, kmax)
         rkm = float(k-1-km)
         rkp = float(kp-k-1)
c NOTE: 'sign(a,b)' returns the value of 'a' with the sign of 'b'
         swkm = half + sign(half, rkm)
         swkp = half + sign(half, rkp)

      if (DEBUG.and.nodeid.eq.0.and.k.le.10) then
       print*,'' 
       print*, 'nodeid  j   k     rhs_j1      rhs_j2      rhs_j3'
      endif
c sweep over j
         do 12 j=jbeg,jend
            rhs(j,1) = s(j,k,1)
     &      - ( bkm(j,k,1,1)*dq(j,km,1) + bkm(j,k,1,2)*dq(j,km,2)
     &        + bkm(j,k,1,3)*dq(j,km,3) )*swkm
     &      - ( bkp(j,k,1,1)*dq(j,kp,1) + bkp(j,k,1,2)*dq(j,kp,2)
     &        + bkp(j,k,1,3)*dq(j,kp,3) )*swkp
            rhs(j,2) = s(j,k,2)
     &      - ( bkm(j,k,2,1)*dq(j,km,1) + bkm(j,k,2,2)*dq(j,km,2)
     &        + bkm(j,k,2,3)*dq(j,km,3) )*swkm
     &      - ( bkp(j,k,2,1)*dq(j,kp,1) + bkp(j,k,2,2)*dq(j,kp,2)
     &        + bkp(j,k,2,3)*dq(j,kp,3) )*swkp
c            if (DEBUG.and.nodeid.eq.0) then
c             print*, 'DWB20: ------------------------'
c             print*, ' Iteration #: ',nj
c             print*, ' Node 8: j,k = ',j,k
c             print*, ' rhs(j,1),rhs(j,2) = ',rhs(j,1),rhs(j,2)
c             print*, ' s(j,k,1) = ',s(j,k,1)
c             print*, ' s(j,k,2) = ',s(j,k,2)
c             print*, ' bkm(j,k,2,1) = ',bkm(j,k,2,1)
c             print*, ' dq(j,km,1) = ',dq(j,km,1)
c             print*, ' bkm(j,k,2,2) = ',bkm(j,k,2,2)
c             print*, ' dq(j,km,2) = ',dq(j,km,2)
c             print*, ' bkm(j,k,2,3) = ',bkm(j,k,2,3)
c             print*, ' dq(j,km,3) = ',dq(j,km,3)
c             print*, ' swkm = ',swkm
c             print*, ' bkp(j,k,2,1) = ',bkp(j,k,2,1)
c             print*, ' dq(j,kp,1) = ',dq(j,kp,1)
c             print*, ' bkp(j,k,2,2) = ',bkp(j,k,2,2)
c             print*, ' dq(j,kp,2) = ',dq(j,kp,2)
c             print*, ' bkp(j,k,2,3) = ',bkp(j,k,2,3)
c             print*, ' dq(j,kp,3) = ',dq(j,kp,3)
c             print*, ' swkp = ',swkp
c           endif
            rhs(j,3) = s(j,k,3)
     &      - ( bkm(j,k,3,1)*dq(j,km,1) + bkm(j,k,3,2)*dq(j,km,2)
     &        + bkm(j,k,3,3)*dq(j,km,3) )*swkm
     &      - ( bkp(j,k,3,1)*dq(j,kp,1) + bkp(j,k,3,2)*dq(j,kp,2)
     &        + bkp(j,k,3,3)*dq(j,kp,3) )*swkp
      
      if (DEBUG.and.nodeid.eq.0.and.j.le.17.and.k.le.10) then
       write(*,999) nodeid,j,k,rhs(j,1),rhs(j,2),rhs(j,3)
999    format(3i5,1p3e13.5)
      endif

12       continue

      if (DEBUG.and.nodeid.eq.0.and.k.le.2) then
       print*,''
       print*, ' After do 12 loop in lrelax.f/ljsweep'
       print*, ' k = ',k
       print*, ' j  k    s(j,k,1)        s(j,k,2)       s(j,k,3)'
       do 961 j=1,17
        write(*,962) j,k,s(j,k,1),s(j,k,2),s(j,k,3)
962     format(2i3,1p3e15.7)
961    continue
       if (k.eq.2) then
        print*,''
c        print*, ' stop: after do 12 loop in lrelax.f/ljsweep'
c        stop 'stop: after do 12 loop in lrelax.f/ljsweep'
       endif
      endif

      if (DEBUG.and.nodeid.eq.0) then
       call flush(6)
       call flush(istdout)
c       stop 'stop: before do 15 loop in lrelax.f'
      endif

      if (DEBUG.and.nodeid.eq.0.and.k.le.10) then
       print*,''
       print*, ' j   k   m     ajm3       bjm3      cjm3'
       print*, ' jbeg,jend = ',jbeg,jend
      endif

         do 15 m=1,3
         do 15 j=jbeg,jend
            a(j,m,1) = bjm(j,k,m,1)
            b(j,m,1) = btc(j,k,m,1)
            c(j,m,1) = bjp(j,k,m,1)
            a(j,m,2) = bjm(j,k,m,2)
            b(j,m,2) = btc(j,k,m,2)
            c(j,m,2) = bjp(j,k,m,2)
            a(j,m,3) = bjm(j,k,m,3)
            b(j,m,3) = btc(j,k,m,3)
            c(j,m,3) = bjp(j,k,m,3)
      if (DEBUG.and.nodeid.eq.0.and.j.le.17.and.k.le.10) then
       write(*,990) j,k,m,a(j,m,3),b(j,m,3),c(j,m,3)
990    format(3i5,1p3e13.5)
      endif

15       continue

      if (DEBUG.and.nodeid.eq.0) then
       call flush(6)
       call flush(istdout)
c       stop 'stop: after do 15 loop in lrelax.f/lrelax'
      endif
c
c
c  Inlining of btso3j
c
c-----
c  Forward sweep
c-----
         j = jbeg
         d1 = b(j,1,1)* rhs(j,1)
         d2 = b(j,2,2)*(rhs(j,2) - b(j,2,1)*d1)
         d3 = b(j,3,3)*(rhs(j,3) - b(j,3,1)*d1 - b(j,3,2)*d2)
         rhs(j,3) = d3
         rhs(j,2) = d2 - b(j,2,3)*rhs(j,3)
         rhs(j,1) = d1 - b(j,1,3)*rhs(j,3) - b(j,1,2)*rhs(j,2)
c
         do 120 j=jbeg+1,jend
            do 110 n=1,3
               rhs(j,n) = rhs(j,n) - a(j,n,1)*rhs(j-1,1)
     &                             - a(j,n,2)*rhs(j-1,2)
     &                             - a(j,n,3)*rhs(j-1,3)
110         continue
            d1 = b(j,1,1)* rhs(j,1)
            d2 = b(j,2,2)*(rhs(j,2) - b(j,2,1)*d1)
            d3 = b(j,3,3)*(rhs(j,3) - b(j,3,1)*d1 - b(j,3,2)*d2)
            rhs(j,3) = d3
            rhs(j,2) = d2 - b(j,2,3)*rhs(j,3)
            rhs(j,1) = d1 - b(j,1,3)*rhs(j,3) - b(j,1,2)*rhs(j,2)
120      continue

      if (DEBUG.and.nodeid.eq.0.and.k.le.2) then
       print*,''
       print*, 'Forward sweep in lrelax.f/ljsweep'
       print*, '  jbeg,jend = ',jbeg,jend
       print*, '  k = ',k
       print*, '  nodeid = ',nodeid
       print*, ' j   k      rhs1       rhs2        rhs3'
       do 121 j=1,17
        write(*,122) j,k,rhs(j,1),rhs(j,2),rhs(j,3)
122     format(2i3,1p3e15.7)
121    continue
c       print*, ' stop: after do 120 in lrelax.f/ljsweep'
c       stop 'stop: after do 120 in lrelax.f/ljsweep'
      endif

c-----
c  Backward substitution
c-----
      if (DEBUG.and.nodeid.eq.0.and.k.le.2) then
       print*,''
       print*, ' j  k  n  rhs_jn'
      endif

         do 140 j=jend-1,jbeg,-1
            do 130 n=1,3
               rhs(j,n) = rhs(j,n) - c(j,n,1)*rhs(j+1,1)
     &                             - c(j,n,2)*rhs(j+1,2)
     &                             - c(j,n,3)*rhs(j+1,3)
               if (DEBUG.and.nodeid.eq.0.and.j.le.17.and.k.le.2) then
                write(*,131) j,k,n,rhs(j,n)
131             format(3i3,1pe13.5)
               endif
130         continue
140      continue

      if (DEBUG.and.nodeid.eq.0.and.k.ge.2) then
c       stop 'stop: after do 140 in lrelax.f/ljsweep'
      endif

c
c      if (nodeid.eq.8) then
c       print*,''
c       print*, 'DWB12: dq(2,1,1) = ',dq(2,1,1)
c       print*, 'DWB12: dq(10,1,2) = ',dq(10,1,2)
c       print*, 'DWB12: rhs(2,1) = ',rhs(2,1)
c       print*, 'DWB12: rhs(10,2) = ',rhs(10,2)
c       stop 'stop: before do 10 loop in lrelax.f'
c      endif
c
         do 10 j=jbeg,jend
            dq(j,k,1) = dq(j,k,1) + underr*(rhs(j,1) - dq(j,k,1))
            dq(j,k,2) = dq(j,k,2) + underr*(rhs(j,2) - dq(j,k,2))
            dq(j,k,3) = dq(j,k,3) + underr*(rhs(j,3) - dq(j,k,3))
10       continue

      if (DEBUG.and.nodeid.eq.0.and.k.le.2) then
       print*,''
       print*, ' After do 10 loop in lrelax.f/ljsweep'
       print*, 'underr = ',underr
       print*, ' nodeid = ',nodeid
       print*, ' j    k     dq1          dq2          dq3'
       do 222 j=1,17
        write(*,230) j,k,dq(j,k,1),dq(j,k,2),dq(j,k,3)
230     format(2i3,1p3e15.7)
222    continue
c       print*, ' stop: after do 10 loop in lrelax.f/ljsweep'
c       stop 'stop: after do 10 loop in lrelax.f/ljsweep'
      endif

20    continue

      if (DEBUG.and.nodeid.eq.0) then
       print*,''
       print*, ' nodeid  j   k   dq(jk1)     dq(jk2)      dq(jk3)'
       if (kmax.gt.10) then
        kmaxx = 10
       else
        kmaxx = kmax
      endif
      if (jmax.gt.17) then
       jmaxx = 17
      else
       jmaxx = jmax
      endif
       do 1011 k=1,kmaxx
       do 1011 j=1,jmaxx
        write(*,1012) nodeid,j,k,dq(j,k,1),dq(j,k,2),dq(j,k,3)
1012    format(3i5,1p3e13.5)
1011   continue
       call flush(6)
       call flush(istdout)
c       stop 'stop: at end of lrelax.f/ljsweep'
      endif

c-----
c  End of ljsweep
c-----
      return
      end
c
c
c****************************************************************
      subroutine lksweep(jmax,kmax,jbeg,jend,kbeg,kend,nz,nk,dq,s,
     &                   btc,bjm,bjp,bkm,bkp)
c****************************************************************
#include "common.f"
#include "mpif.h"
#include "mpi_params.f"
      __REAL btc(jmax,kmax,3,3), bjm(jmax,kmax,3,3),
     &     bjp(jmax,kmax,3,3), bkm(jmax,kmax,3,3),
     &     bkp(jmax,kmax,3,3), dq(jmax,kmax,3),
     &       s(jmax,kmax,3)
      dimension a(jkmax,3,3),b(jkmax,3,3),c(jkmax,3,3),rhs(jkmax,3),
     &            d(jkmax,3,3),e(jkmax,3,3)
      common/btri/a,b,c,rhs,
     &            d,e
      double precision, parameter :: half=0.5
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.
c-----
c  Non-periodic
c-----
c-----
c  Form forcing function by moving off-line implicit terms*delta q
c  to rhs and then solve system one line at a time.
c-----
      if( (nk/2)*2 .eq. nk) then
         jb = jend
         je = jbeg
         ji = -1
       else
         jb = jbeg
         je = jend
         ji = 1
      endif
      do 40 j=jb,je,ji
         jm = max(j-1, 1)
         jp = min(j+1, jmax)
         rjm = float(j-1-jm)
         rjp = float(jp-j-1)
         swjm = half + sign(half, rjm)
         swjp = half + sign(half, rjp)
         do 26 k=kbeg,kend
            rhs(k,1) = s(j,k,1)
     &      - ( bjm(j,k,1,1)*dq(jm,k,1) + bjm(j,k,1,2)*dq(jm,k,2)
     &        + bjm(j,k,1,3)*dq(jm,k,3) )*swjm
     &      - ( bjp(j,k,1,1)*dq(jp,k,1) + bjp(j,k,1,2)*dq(jp,k,2)
     &        + bjp(j,k,1,3)*dq(jp,k,3) )*swjp
            rhs(k,2) = s(j,k,2)
     &      - ( bjm(j,k,2,1)*dq(jm,k,1) + bjm(j,k,2,2)*dq(jm,k,2)
     &        + bjm(j,k,2,3)*dq(jm,k,3) )*swjm
     &      - ( bjp(j,k,2,1)*dq(jp,k,1) + bjp(j,k,2,2)*dq(jp,k,2)
     &        + bjp(j,k,2,3)*dq(jp,k,3) )*swjp
            rhs(k,3) = s(j,k,3)
     &      - ( bjm(j,k,3,1)*dq(jm,k,1) + bjm(j,k,3,2)*dq(jm,k,2)
     &        + bjm(j,k,3,3)*dq(jm,k,3) )*swjm
     &      - ( bjp(j,k,3,1)*dq(jp,k,1) + bjp(j,k,3,2)*dq(jp,k,2)
     &        + bjp(j,k,3,3)*dq(jp,k,3) )*swjp
26       continue
c
      if (DEBUG.and.nodeid.eq.0) then
       print*,''
       print*, ' m  k     akm1       bkm1       ckm1'
      endif
         do 35 m=1,3
         do 35 k=kbeg,kend
            a(k,m,1) = bkm(j,k,m,1)
c DWB: ERROR
c b(k,3,1) = 1, but should be 0 (means btc(j,k,m,1) is incorrect
            b(k,m,1) = btc(j,k,m,1)
            c(k,m,1) = bkp(j,k,m,1)
            a(k,m,2) = bkm(j,k,m,2)
            b(k,m,2) = btc(j,k,m,2)
            c(k,m,2) = bkp(j,k,m,2)
            a(k,m,3) = bkm(j,k,m,3)
            b(k,m,3) = btc(j,k,m,3)
            c(k,m,3) = bkp(j,k,m,3)
      if (DEBUG.and.nodeid.eq.0.and.m.eq.3) then
       write(*,2090) m,k,a(k,m,1),b(k,m,1),c(k,m,1)
2090   format(i3,2x,i3,3e13.5)
      endif
35       continue
c
c
c  Inlining of btso3k
c

      if (DEBUG.and.nodeid.eq.0) then
       print*,''
c       print*,'nt node  k    b(k11)    b(k12)   b(k13)     b(k21)    b(k22)
c     &     b(k23)    b(k31)    b(k32)    b(k33)'   
c       print*,'nt  node  k    rhs(k1)    rhs(k2)   rhs(k3)'
c       print*, ' time step = ',nt
       call flush(istdout)
       call flush(6)
      endif

         k = kbeg

      if (DEBUG.and.nodeid.eq.0) then
c       write(*,999) nt,nodeid,k,b(k,1,1),b(k,1,2),b(k,1,3),b(k,2,1),
c     & b(k,2,2),b(k,2,3),b(k,3,1),b(k,3,2),b(k,3,3) 
999   format(i3,i4,i5,1p9e12.4)
       print*,'nt  node  k    rhs(k1)    rhs(k2)   rhs(k3)'
       print*, ' time step = ',nt
       write(*,998) 120,nt,nodeid,k,rhs(k,1),rhs(k,2),rhs(k,3)
998   format(i3,2x,i3,i4,i5,1p3e12.4)
      call flush(istdout)
      call flush(6)
c       stop 'stop: in lrelax.f/lksweep'
      endif

         d1 = b(k,1,1)* rhs(k,1)
         d2 = b(k,2,2)*(rhs(k,2) - b(k,2,1)*d1)
         d3 = b(k,3,3)*(rhs(k,3) - b(k,3,1)*d1 - b(k,3,2)*d2)
         rhs(k,3) = d3
         rhs(k,2) = d2 - b(k,2,3)*rhs(k,3)
         rhs(k,1) = d1 - b(k,1,3)*rhs(k,3) - b(k,1,2)*rhs(k,2)
c
         do 120 k=kbeg+1,kend
            do 110 n=1,3
               rhs(k,n) = rhs(k,n) - a(k,n,1)*rhs(k-1,1)
     &                             - a(k,n,2)*rhs(k-1,2)
     &                             - a(k,n,3)*rhs(k-1,3)
110          continue
            d1 = b(k,1,1)* rhs(k,1)
            d2 = b(k,2,2)*(rhs(k,2) - b(k,2,1)*d1)
            d3 = b(k,3,3)*(rhs(k,3) - b(k,3,1)*d1 - b(k,3,2)*d2)
            rhs(k,3) = d3
            rhs(k,2) = d2 - b(k,2,3)*rhs(k,3)
            rhs(k,1) = d1 - b(k,1,3)*rhs(k,3) - b(k,1,2)*rhs(k,2)
      if (DEBUG.and.nodeid.eq.0) then
       write(*,997) 120,nt,nodeid,k,rhs(k,1),rhs(k,2),rhs(k,3)
997   format(i3,2x,i3,i4,i5,1p3e12.4)
      endif
120       continue

      if (DEBUG) then
       stop 'stop: after do 120 in lrelax.f/lksweep'
      endif
c-----
c  Backward substitution
c-----
         do 140 k=kend-1,kbeg,-1
            do 130 n=1,3
               rhs(k,n) = rhs(k,n) - c(k,n,1)*rhs(k+1,1)
     &                             - c(k,n,2)*rhs(k+1,2)
     &                             - c(k,n,3)*rhs(k+1,3)
130          continue
140       continue
c
c End of btso3k inlining
c
         do 30 k=kbeg,kend
            dq(j,k,1) = dq(j,k,1) + underr*(rhs(k,1) - dq(j,k,1))
            dq(j,k,2) = dq(j,k,2) + underr*(rhs(k,2) - dq(j,k,2))
            dq(j,k,3) = dq(j,k,3) + underr*(rhs(k,3) - dq(j,k,3))
30       continue
40    continue
c-----
c  End of lksweep
c-----
      return
      end
