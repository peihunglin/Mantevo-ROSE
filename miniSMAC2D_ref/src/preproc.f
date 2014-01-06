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
c  1. readgo
c  2. switchzy
c  4. addk2d
c  5. cube
c  6. partishn
c  7. average
c  8. break
c  9. patch_2d
c 10. min_surf_vol_ratio
c 11. subgrid_dimensions
c 12. grid_point_comparison
c 13. global_indices
c 14. write_subgrids
c 15. get_global_index
c 16. out1planetoplt3d_break
c
c************************************************************
c
c
c************************************************************************
      subroutine readgo(iunit,x,y,z,jdim,kdim,ldim,lindex)
c************************************************************************
#include "forttype.h"
c
c Purpose: This subroutine reads the original x,y,z grid values 
c
c Author: D. W. Barnette, SNL
c
      __REAL x(jdim,kdim,ldim),y(jdim,kdim,ldim),z(jdim,kdim,ldim)
      logical DEBUG
c
c   ^ y,k
c   |
c   |
c   |
c   -----------> x,j
c   \
c    \
c     \
c      V  z,l
c where x-y is the geometry plane

      DEBUG = .false.
c      DEBUG = .true.
c
c read 2d array into x and y, set z to zero
       rewind(iunit)
       read(iunit,*) idum
       read(iunit,*) idum1,idum2
       read(iunit,*)
     1  (((x(j,k,l),j=1,jdim),k=1,kdim),l=1,lindex),
     2  (((y(j,k,l),j=1,jdim),k=1,kdim),l=1,lindex)

c set z equal to 0
       do 5 j=1,jdim
       do 5 k=1,kdim
       do 5 l=1,ldim
         z(j,k,l)=0.0
5      continue
       z = 0.0

c assign values to x and y for l=2,3
       do 10 j=1,jdim
       do 10 k=1,kdim
       do 10 l=2,ldim
       x(j,k,l) = x(j,k,l-1)
       y(j,k,l) = y(j,k,l-1)
10     continue
c
      if (DEBUG) then
        print*,''
        print*, '=================================================='
        print*,''
        print*, ' In sub readgo, grid corner values are:'
        print*, ' jdim,kdim,ldim,lindex = ',jdim,kdim,ldim,lindex
        print*, ' 1. xyz(1,1,1) = ',
     1    x(1,1,1), y(1,1,1),z(1,1,1)

        print*, ' 2. xyz(jmx,1,1), y(jmx,1,1) = ',
     1    x(jdim,1,1),y(jdim,1,1),z(jdim,1,1)

        print*, ' 3. xyz(1,kmx,1), y(1,kmax,1) = ',
     1    x(1,kdim,1),y(1,kdim,1),z(1,kdim,1)

        print*, ' 4. xyz(1,1,ldim) = ',
     1    x(1,1,ldim),y(1,1,ldim),z(1,1,ldim)

        print*, ' 5. xyz(jmx,kmx,1) = ',
     1    x(jdim,kdim,1),y(jdim,kdim,1),z(jdim,kdim,1)

        print*, ' 6. xyz(jmx,kmx,lmx) = ',
     1    x(jdim,kdim,ldim),y(jdim,kdim,ldim),z(jdim,kdim,ldim)

        print*,''
        print*, '=================================================='
      endif
c
      return
      end  
c
c
c======================================================================|
      subroutine switchzy(x,y,z,jkl,jdim,kdim,ldim)
c
c  Purpose: This subroutine switches y and z as follows. 
c               z(new) =  y(old)
c               y(new) = -z(old)
c           This is equivalent to a CLOCKWISE rotation
c           about the +x axis
c   Author: D. W. Barnette, Org. 1556, SNL
c  Written: 1-28-91
c Comments:
c
#include "forttype.h"

      __REAL x(jkl),y(jkl),z(jkl)
c switch k and l dimensions
      __REAL xtemp(jdim,ldim,kdim),ytemp(jdim,ldim,kdim),
     1 ztemp(jdim,ldim,kdim)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.
      
      if (DEBUG) then
        print*,''
        print*, ' >> In sub. switchzy in preproc.f'
        print*,''
      endif
c
c
c      print *,'  Initial block dimensions are: ',jdim,kdim,ldim
c
c Switch y and z
      do 10 j=1,jdim
      do 10 k=1,kdim
      do 10 l=1,ldim
c goes from 1 to lmax
        index1 = (l-1)*kdim*jdim + (k-1)*jdim + j 
        xtemp(j,l,k)=x(index1)
        ytemp(j,l,k)=-z(index1)
        ztemp(j,l,k)=y(index1)
10    continue

      ksave = kdim
      kdim = ldim
      ldim = ksave
c
      do 12 j=1,jdim
      do 12 k=1,kdim
      do 12 l=1,ldim
        index1 = (l-1)*kdim*jdim + (k-1)*jdim + j 
        x(index1) = xtemp(j,k,l)
        y(index1) = ytemp(j,k,l)
        z(index1) = ztemp(j,k,l)
12    continue
c
      print*,''
      print*,'  Final block dimensions after call switchzy are: ',
     1 jdim,kdim,ldim
      print*,''
c      
      return
      end
c
c
c-----------------------------------------------------------------
      subroutine addk2d(x,y,z,jkl,jdim,kdim,ldim)
c
c  Purpose: This subroutine adds two additional k planes, one at
c           k=1 and one at k=kmax. It is assumed that j,k,l and
c           x,y,z have been defined upon entry into this subroutine
c           for use by the Navier-Stokes code F3D. The subroutine
c           is used to add additional planes of data for two-dimensional
c           flows. This subroutine should not be called for three-
c           dimensional flows.
c   Author: D. W. Barnette, Org. 1556, SNL
c  Written: 1-28-91
c Comments:
c
c
#include "forttype.h"

c      dimension x(jkl),y(jkl),z(jkl)
      __REAL x(jdim,kdim,ldim),y(jdim,kdim,ldim),
     1 z(jdim,kdim,ldim)

c
c      print *,'  Initial block dimensions are: ',jdim,kdim,ldim
c
c Add extra planes
c Following is for 2d type flows.
c
c Scoot everything up one indici on k
      do 90 j=1,jdim
      do 90 k=2,kdim
      do 90 l=1,ldim
c      index1 = (l-1)*kdim*jdim + (k-1)*jdim + j
c      index2 = (l-1)*kdim*jdim + (k)*jdim + j
c      x(index2)=x(index1)
c      y(index2)=y(index1) + 1.0
c      z(index2)=z(index1)
      x(j,k,l) = x(j,k-1,l)
      y(j,k,l) = y(j,k-1,l) + 1.
      z(j,k,l) = z(j,k-1,l)
90    continue

c Bump kdim
c      kdim=kdim+2
c
      print *,''
      print *,'  Additional k = 1 & k =',kdim,
     1 ' planes defined for original grid'
c
      print *,' sub addk2d: Final block dimensions are: ',
     1 jdim,kdim,ldim
c      
      return
      end
c
c
c-----------------------------------------------------------------
      subroutine cube(x,y,z,jdim,kdim,ldim)
c
c  Purpose: This subroutine writes out select grid points
c           for verification
c   Author: D. W. Barnette, Org. 1556, SNL
c  Written: 1-28-91
c Comments:
c
#include "forttype.h"
c
      __REAL x(jdim,kdim,ldim),
     1 y(jdim,kdim,ldim),z(jdim,kdim,ldim)
c
c
      print *,''
      print *,'Block corner locations with/without surrounding points:'
      print *,''
      print *,'     Jmax, Kmax, Lmax = ',jdim, kdim, ldim
      print *,''
      print *,'               6____________________7'
      print *,'               /|                 /|'
      print *,'          (K) / |                / |'
      print *,'             /  |               /  |'
      print *,'          2 /__________________/3  |'
      print *,'           |    |             |    |'
      print *,'           |   5|-------------|----|8'
      print *,'           |   /              |   /'
      print *,'       (L) ^  /               |  /'
      print *,'           | /                | /'
      print *,'           |/                 |/'
      print *,'         1 |------>-----------/ 4'
      print *,'                   (J)'
      print *,''
      print *,'  point  j   k   l      x          y          z'
      print *,' ------ --- --- --- ---------- ---------- ----------'
c
c Outer block indices
      j0=1
      k0=1
      l0=1
      jdimm=jdim
      kdimm=kdim
      ldimm=ldim
c
      print *,'Corner grid points:'
      write(6,500) 1,j0,k0,l0,
     1 x(j0,k0,l0),y(j0,k0,l0),z(j0,k0,l0)
      write(6,500) 2,j0,k0,ldimm,
     1 x(j0,k0,ldimm),y(j0,k0,ldimm),
     2 z(j0,k0,ldimm)
      write(6,500) 3,jdimm,k0,ldimm,
     1 x(jdimm,k0,ldimm),
     2 y(jdimm,k0,ldimm),
     3 z(jdimm,k0,ldimm)
      write(6,500) 4,jdimm,k0,l0,
     1 x(jdimm,k0,l0),y(jdimm,k0,l0),
     2 z(jdimm,k0,l0)
      write(6,500) 5,j0,kdimm,l0,
     1 x(j0,kdimm,l0),y(j0,kdimm,l0),
     2 z(j0,kdimm,l0)
      write(6,500) 6,j0,kdimm,ldimm,
     1 x(j0,kdimm,ldimm),
     2 y(j0,kdimm,ldimm),
     2 z(j0,kdimm,ldimm)
      write(6,500) 7,jdimm,kdimm,ldimm,
     1 x(jdimm,kdimm,ldimm),
     2 y(jdimm,kdimm,ldimm),
     3 z(jdimm,kdimm,ldimm)
      write(6,500) 8,jdimm,kdimm,l0,
     1 x(jdimm,kdimm,l0),
     2 y(jdimm,kdimm,l0),
     3 z(jdimm,kdimm,l0)
c
      print *,' MID points on each surface:'
      jmid=(jdim+1)/2
      kmid=(kdim+1)/2
      lmid=(ldim+1)/2
      print *,' Jmid, Kmid, Lmid = ',jmid,kmid,lmid
c 1,kmid,lmid
      write(6,500) 1,1,kmid,lmid,
     1 x(1,kmid,lmid),
     2 y(1,kmid,lmid),
     3 z(1,kmid,lmid)
c jmax,kmid,lmid
      write(6,500) 2,jdim,kmid,lmid,
     1 x(jdim,kmid,lmid),
     2 y(jdim,kmid,lmid),
     3 z(jdim,kmid,lmid)
c jmid,1,lmid
      write(6,500) 3,jmid,1,lmid,
     1 x(jmid,1,lmid),
     2 y(jmid,1,lmid),
     3 z(jmid,1,lmid)
c jmid,kmax,lmid
      write(6,500) 4,jmid,kdim,lmid,
     1 x(jmid,kdim,lmid),
     2 y(jmid,kdim,lmid),
     3 z(jmid,kdim,lmid)
c jmid,kmid,1
      write(6,500) 5,jmid,kmid,1,
     1 x(jmid,kmid,1),
     2 y(jmid,kmid,1),
     3 z(jmid,kmid,1)
c jmid,kmid,lmax
      write(6,500) 6,jmid,kmid,ldim,
     1 x(jmid,kmid,ldim),
     2 y(jmid,kmid,ldim),
     3 z(jmid,kmid,ldim)
c
      print *,''
      print *,' POINT VALUES ALONG BLOCK EDGES'
      print *,'         j   k   l      x          y          z'
      print *,' ------ --- --- --- ---------- ---------- ----------'
      print *,'All points from 1 to 2:'
      do 10 j=1,1
      do 10 k=1,1
      do 10 l=1,ldim
      write(6,502) j,k,l,
     1 x(j,k,l),
     2 y(j,k,l),
     3 z(j,k,l)
10    continue
c
      print *,'All points from 2 to 3:'
      do 11 j=1,jdim
      do 11 k=1,1
      do 11 l=ldim,ldim
      write(6,502) j,k,l,
     1 x(j,k,l),
     2 y(j,k,l),
     3 z(j,k,l)
11    continue
c
      print *,''
      print *,'All points from 4 to 3:'
      do 12 j=jdim,jdim
      do 12 k=1,1
      do 12 l=1,ldim
      write(6,502) j,k,l,
     1 x(j,k,l),
     2 y(j,k,l),
     3 z(j,k,l)
12    continue
c
      print *,''
      print *,'All points from 1 to 4:'
      do 13 j=1,jdim
      do 13 k=1,1
      do 13 l=1,1
      write(6,502) j,k,l,
     1 x(j,k,l),
     2 y(j,k,l),
     3 z(j,k,l)
13    continue
c
      print *,''
      print *,'All points from 1 to 5:'
      do 14 j=1,1
      do 14 k=1,kdim
      do 14 l=1,1
      write(6,502) j,k,l,
     1 x(j,k,l),
     2 y(j,k,l),
     3 z(j,k,l)
14    continue
c
      print *,''
      print *,'All points from 2 to 6:'
      do 15 j=1,1
      do 15 k=1,kdim
      do 15 l=ldim,ldim
      write(6,502) j,k,l,
     1 x(j,k,l),
     2 y(j,k,l),
     3 z(j,k,l)
15    continue
c
      print *,''
      print *,'All points from 3 to 7:'
      do 16 j=jdim,jdim
      do 16 k=1,kdim
      do 16 l=ldim,ldim
      write(6,502) j,k,l,
     1 x(j,k,l),
     2 y(j,k,l),
     3 z(j,k,l)
16    continue
      print *,''
      print *,'All points from 4 to 8:'
      do 17 j=jdim,jdim
      do 17 k=1,kdim
      do 17 l=1,1
      write(6,502) j,k,l,
     1 x(j,k,l),
     2 y(j,k,l),
     3 z(j,k,l)
17    continue
c
      print *,''
      print *,'All points from 6 to 7:'
      do 18 j=1,jdim
      do 18 k=kdim,kdim
      do 18 l=ldim,ldim
      write(6,502) j,k,l,
     1 x(j,k,l),
     2 y(j,k,l),
     3 z(j,k,l)
18    continue
c
      print *,''
      print *,'All points from 5 to 8:'
      do 19 j=1,jdim
      do 19 k=kdim,kdim
      do 19 l=1,1
      write(6,502) j,k,l,
     1 x(j,k,l),
     2 y(j,k,l),
     3 z(j,k,l)
19    continue
c
      print *,''
      print *,'All points from 5 to 6:'
      do 20 j=1,1
      do 20 k=kdim,kdim
      do 20 l=1,ldim
      write(6,502) j,k,l,
     1 x(j,k,l),
     2 y(j,k,l),
     3 z(j,k,l)
20    continue
      print *,''
      print *,'All points from 8 to 7:'
      do 21 j=jdim,jdim
      do 21 k=kdim,kdim
      do 21 l=1,ldim
      write(6,502) j,k,l,
     1 x(j,k,l),
     2 y(j,k,l),
     3 z(j,k,l)
21    continue
c
c
500   format(i3,5x,i3,1x,i3,1x,i3,1x,g10.3,1x,g10.3,1x,g10.3)
502   format(8x,i3,1x,i3,1x,i3,1x,g10.3,1x,g10.3,1x,g10.3)
c
      return
      end
c
c
c-----------------------------------------------------------------
      subroutine partishn(iunit,form,mzone,xorig,yorig,zorig,
     & jdim,kdim,ldim,idimj,idimk,idiml,jkl,nodeid,numprocs,j_subgrid,
     & k_subgrid,l_subgrid,jindex_global,kindex_global,lindex_global,
     & ioverlap)
c
c Purpose:
c     partitions original grid into user-specified (numprocs)
c     number of subgrids which are load balanced and optimized for maximum
c     volume and minimum surface area; subrids are overlapped for
c     subgrid-to-subgrid communication
c
c Output:
c     file bcpzn.dat which tells which subgrid communicates with which
c     other subgrid
c
#include "forttype.h"
c     
      __REAL xorig(jkl),yorig(jkl),zorig(jkl)
      dimension j_subgrid(numprocs),k_subgrid(numprocs),
     1 l_subgrid(numprocs)
      dimension jindex_global(numprocs),kindex_global(numprocs),
     &  lindex_global(numprocs)
      dimension ipnum(numprocs)
      character*10 meshname
c
c
c calculate total grid points, average number of grid points per 
c  grid, etc.
      nzonetot = 1
      nzne = 1
      meshname='2D_mesh'
      call average(nzonetot,nzne,meshname,numprocs,jdim,kdim,ldim)
c
c breakup grids
      call break(nzonetot,meshname,nodeid,numprocs,iunit,form,
     1 ioverlap,ipnum(1),jdim,kdim,ldim,idimj,idimk,idiml,
     2 j_subgrid(1),k_subgrid(1),l_subgrid(1),jkl,
     3 xorig(1),yorig(1),zorig(1),
     4 jindex_global(1),kindex_global(1),lindex_global(1)
     5 )
c   
c link patched grids to subgrids
       call link_overlap(ioverlap,ipnum(1),nodeid,numprocs,
     1 jdim,kdim,ldim,idimj,idimk,idiml,
     2 j_subgrid(1),k_subgrid(1),l_subgrid(1))
c

      write(6,*) '                                                     '
      write(6,*) '------------------------------------------------     '
      write(6,*) '  Subgrids output to file: GRIDS2D (PLOT3D format)   '
      write(6,*) '     and boundary condition link file: 2D_PATCH_TABLE'
      write(6,*) '  in PLOT3D format.                                  '
      write(6,*) '------------------------------------------------     '
      write(6,*) '                                                     '
      write(6,*) '   >> Grid partitioning finished <<                  '
      write(6,*) '                                                     '
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c NONE                                                                  |
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end

c======================================================================|
      subroutine average(nzone,ndim,meshname,numprocs,jdim,kdim,ldim)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c       This subroutine takes in all grids and computes average number
c       of grid points per subgrid.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c       D. W. Barnette
c       Sandia National Laboratories
c       Org. 9221, MS 1111
c       Albuquerque, NM 87112
c       email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c       breakup_for_fun
c       breakup_grids
c       chooser
c
c Called routines, in order of appearance:
c       NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c       02-13-96        Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
#include "forttype.h"
c
      common/averages/isum,gpavg,igptotal
c
      character*10 meshname
c
c
      igptotal=jdim*kdim*ldim
      isum = igptotal
c
      gpavg=float(igptotal)/float(numprocs)
      write(6,100)
c
      write(6,110) meshname,jdim,kdim,ldim,igptotal
c
      write(6,115) isum,numprocs,gpavg
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
100   format( //,3x,'Grid',5x,'Jmax',3x,'Kmax',3x,'Lmax',6x,'Total')    
110   format(a10,2x,i4,3x,i4,3x,i4,3x,i8)                               
115   format(/,'Total number of grid points = ',i10,//,                 
     1 'Total number of final zones (or subgrids) = ',i3,//,            
     2 'Average number of grid points/grid = ',f10.2,//)                
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c
c=======================================================================
      subroutine break(nzone,meshname,nodeid,numprocs,iunit,form,
     1 ioverlap,ipnum,jdim,kdim,ldim,idimj,idimk,idiml,
     2 j_subgrid,k_subgrid,l_subgrid,jklmax,xorig,yorig,zorig,
     3 jindex_global,kindex_global,lindex_global
     4 )
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c       This subroutine breaks up grids into a user-specified number
c       of zones.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c       D. W. Barnette
c       Sandia National Laboratories
c       Org. 9221, MS 1111
c       Albuquerque, NM 87112
c       email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c       partishn 
c
c Called routines, in order of appearance:
c       min_surf_vol_ratio
c       subgrid_dimensions
c       grid_point_comparison
c       global_indices
c       get_global_index
c       write_subgrids
c       get_global_index
c       write_subgrids
c       out1planetoplt3d_break
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c       02-13-96        Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   Guess at these; error catching will inform user if not correct
c      parameter(iicntr_max=maxdim*500,
c     1 ibcntr_max=iicntr_max)
c
#include "forttype.h"

      dimension j_subgrid(numprocs),k_subgrid(numprocs),
     1 l_subgrid(numprocs)
      dimension ipnum(numprocs)
      dimension jindex_global(numprocs),kindex_global(numprocs),
     1  lindex_global(numprocs)
      __REAL xorig(jklmax),yorig(jklmax),zorig(jklmax)

      common/averages/isum,gpavg,igptotal
      common/parallel/ippzone
c
      character*10 meshname
c number of subgrids (or cubes) formed in original grid 
      integer idimj,idimk,idiml
c
      logical form
c
c original grid dimensions
      jmax_compout = jdim
      kmax_compout = kdim
      lmax_compout = ldim

      write(6,*) numprocs
      write(6,*) '                                                     '
c Calculate grid points per processor node
      gppproc=isum/numprocs
c Print parameters
      do 41 i=1,nzone
      write(6,42) i,jmax_compout*kmax_compout*lmax_compout
41    continue
      write(6,43) isum,gpavg,gppproc

c
c Calculate processors per zone
      isump=0
      do 15 i=1,nzone
       ippzone=nint(float(igptotal)/real(gppproc))
       if(ippzone.eq.0) ippzone=1
       isump=isump+ippzone
       write(6,25) i,ippzone
15    continue
      write(6,80) numprocs,isump

c
c Determine exact number jxkxl per processor for each zone.
c Optimize surface to volume ratio for efficient 'communications
c  to compute' ratio, given the number of processors for each zone.
c
      write(6,*) '                                                     '
      write(6,*) ' Calculate j*k*l values for subgrids for minimum     '
      write(6,*) '   area-to-volume ratio.                             '
c
c
c The following logic is applicable to dividing the grid into cubes.
c  This is optimal for load-balancing 3D grids.
c
c
       write(6,301) meshname
c
c Calculate minimum surface/volume ratio, given the number of desired 
c  processors
c   (determines idimj, idimk, etc.)
      call min_surf_vol_ratio(1,ippzone,jmax_compout,
     1 kmax_compout,lmax_compout,
     2 jpts_per_subgrid,kpts_per_subgrid,lpts_per_subgrid,nzone,
     3 nodeid,numprocs,idimj,idimk,idiml,jdim,kdim,ldim)
c
c Divide the grid and determine dimensions
c  (determines j_subgrid, k_subgrid, etc.
       call subgrid_dimensions(1,jmax_compout,kmax_compout,
     1 lmax_compout,jpts_per_subgrid,kpts_per_subgrid,
     2 lpts_per_subgrid,nodeid,numprocs,idimj,idimk,idiml,
     3 j_subgrid(1),k_subgrid(1),l_subgrid(1))
c
c Sum of grid points in the subdivided grids must add up to number in 
c  original grid before subdividing
       call grid_point_comparison(1,jmax_compout,kmax_compout,
     1 lmax_compout,nodeid,numprocs,idimj,idimk,idiml,
     2 j_subgrid(1),k_subgrid(1),l_subgrid(1))
c
300   continue
c
c assign processor numbers for all zones in the array ipnum()
      itemp=0
      do 302 j=1,idimj
      do 302 k=1,idimk
      do 302 l=1,idiml
       index=(j-1)*idimk*idiml+(k-1)*idiml+l
       ipnum(index)=itemp+1
       itemp=ipnum(index)
302   continue
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Set global indices
      call global_indices(nzone,meshname,nodeid,numprocs,
     & j_subgrid(1),k_subgrid(1),l_subgrid(1),idimj,idimk,idiml,
     2 jindex_global(1),kindex_global(1),lindex_global(1))
c
c establish new grids for each zone, based on new dimensions
c
      write(6,*) '                                                     '
      write(6,*) '                                                     '
      write(6,*) ' Breaking up grids into subgrids WITH NO overlap:    '
      write(6,295)
c
      iproc=0
c
       ii=0
       do 325 j=1,idimj
       do 325 k=1,idimk
       do 325 l=1,idiml
        index=(j-1)*idimk*idiml+(k-1)*idiml+l
c
c set ioverlap
      ioverlap=0
c
c get global indices
        call get_global_index(1,j,k,l,1,1,1,
     1  jp3d_min,jp3d_max,kp3d_min,kp3d_max,lp3d_min,lp3d_max,
     2  ioverlap,nodeid,numprocs,
     3  j_subgrid(1),k_subgrid(1),l_subgrid(1),
     4  idimj,idimk,idiml,jdim,kdim,ldim,
     5  jindex_global(1),kindex_global(1),lindex_global(1))
c
        itotal=j_subgrid(index)*k_subgrid(index)*
     1  l_subgrid(index)
c
c
c number the subgrids
        iproc=iproc+1
        ii=ii+1
c
        write(6,333) 1,ii,iproc,j_subgrid(index),
     1  k_subgrid(index),l_subgrid(index),itotal,jp3d_min,
     2  jp3d_max,kp3d_min,kp3d_max,lp3d_min,lp3d_max
325    continue
c
       iproc_total=iproc+1
c
       if(1.lt.nzone) then
        write(6,*) '----------------------------------------------------
     1----------------------------'
       endif
c
320   continue
c
      iproc_total=iproc_total-1
c
c
c Write out Zone, Subzone, J-range, K-range, L-range (relative values)
c   for writing grids to file GRIDS
c
      write(6,*) '                                                     '
      write(6,*) '                                                     '
      write(6,*) ' Breaking up grids into subgrids WITH overlap:       '
      write(6,295)
c
      iproc=0
c
      ioverlap=1
c
      ii=0
      do 425 j=1,idimj
      do 425 k=1,idimk
      do 425 l=1,idiml
       index=(j-1)*idimk*idiml+(k-1)*idiml+l
c
c get global indices
       call get_global_index(1,j,k,l,1,1,1,jp3d_min,jp3d_max,
     1 kp3d_min,kp3d_max,lp3d_min,lp3d_max,ioverlap,nodeid,numprocs,
     2  j_subgrid(1),k_subgrid(1),l_subgrid(1),idimj,idimk,idiml,
     3  jdim,kdim,ldim,
     4  jindex_global(1),kindex_global(1),lindex_global(1))


c j values
       if (j.lt.idimj) then
        j_subgrid(index)=j_subgrid(index)+ioverlap
       endif
c k values
       if (k.lt.idimk) then
        k_subgrid(index)=k_subgrid(index)+ioverlap
       endif
c l values
       if (l.lt.idiml) then
        l_subgrid(index)=l_subgrid(index)+ioverlap
       endif
c
       itotal=(jp3d_max - jp3d_min + 1)*
     1        (kp3d_max - kp3d_min + 1)*
     2        (lp3d_max - lp3d_min + 1)
c
c
c number the subgrids
       iproc=iproc+1
       ii=ii+1
c
       write(6,433) 1,ii,iproc,j_subgrid(index),
     1 k_subgrid(index),l_subgrid(index),itotal,jp3d_min,
     2 jp3d_max,kp3d_min,kp3d_max,lp3d_min,lp3d_max
425     continue
c
       iproc_total=iproc+1
c
       if(1.lt.nzone) then
        write(6,*) '----------------------------------------------------
     1---------------------------'
       endif
c
420   continue

c
      iproc_total=iproc_total-1
c 
c+++++++++++++++++++++
c
c Write new subgrid gridpoint values x,y,z to one file in PLOT3D format
       write(6,304)
c
       call write_subgrids(nzone,iunit,form,nodeid,numprocs,
     1  ioverlap,idimj,idimk,idiml,jdim,kdim,ldim,
     2  j_subgrid(1),k_subgrid(1),l_subgrid(1),
     3  ippzone,jklmax,
     4  xorig(1),yorig(1),zorig(1),
     5  jindex_global(1),kindex_global(1),lindex_global(1))
c
c
       write(6,*) '                                                    '
c+++++++++++++++++++++
c
c
c Now, do some output
c
      write(6,*) '                                                     '
      write(6,*) '       Output ONE plane of data for each subgrid for '
      write(6,*) '       two-dimensional representation. Each subgrid  '
      write(6,*) '       plane will be written to one file, >GRIDS2D<. '
      write(6,*) '       This file will have multiblock 2-D PLOT3D     '
      write(6,*) '       format.                                       '
      write(6,*) '                                                     '
c
       call out1planetoplt3d_break(
     1 nzone,ioverlap,nodeid,numprocs,
     2 j_subgrid(1),k_subgrid(1),l_subgrid(1),
     3 idimj,idimk,idiml,jdim,kdim,ldim,
     4 ippzone,jklmax,
     5 xorig(1),yorig(1),zorig(1),
     6 jindex_global(1),kindex_global(1),lindex_global(1)
     7 )
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
5     format(' >> Number of processors to be used: ',i6)                
25    format(' For zone ',i3,' use ',i3,' processors.')                 
40    format(//,' Parallel processing parameters:',//,                  
     1'    Number of processors to be used: ',i9,/)                     
42    format(/,' Number of grid points in zone ',i3,': ',i9)            
43    format(/' Total number of grid points over all grids: ',i9,/,     
     1'     Average number of grid points per grid: ',f9.1,/,           
     2'Average number of grid points per processor: ',f9.1,/)           
80    format(/,' User-specified processorss to be used = ',i5,/,        
     1' Nodes calculated to be used = ',i5,/)                           
100   format(//,' Is this ok? If not, you will be asked to enter',/     
     110x,'a different number of processors. Enter (y/n): ')          
200   format(' Zone = ',i3,', Total grid points = ',i9)                 
205   format(/,' Warning: ',i3,' grids will not fit optimally on ',     
     1i5,' processor(s)!',/)                                            
210   format(' Optimal grid points per processor = ',f15.1,/)           
295   format(/,' Zone',' Subzone',' Prc#','  Jmax','  Kmax','  Lmax',   
     1'  Points','   J-rel range','   K-rel range','   L-rel range')    
301   format(//,'==========================================='           
     1,/,' For mesh ',a,':',/)                                          
303   format(/,' Subgrids written to file GRID_DIVIDED_3D.G in PLOT3D fo
     1rmat.',/)                                                         
304   format(///,' Writing subgrids to file.',//)                       
333   format(i5,i6,i6,i6,i6,i6,i7,5x,i3,' to ',i3,4x,i3,' to ',i3,      
     14x,i3,' to ',i3)                                                  
433   format(i5,i6,i6,i6,i6,i6,i7,5x,i3,' to ',i3,4x,i3,' to ',i3,      
     14x,i3,' to ',i3)                                                  
c660   format('     Enter choice (1-5): ',$)                             
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c
c======================================================================|
      subroutine link_overlap(ioverlap,ipnum,nodeid,numprocs,
     1 jdim,kdim,ldim,idimj,idimk,idiml,
     2 j_subgrid,k_subgrid,l_subgrid)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c       This subroutine links overlapped (patched) subgrids in a zone.
c       This subroutine is not called for overlapped regions for which
c       interpolation coefficients have been generated. The subgrids
c       are assumed to align point-to-point.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c       D. W. Barnette
c       Sandia National Laboratories
c       Org. 9221, MS 1111
c       Albuquerque, NM 87112
c       email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c       chooser
c
c Called routines, in order of appearance:
c       write_unravel
c       patch_2d
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c       02-13-96        Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      dimension ipnum(numprocs)
      dimension j_subgrid(numprocs),k_subgrid(numprocs),
     1 l_subgrid(numprocs)

      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.
c
      if (DEBUG) then
       print*,''
       print*, ' >>> In preproc.f/link_overlap <<<'
       print*,''
       print*, ' ipnum(1), ipnum(2) = ',ipnum(1),ipnum(2)
       print*, ' idimj,idimk,idiml = ',idimj,idimk,idiml
       print*, ' numprocs = ',numprocs
       print*, ' ioverlap = ',ioverlap
      endif
c
c link overlapped portions of single grids
c
c Notes:
c       - Grids assumed to consist of right-handed coord. system
c       - Grids in this section have point-to-point matchup
c
c
      jklmax = jdim*kdim*ldim 

      write(6,*) '                                                     '
      write(6,*) ' >> Calculating LINK Table for adjacent sub-grids in e
     1ach zone.'
      write(6,*) '       Output file will be >2D_PATCH_TABLE<, which   '
      write(6,*) '        lists which subgrids are linked to which     '
      write(6,*) '        other neighboring subgrids.                  '
c
c
c write to file >3D_PATCH_TABLE<
c
c OUTPUT: open the >3D_PATCH_TABLE< output file used to store link info
      open(21,file='3D_PATCH_TABLE',iostat=ierr_unit21,form='formatted',
     1 status='unknown')
c     1 status='unknown')
      if(ierr_unit21.ne.0) then
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine LINK_OVERLAP                     '
       write(6,*) '   File 3D_PATCH_TABLE could not be opened.         '
       write(6,*) '   PROBABLE REASON: File exists and cannot be       '
       write(6,*) '                    overwritten.                    '
       write(6,*) '   FIX: If file exists, rename or delete old file.  '
       write(6,*) '                                                    '
       write(6,*) '                                                    '
       write(6,*) '   PROGRAM STOPPED: Problem with i/o file.          '
       write(6,*) '                                                    '
       stop '1: Subroutine LINK_OVERLAP '
      endif
c
c write header for standard output
      write(6,*) '                                                     '
      write(6,*) ' Print first few lines of file >3D_PATCH_TABLE<      '
      write(6,*) '                                                     '
c
c write heading
      write(6,540)
      write(21,540)
c
c start counter
      icnt = 0
c
c cycle thru cube indices (1 cube = 1 subgrid)
      do 20 j=1,idimj
      do 20 k=1,idimk
      do 20 l=1,idiml
c
      index    =(j-1)*idimk*idiml+(k-1)*idiml+l
      index_jm1=(j-2)*idimk*idiml+(k-1)*idiml+l
      index_jp1=(j  )*idimk*idiml+(k-1)*idiml+l
      index_km1=(j-1)*idimk*idiml+(k-2)*idiml+l
      index_kp1=(j-1)*idimk*idiml+(k  )*idiml+l
      index_lm1=(j-1)*idimk*idiml+(k-1)*idiml+l-1
      index_lp1=(j-1)*idimk*idiml+(k-1)*idiml+l+1
c
c do j faces
c
c at j=1
      if(j.eq.1) then
       if(j.ne.idimj) then
        jtarget_min=j_subgrid(index)
        jtarget_max=jtarget_min
        ktarget_min=1
        ktarget_max=k_subgrid(index)
        ltarget_min=1
        ltarget_max=l_subgrid(index)
        jbase_min=1+ioverlap
        jbase_max=1+ioverlap
        kbase_min=1
        kbase_max=k_subgrid(index_jp1)
        lbase_min=1
        lbase_max=l_subgrid(index_jp1)
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(21,15) ipnum(index_jp1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
        if(icnt.le.10) then
         write(6,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
         write(6,15) ipnum(index_jp1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
        endif
       endif
c at j=jmax
      else if(j.eq.idimj)then
       if(j.ne.1) then
        jtarget_min=1
        jtarget_max=1
        ktarget_min=1
        ktarget_max=k_subgrid(index)
        ltarget_min=1
        ltarget_max=l_subgrid(index)
        jbase_min=j_subgrid(index_jm1)-ioverlap
        jbase_max=jbase_min
        kbase_min=1
        kbase_max=k_subgrid(index_jm1)
        lbase_min=1
        lbase_max=l_subgrid(index_jm1)
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(21,15) ipnum(index_jm1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(6,15) ipnum(index_jm1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       endif
       endif
c everything in between j=1 and j=jmax
       else if((j.gt.1).and.(j.lt.idimj)) then
        jtarget_min=j_subgrid(index)
        jtarget_max=jtarget_min
        ktarget_min=1
        ktarget_max=k_subgrid(index)
        ltarget_min=1
        ltarget_max=l_subgrid(index)
        jbase_min=1+ioverlap
        jbase_max=1+ioverlap
        kbase_min=1
        kbase_max=k_subgrid(index_jp1)
        lbase_min=1
        lbase_max=l_subgrid(index_jp1)
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(21,15) ipnum(index_jp1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(6,15) ipnum(index_jp1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       endif
        jtarget_min=1
        jtarget_max=1
        ktarget_min=1
        ktarget_max=k_subgrid(index)
        ltarget_min=1
        ltarget_max=l_subgrid(index)
        jbase_min=j_subgrid(index_jm1)-ioverlap
        jbase_max=jbase_min
        kbase_min=1
        kbase_max=k_subgrid(index_jm1)
        lbase_min=1
        lbase_max=l_subgrid(index_jm1)
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(21,15) ipnum(index_jm1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
c
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(6,15) ipnum(index_jm1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       endif
      endif
c
c        
c do k faces
c
c at k=1
      if(k.eq.1) then
       if(k.ne.idimk) then
        jtarget_min=1
        jtarget_max=j_subgrid(index)
        ktarget_min=k_subgrid(index)
        ktarget_max=ktarget_min
        ltarget_min=1
        ltarget_max=l_subgrid(index)
        jbase_min=1
        jbase_max=j_subgrid(index_kp1)
        kbase_min=1+ioverlap
        kbase_max=1+ioverlap
        lbase_min=1
        lbase_max=l_subgrid(index_kp1)
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(21,15) ipnum(index_kp1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(6,15) ipnum(index_kp1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       endif
       endif
c at k=kmax
      else if(k.eq.idimk)then
       if(k.ne.1) then
        jtarget_min=1
        jtarget_max=j_subgrid(index)
        ktarget_min=1
        ktarget_max=1
        ltarget_min=1
        ltarget_max=l_subgrid(index)
        jbase_min=1
        jbase_max=j_subgrid(index_km1)
        kbase_min=k_subgrid(index_km1)-ioverlap
        kbase_max=kbase_min
        lbase_min=1
        lbase_max=l_subgrid(index_km1)
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(21,15) ipnum(index_km1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(6,15) ipnum(index_km1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       endif
       endif
c everything in between k=1 and k=kmax
       else if((k.gt.1).and.(k.lt.idimk)) then
        jtarget_min=1
        jtarget_max=j_subgrid(index)
        ktarget_min=k_subgrid(index)
        ktarget_max=ktarget_min
        ltarget_min=1
        ltarget_max=l_subgrid(index)
        jbase_min=1
        jbase_max=j_subgrid(index_kp1)
        kbase_min=1+ioverlap
        kbase_max=kbase_min
        lbase_min=1
        lbase_max=l_subgrid(index_kp1)
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(21,15) ipnum(index_kp1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(6,15) ipnum(index_kp1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       endif
        jtarget_min=1
        jtarget_max=j_subgrid(index)
        ktarget_min=1
        ktarget_max=ktarget_min
        ltarget_min=1
        ltarget_max=l_subgrid(index)
        jbase_min=1
        jbase_max=j_subgrid(index_km1)
        kbase_min=k_subgrid(index_km1)-ioverlap
        kbase_max=kbase_min
        lbase_min=1
        lbase_max=l_subgrid(index_km1)
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(21,15) ipnum(index_km1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
c
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(6,15) ipnum(index_km1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       endif
      endif
c
c        
c do l faces
c
c at l=1
      if(l.eq.1) then
       if(l.ne.idiml) then
        jtarget_min=1
        jtarget_max=j_subgrid(index)
        ktarget_min=1
        ktarget_max=k_subgrid(index)
        ltarget_min=l_subgrid(index)
        ltarget_max=ltarget_min
        jbase_min=1
        jbase_max=j_subgrid(index_lp1)
        kbase_min=1
        kbase_max=k_subgrid(index_lp1)
        lbase_min=1+ioverlap
        lbase_max=lbase_min
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(21,15) ipnum(index_lp1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(6,15) ipnum(index_lp1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       endif
       endif
c at l=lmax
      else if(l.eq.idiml)then
       if(l.ne.1) then
        jtarget_min=1
        jtarget_max=j_subgrid(index)
        ktarget_min=1
        ktarget_max=k_subgrid(index)
        ltarget_min=1
        ltarget_max=ltarget_min
        jbase_min=1
        jbase_max=j_subgrid(index_lm1)
        kbase_min=1
        kbase_max=k_subgrid(index_lm1)
        lbase_min=1+ioverlap
        lbase_max=lbase_min
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(21,15) ipnum(index_lm1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
      if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(6,15) ipnum(index_lm1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
      endif
      endif
c everything in between l=1 and l=lmax
       else if((l.gt.1).and.(l.lt.idiml)) then
        jtarget_min=1
        jtarget_max=j_subgrid(index)
        ktarget_min=1
        ktarget_max=k_subgrid(index)
        ltarget_min=l_subgrid(index)
        ltarget_max=ltarget_min
        jbase_min=1
        jbase_max=j_subgrid(index_lp1)
        kbase_min=1
        kbase_max=k_subgrid(index_lp1)
        lbase_min=1+ioverlap
        lbase_max=lbase_min
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(21,15) ipnum(index_lp1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(6,15) ipnum(index_lp1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       endif
        jtarget_min=1
        jtarget_max=j_subgrid(index)
        ktarget_min=1
        ktarget_max=k_subgrid(index)
        ltarget_min=1
        ltarget_max=ltarget_min
        jbase_min=1
        jbase_max=j_subgrid(index_lm1)
        kbase_min=1
        kbase_max=k_subgrid(index_lm1)
        lbase_min=l_subgrid(index_lm1)-ioverlap
        lbase_max=lbase_min
        icnt=icnt+1
        write(21,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(21,15) ipnum(index_lm1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
c
       if(icnt.le.10) then
        write(6,15) ipnum(index),jtarget_min,jtarget_max,1,
     1                               ktarget_min,ktarget_max,1,
     2                                ltarget_min,ltarget_max,1
        write(6,15) ipnum(index_lm1),jbase_min,jbase_max,1,
     1                                 kbase_min,kbase_max,1,
     2                                  lbase_min,lbase_max,1
       endif
      endif
c
20    continue
c
c        
      write(6,*) '                                                     '
      write(6,*) '  Link table finished. File >3D_PATCH_TABLE< written.'
      write(6,*) '  Number of points written to file: ',icnt
      write(6,*) '   This file corresponds to file bcpzn.dat for 3D.   '

c
c 3D_PATCH_TABLE file will now be reduced to the 2D file bcpzn.dat 
       call patch_2d
c
c close 3D_PATCH_TABLE file
       close(21)
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
15    format(2x,i4,2x,i5,8(1x,i5))                                      
540   format('c This file corresponds to file bcpzn.dat for 3D codes.   
     1',/,                                                              
     2'c Zone#  J_beg J_end J_inc K_beg K_end K_inc L_beg L_end L_inc') 
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end

c======================================================================|
      subroutine patch_2d
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c       Reduces the table of 3D patched, i.e. point-to-point overlapped
c       grids, to a 2D table for 2D flow solvers.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c       D. W. Barnette
c       Sandia National Laboratories
c       Org. 9221, MS 1111
c       Albuquerque, NM 87112
c       email: dwbarne@sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c       link_overlap
c
c Called routines, in order of appearance:
c       NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c       02-13-96        Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
c      include "PARAMETERS_DB"
c      parameter(jklmax=jmax*kmax*lmax)
c
      dimension idum(10)
c
      character*1 icolumn1,jkl
      character*80 line
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.
c
c
c write to file >2D_PATCH_TABLE<
c read in 1 line at a time, delete the requested index and list out.
c as a final option, ask to switch indices.
c
c
c OUTPUT: open the >bcpzn.dat< output file used to store link info
      open(22,file='bcpzn.dat',form='formatted',iostat=ierr_unit22,
     1status='unknown')
      if(ierr_unit22.ne.0) then
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine PATCH_2D                         '
       write(6,*) '   File bcpzn.dat could not be opened.         '
       write(6,*) '   PROBABLE REASON: File exists and cannot be       '
       write(6,*) '                    overwritten.                    '
       write(6,*) '   FIX: If file exists, rename or delete old file.  '
       write(6,*) '                                                    '
       write(6,*) '                                                    '
       write(6,*) '   PROGRAM STOPPED: Problem with i/o file.          '
       write(6,*) '                                                    '
       stop ' stop: preproc.f/patch_2d - cannot open file bcpzn.dat    '
      endif
c
      rewind(21)
      read(21,'(a1)',err=12,end=13) icolumn1
      read(21,'(a1)',err=12,end=13) icolumn1
      goto 14
c
12    continue
      write(6,*) '                                                     '
      write(6,*) ' WARNING: Subroutine PATCH_2D                        '
      write(6,*) '   Problem reading file >3D_PATCH_TABLE<.            '
      write(6,*) '   Construction of file >bcpzn.dat< terminated.      '
      write(6,*) '   Program will stop.                                '
      write(6,*) '                                                     '
      close(22)
      stop ' stop: cannot construct bcpzn.dat file.'
c
13    continue
      write(6,*) '                                                     '
      write(6,*) ' WARNING: Subroutine PATCH_2D                        '
      write(6,*) '   Premature end-of-file encountered in file         '
      write(6,*) '   >3D_PATCH_TABLE<. May not be able to find file.   '
      write(6,*) '   Construction of file >bcpzn.dat< terminated.      '
      write(6,*) '   Program will stop.                                '
      write(6,*) '                                                     '
      close(22)
      stop ' stop: cannot construct bcpzn.dat file.'
c
14    continue
c
      jkl = 'k'
      write(6,*) ' Deleting index=',jkl
      write(6,*) '                                                    '
c
c write heading
      write(22,5)
c
c start a counter for the number of lines on the file
      icount=0
c
30     continue
       read(21,*,end=40,err=25) (idum(i),i=1,10)
       write(22,60) idum(1),(idum(i),i=2,4),(idum(i),i=8,10)
       icount=icount+1
       goto 30
25     continue
       write(6,*) ' '
       write(6,*) ' WARNING: subroutine patch_2d'
       write(6,*) '   Problem reading scratch file >3D_PATCH_TABLE<.'
       write(6,*) '   Construction of file >bcmain.dat< terminated.'
       write(6,*) '   Program will terminate.'
       write(6,*) ' '
       call exit(1)
40     continue
c
      write(6,*) '                                                     '
      write(6,*) '  File >bcpzn.dat< written.                          '
      write(6,*) '  Total lines written to file >bcpzn.dat<:',icount
      write(6,*) '                                                     '
c
c print contents of file >bcpzn.dat<
500   continue

      if (DEBUG) then
       write(6,*) '                                                    '
       write(6,*) ' Contents of file >bcpzn.dat<                       '
       write(6,*) '                                                    '
       rewind(22)

c skip over first three lines
       read(22,'(a80)') line
       write(6,'(a80)') line
       read(22,'(a80)') line
       write(6,'(a80)') line
       read(22,'(a80)') line
       write(6,'(a80)') line

       icnt=0
100    continue
c       icnt = icnt + 1
       read(22,*,end=110) (idum(i),i=1,7)
       if(icnt.le.10) write(6,60) (idum(i),i=1,7)
       read(22,*) (idum(i),i=1,7)
       if(icnt.le.10) write(6,60) (idum(i),i=1,7)
       goto 100
c
110    continue
      endif

      close(22)
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
5     format('c This file corresponds to file bcpzn.dat for 2D codes.   
     1',/,                                                              
     2'c first line is target (recipient); next line is base (donor)',/,
     3'c Zone#  J_beg J_end J_inc K_beg K_end K_inc')                   
60    format(2x,i4,2x,i5,5(1x,i5))                                      
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end

c======================================================================|
      subroutine min_surf_vol_ratio
     1 (nz,ippzone,jmaxx,kmaxx,lmaxx,jpts_per_subgrid,kpts_per_subgrid,
     2  lpts_per_subgrid,nzone,nodeid,numprocs,idimj,idimk,idiml,
     3  jdim,kdim,ldim)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c 	This subroutine calculates the 3-factors for a given number of
c   	subgrids, then calculates the surface to volume ratio for each 
c	set of 3-factors given the maximum j,k,l dimensions of each grid. 
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c	break
c	breakup_for_fun
c
c Called routines, in order of appearance:
c	NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
      parameter (imax=300)
c      
      dimension idim1(imax),idim2(imax),idim3(imax)  
c
      jklmax=jdim*kdim*ldim
      max_subs = numprocs
c
      write(6,*) ' List all possible factors for number of subgrids    '
      write(6,*) '       for this zone.                                '
      write(6,5)
      fnumprocs=real(numprocs)     
      itag=0
      do 10 i=1,numprocs
       xtry1=fnumprocs/float(i)
       check1=xtry1-int(xtry1)
       if(check1.eq.0) then   
        do 15 j=1,int(xtry1)
         xtry2=xtry1/float(j)
         check2=xtry2-int(xtry2)
         if(check2.eq.0) then
          do 20 k=1,int(xtry2)
           xtry3=xtry2/float(k)
           check3=xtry3-int(xtry3)
           check_answer=i*j*k
           if((check3.eq.0.0).and.(check_answer.eq.numprocs)) then
            itag=itag+1
            if(itag.gt.imax) then
             write(6,*) '                                              '
             write(6,*) ' ERROR: Subroutine MIN_SURF_VOL_RATIO         '
             write(6,*) '    ITAG is greater than IMAX!                '
             write(6,*) '      ITAG = ',itag
             write(6,*) '      IMAX = ',imax
             write(6,*) '        Change parameter statement!           '
             write(6,*) '                                              '
             write(6,*) '  PROGRAM STOPPED.                            '
             write(6,*) '                                              '
             stop 'stopping - 1: Subroutine MIN_SURF_VOL_RATIO '
            endif
            idim1(itag)=i
            idim2(itag)=j
            idim3(itag)=k
            write(6,25) itag,i,j,k
           endif
20        continue
         endif
15      continue
       endif
10    continue    
c
c calculate the min(surf/vol)
      write(6,*) '                                                     '
      write(6,*) '                                                     '
      write(6,*) ' Finding min(surf/vol) for:                          '
      write(6,*) '      Jmax = ',jmaxx
      write(6,*) '      Kmax = ',kmaxx
      write(6,*) '      Lmax = ',lmaxx    
c 
c write header for output
      write(6,30)
      do 35 i=1,itag
      surf_over_vol = 2*(float(idim1(i))/float(jmaxx)
     1  + float(idim2(i))/float(kmaxx) + float(idim3(i))/float(lmaxx)) 
      if(i.gt.1) then
       xmin_save=xmin_surf_over_vol
       xmin_surf_over_vol = amin1(xmin_surf_over_vol, surf_over_vol) 
       if (xmin_surf_over_vol.ne.xmin_save) then 
        imin=i
       endif    
      else
       xmin_surf_over_vol = surf_over_vol
       imin=1    
      endif  
      x1=surf_over_vol
      x2=xmin_surf_over_vol
      write(6,40) i,x1,x2
35    continue 
      idimj=idim1(imin)
      idimk=idim2(imin)
      idiml=idim3(imin)
c
c error check
      idim_check=idimj*idimk*idiml
      if(idim_check.gt.max_subs) then
       imax_zone=0
       do 200 i=1,nzone
        imax_zone=max(ippzone,imax_zone)
200    continue
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine MIN_SURF_VOL_RATIO               '
       write(6,*) '   More subgrids are needed than allowed for in     '
       write(6,*) '   parameter statement.                             '
       write(6,100) nz,idim_check,max_subs,imax_zone,imax_zone
       write(6,*) '                                                    '
       write(6,*) '   PROGRAM STOPPED.                                 '
       write(6,*) '                                                    '
       stop 'stopping - 2: Subroutine MIN_SURF_VOL_RATIO '
      endif
c
c output grid values
c   jpts_per_subgrid, etc: approx number of grid points per subgrid

      jpts_per_subgrid=(jmaxx-1)/idimj + 1
      kpts_per_subgrid=(kmaxx-1)/idimk + 1
      lpts_per_subgrid=(lmaxx-1)/idiml + 1
      write(6,*) '                                                     '
      write(6,*) ' Final values for min(surf/vol) for this zone are:   '
      write(6,*) '    itag = ',imin
      write(6,50) xmin_surf_over_vol
      write(6,60) idimj,jmaxx,jpts_per_subgrid
      write(6,70) idimk,kmaxx,kpts_per_subgrid
      write(6,80) idiml,lmaxx,lpts_per_subgrid
      write(6,*) '                                                     '
c
c  
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
5     format(/,' ITAG',7x,'J',7x,'K',7x,'L',/                           
     1         ' ----',6x,'---',5x,'---',5x,'---')                      
25    format(i5,6x,i3,5x,i3,5x,i3)                                      
30    format(/,' itag',7x,'surf/vol',5x,'min(surf/vol)')                
40    format(i4,5x,f11.6,5x,f11.6)                                      
50    format(5x,'min(surf/vol) = ',f11.6)                               
60    format(5x,'N/Jmax = ',i4,' /',i4,' ==> approx. ',i3,' pts/subgrid 
     1in J direction')                                                  
70    format(5x,'M/Kmax = ',i4,' /',i4,' ==> approx. ',i3,' pts/subgrid 
     1in K direction')                                                  
80    format(5x,'P/Lmax = ',i4,' /',i4,' ==> approx. ',i3,' pts/subgrid 
     1in L direction')                                                  
100   format(/,                                                         
     1'              Subgrids needed for zone ',i5,': ',i5,/            
     2'       Number of subgrids currently allowed : ',i5,/             
     3'              Max number of subgrids needed : ',i5,//            
     3'   Change max_subs in preproc.f/min_surf_vol_ratio to ',i5/)     
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'min_surf_vol_ratio'


c======================================================================|
      subroutine subgrid_dimensions(nz,jmaxx,kmaxx,lmaxx,
     1 jpts_per_subgrid,kpts_per_subgrid,lpts_per_subgrid,
     2 nodeid,numprocs,idimj,idimk,idiml,j_subgrid,k_subgrid,l_subgrid)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c 	This subroutine calculates the dimensions of the subgrids.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c	break
c
c Called routines, in order of appearance:
c	NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c      
      dimension j_subgrid(numprocs),k_subgrid(numprocs),
     & l_subgrid(numprocs)
c
c check if grids are evenly divisible or if remainder is left over
      j_remainder = mod(jmaxx-1,idimj)
      k_remainder = mod(kmaxx-1,idimk)
      l_remainder = mod(lmaxx-1,idiml)
      if(j_remainder.ne.0) then
       write(6,*) ' Warning: Jmaxx/N not evenly divisible!'
       write(6,5) j_remainder 
      endif          
      if(k_remainder.ne.0) then
       write(6,*) ' Warning: Kmaxx/M not evenly divisible!'
       write(6,10) k_remainder
      endif          
      if(l_remainder.ne.0) then
       write(6,*) ' Warning: Lmaxx/P not evenly divisible!'
       write(6,15) l_remainder
      endif 
c
      if((j_remainder.ne.0).or.(k_remainder.ne.0).or.
     1 (l_remainder.ne.0)) then
       write(6,*) '                                                    '
       write(6,*) '  BREAKUP will make adjustments in the size of the  '
       write(6,*) '    subgrids to compensate for the above remainders.'
       write(6,*) '    Adjustments will start at the boundaries and    '
       write(6,*) '    progress inward. Hence, when remainders appear, '
       write(6,*) '    outer subgrids will have their dimensions       '
       write(6,*) '    increased in the corresponding direction.       '
      endif
c 
c fill in approximate values      
      do 20 j=1,idimj
      do 20 k=1,idimk
      do 20 l=1,idiml 
       index=(j-1)*idimk*idiml+(k-1)*idiml+l
       j_subgrid(index)=jpts_per_subgrid
       k_subgrid(index)=kpts_per_subgrid
       l_subgrid(index)=lpts_per_subgrid
20    continue
c now include remainders
      j_pool=j_remainder
      k_pool=k_remainder
      l_pool=l_remainder
c work on j values 
      if(j_remainder.ne.0) then
       do 30 k=1,idimk
       do 30 l=1,idiml
        index_j1=(1-1)*idimk*idiml+(k-1)*idiml+l
        j_subgrid(index_j1)=j_subgrid(index_j1)+1  
30     continue
       j_pool=j_pool-1
       if(j_pool.ne.0) then
        do 35 k=1,idimk
        do 35 l=1,idiml
         index_jmax=(idimj-1)*idimk*idiml
     1   +(k-1)*idiml+l
         j_subgrid(index_jmax)=j_subgrid(index_jmax)+1
35      continue
        j_pool=j_pool-1
        if(j_pool.ne.0) then
         jmaxcp1=idimj+1
         jmaxcp1_half=jmaxcp1/2
         do 40 j=2,jmaxcp1_half
          do 45 k=1,idimk
          do 45 l=1,idiml
           index=(j-1)*idimk*idiml+(k-1)*idiml+l
           j_subgrid(index)=j_subgrid(index)+1
45        continue
          j_pool=j_pool-1
          if(j_pool.eq.0) goto 55
          jtemp=idimj-j+1
          do 50 k=1,idimk
          do 50 l=1,idiml
           index_jtemp=(jtemp-1)*idimk*idiml
     1     +(k-1)*idiml+l
           j_subgrid(index_jtemp)=j_subgrid(index_jtemp)+1  
50        continue
          j_pool=j_pool-1
          if(j_pool.eq.0) goto 55 
40       continue 
55       continue
        endif
       endif
      endif 
c      
c work on k values
      if(k_remainder.ne.0) then
       do 60 j=1,idimj
       do 60 l=1,idiml
        index_k1=(j-1)*idimk*idiml+(1-1)*idiml+l
        k_subgrid(index_k1)=k_subgrid(index_k1)+1  
60     continue
       k_pool=k_pool-1
       if(k_pool.ne.0) then
        do 65 j=1,idimj
        do 65 l=1,idiml
         index_kmax=(j-1)*idimk*idiml
     1   +(idimk-1)*idiml+l
         k_subgrid(index_kmax)=k_subgrid(index_kmax)+1
65      continue
        k_pool=k_pool-1
        if(k_pool.ne.0) then
         kmaxcp1=idimk+1
         kmaxcp1_half=kmaxcp1/2
         do 70 k=2,kmaxcp1_half
          do 75 j=1,idimj
          do 75 l=1,idiml
           index=(j-1)*idimk*idiml+(k-1)*idiml+l
           k_subgrid(index)=k_subgrid(index)+1
75        continue
          k_pool=k_pool-1
          if(k_pool.eq.0) goto 85
          ktemp=idimk-k+1
          do 80 j=1,idimj
          do 80 l=1,idiml
           index_ktemp=(j-1)*idimk*idiml
     1     +(ktemp-1)*idiml+l
           k_subgrid(index_ktemp)=k_subgrid(index_ktemp)+1  
80        continue
          k_pool=k_pool-1
          if(k_pool.eq.0) goto 85 
70       continue 
85       continue
        endif
       endif
      endif 
c
c work on l values
      if(l_remainder.ne.0) then
       do 90 j=1,idimj
       do 90 k=1,idimk
        index_l1=(j-1)*idimk*idiml+(k-1)*idiml+1
        l_subgrid(index_l1)=l_subgrid(index_l1)+1  
90     continue
       l_pool=l_pool-1
       if(l_pool.ne.0) then
        do 95 j=1,idimj
        do 95 k=1,idimk
         index_lmax=(j-1)*idimk*idiml
     1   +(k-1)*idiml+idiml
         l_subgrid(index_lmax)=l_subgrid(index_lmax)+1
95      continue
        l_pool=l_pool-1
        if(l_pool.ne.0) then
         lmaxcp1=idiml+1
         lmaxcp1_half=lmaxcp1/2
         do 100 l=2,lmaxcp1_half
          do 105 j=1,idimj
          do 105 k=1,idimk
           index=(j-1)*idimk*idiml+(k-1)*idiml+l
           l_subgrid(index)=l_subgrid(index)+1
105       continue
          l_pool=l_pool-1
          if(l_pool.eq.0) goto 115
          ltemp=idiml-l+1
          do 110 j=1,idimj
          do 110 k=1,idimk
           index_ltemp=(j-1)*idimk*idiml
     1     +(k-1)*idiml+ltemp
           l_subgrid(index_ltemp)=l_subgrid(index_ltemp)+1  
110       continue
          l_pool=l_pool-1
          if(l_pool.eq.0) goto 115 
100      continue 
115      continue
        endif
       endif
      endif 
c      
c write out final dimensions for subgrids 
      write(6,*) '                                                     '
      write(6,*) ' Final subgrid dimensions are:                       '
c
c print out final grid point values 
      write(6,120)
      inum=0
      do 125 j=1,idimj
      do 125 k=1,idimk
      do 125 l=1,idiml
       index=(j-1)*idimk*idiml+(k-1)*idiml+l
       inum=inum+1
       itotl=j_subgrid(index)*k_subgrid(index)
     1 *l_subgrid(index)
       write(6,130) inum,j,k,l,j_subgrid(index),k_subgrid(index),
     1 l_subgrid(index),itotl
125   continue
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
5     format(5x,'J-remainder = ',i3)                                    
10    format(5x,'K-remainder = ',i3)                                    
15    format(5x,'L-remainder = ',i3)                                    
120   format(/,' Grid #',4x,'Jcube',4x,'Kcube',4x,'Lcube',4x,'Jmax',4x, 
     1 'Kmax',4x,'Lmax',4x,'Total',/,                                   
     2         '-------',4x,'-----',4x,'-----',4x,'-----',4x,'----',4x, 
     3 '----',4x,'----',3x,'-------')                                   
130   format(i6,5x,i4,5x,i4,5x,i4,5x,i4,4x,i4,4x,i4,2x,i8)              
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'subgrid_dimensions'


c======================================================================|
      subroutine grid_point_comparison(nzone,jmaxx,kmaxx,lmaxx,
     1 nodeid,numprocs,idimj,idimk,idiml,j_subgrid,k_subgrid,l_subgrid)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c 	This subroutine calculates the total number of points for the 
c	original grid, and the total number of grid points for the 
c	subgrids of the original grid, and compares the two. The program 
c	stops if there is a difference. This is to provide some error 
c	checking.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c	break
c	breakup_for_fun
c
c Called routines, in order of appearance:
c	NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c      
c      
      dimension j_subgrid(numprocs),k_subgrid(numprocs),
     1 l_subgrid(numprocs)
      logical DEBUG

      DEBUG = .false.
c      DEBUG = .true.
c 
c calculate total number of grid points and compare with original
      write(6,*) '                                                     '
      write(6,*) ' --- ERROR CHECK ---'
      write(6,*) ' Compare total number of grid points with original   '
      write(6,*) '     in this zone:                                   '
      ipoints_orig=jmaxx*kmaxx*lmaxx 

      if (DEBUG) then
        print*, ' In preproc.f/grid_point_comparison:'
        print*, ' idimj = ',idimj
        print*, ' idimk = ',idimk
        print*, ' idiml = ',idiml
        print*, '         i    j_subgrid   k_subgrid    l_subgrid'
        print*, '    ------    ---------   ---------    ---------'
        do 5 i=1,numprocs
         print*,i,j_subgrid(i),k_subgrid(i),l_subgrid(i)
5       continue
      endif
c      
c First, add all subgrids together
      ipoints_calc=0
      do 10 j=1,idimj
      do 10 k=1,idimk
      do 10 l=1,idiml
       index=(j-1)*idimk*idiml+(k-1)*idiml+l
       ipoints_calc=ipoints_calc+(j_subgrid(index)
     1  *k_subgrid(index)*l_subgrid(index))
10    continue
c
c Now, subtract common faces
c   j-direction 
      jflag=0
      do 20 j=1,idimj-1
      jflag=1
      do 20 k=1,idimk
      do 20 l=1,idiml
       index=(j-1)*idimk*idiml+(k-1)*idiml+l
       ipoints_calc=ipoints_calc-k_subgrid(index)
     1 *l_subgrid(index)
20    continue 
      if(jflag.eq.1) then 
       write(6,*) '   Subtracted common faces in J-direction.          '
      endif
c   k-direction
      kflag=0 
      do 30 k=1,idimk-1
      kflag=1
      do 30 j=1,idimj
      do 30 l=1,idiml
       index=(j-1)*idimk*idiml+(k-1)*idiml+l
       if(j.ne.idimj) then 
        ipoints_calc=ipoints_calc
     1  -(j_subgrid(index)-jflag)*l_subgrid(index)
       else
        ipoints_calc=ipoints_calc
     1  -j_subgrid(index)*l_subgrid(index)
       endif
30    continue
      if(kflag.eq.1) then 
       write(6,*) '   Subtracted common faces in K-direction.          '
      endif
c   l-direction
      lflag=0
      do 35 l=1,idiml-1
      lflag=1
      do 35 j=1,idimj
      do 35 k=1,idimk
       index=(j-1)*idimk*idiml+(k-1)*idiml+l
       if((j.eq.idimj).and.(k.eq.idimk)) then
        ipoints_calc=ipoints_calc
     1  -j_subgrid(index)*k_subgrid(index)
       else if(j.eq.idimj) then
        ipoints_calc=ipoints_calc
     1  -j_subgrid(index)*(k_subgrid(index)-kflag)
       else if(k.eq.idimk) then
        ipoints_calc=ipoints_calc
     1  -(j_subgrid(index)-jflag)*k_subgrid(index)
       else 
        ipoints_calc=ipoints_calc
     1  -(j_subgrid(index)-jflag)*(k_subgrid(index)-kflag)
       endif
35    continue
      if(lflag.eq.1) then 
       write(6,*) '   Subtracted common faces in L-direction.          '
      endif
c      
      write(6,40) ipoints_orig,ipoints_calc
c     
      if(ipoints_orig.ne.ipoints_calc) then
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine GRID_POINT_COMPARISON            '
       write(6,*) '   Original number of grid points does not          '
       write(6,*) '   equal the number summed from the subgrids.       '
       write(6,*) '                                                    '
       write(6,*) '  PROGRAM STOPPED.                                  '
       write(6,*) '                                                    '
       stop '1: Subroutine GRID_POINT_COMPARISON '
      endif
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
40    format(5x,'...Original no. of points in this zone = ',i10,/,      
     15x,'...Calculated no. of points in this zone = ',i10)             
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'grid_point_comparison'


c======================================================================|
      subroutine global_indices(nzone,meshname,nodeid,numprocs,
     & j_subgrid,k_subgrid,l_subgrid,idimj,idimk,idiml,
     2 jindex_global,kindex_global,lindex_global)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c 	This subroutine sets up the global indices needed to break up 
c	the zone. Call this routine after calling subroutine 
c	subgrid_dimensions.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c	break
c
c Called routines, in order of appearance:
c	NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      dimension j_subgrid(numprocs),k_subgrid(numprocs),
     & l_subgrid(numprocs)
      dimension jindex_global(numprocs),kindex_global(numprocs),
     1 lindex_global(numprocs)

c
      integer idimj,idimk,idiml
      character*10 meshname
c
      nz = nzone
      max_subs = numprocs
c
      write(6,6)
c
c do checks on dimensions
      idim_max=0
c
      write(6,12)
c
c write original zone dimensions
       write(6,13) nz,idimj,idimk,idiml
       idim_max=max(idim_max,idimj*idimk*idiml)
c
      if(idim_max.gt.max_subs) then
       write(6,1) idim_max,max_subs
       write(6,14) max_subs,idim_max
       write(6,*) '                                                    '
       write(6,*) '  PROGRAM STOPPED.                                  '
       write(6,*) '                                                    '
       stop '1: Subroutine GLOBAL_INDICES '
      endif
c
c
c
       do 10 j=1,idimj
       do 10 k=1,idimk
       do 10 l=1,idiml
         index=(j-1)*idimk*idiml+(k-1)*idiml+l
         index_lm1=(j-1)*idimk*idiml+(k-1)*idiml+l-1
         if(l.eq.1) then
         lindex_global(index)=l_subgrid(index)
        else
         lindex_global(index)=l_subgrid(index)
     1                   +lindex_global(index_lm1)-1
        endif
10     continue                                               
c                          
       do 15 j=1,idimj
       do 15 l=1,idiml
       do 15 k=1,idimk
         index=(j-1)*idimk*idiml+(k-1)*idiml+l
         index_km1=(j-1)*idimk*idiml+(k-2)*idiml+l
        if(k.eq.1) then
         kindex_global(index)=k_subgrid(index)
        else
         kindex_global(index)=k_subgrid(index)
     1                   +kindex_global(index_km1)-1
        endif
15     continue                                               
c                          
       do 20 k=1,idimk
       do 20 l=1,idiml
       do 20 j=1,idimj
       index=(j-1)*idimk*idiml+(k-1)*idiml+l
       index_jm1=(j-2)*idimk*idiml+(k-1)*idiml+l
        if(j.eq.1) then
         jindex_global(index)=j_subgrid(index)
        else
         jindex_global(index)=j_subgrid(index)
     1                   +jindex_global(index_jm1)-1
        endif
20     continue 
c
       write(6,*) '                                                    '
       write(6,*) '                                                    '
       write(6,21) nzone,meshname
       write(6,*)  ' Grd# Jcube Kcube Lcube Jmax Kmax Lmax Jming Jmaxg K
     1ming Kmaxg Lming Lmaxg'
       write(6,*)  ' ---- ----- ----- ----- ---- ---- ---- ----- ----- -
     1---- ----- ----- -----'
c      
c j,k,l are cube dimensions in original grid
c 
       icount=0      
       do 25 j=1,idimj
       do 25 k=1,idimk
       do 25 l=1,idiml
        index=(j-1)*idimk*idiml+(k-1)*idiml+l
        icount=icount+1
        jp3d_max=jindex_global(index)
        kp3d_max=kindex_global(index)
        lp3d_max=lindex_global(index)
        jp3d_min=jp3d_max-j_subgrid(index)+1
        kp3d_min=kp3d_max-k_subgrid(index)+1 
        lp3d_min=lp3d_max-l_subgrid(index)+1 
        jlocal=j_subgrid(index)
        klocal=k_subgrid(index)
        llocal=l_subgrid(index)
        write(6,26) icount,j,k,l,jlocal,klocal,llocal,jp3d_min,jp3d_max,
     1  kp3d_min,kp3d_max,lp3d_min,lp3d_max
25     continue
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
1     format(/,' ERROR: Subroutine GLOBAL_INDICES',/                    
     1'     Array dimensions exceeded',/                                
     2'         idim_max = ',i4,/                                       
     3'         max_subs = ',i4,/                                       
     4'     Adjust parameter statement for max_subs.')                  
6     format(//,' -----------------------------------------',/          
     1' Local and Global Indices for Each Subgrid',/                    
     2' -----------------------------------------')                     
12    format(//' Summary of zone dimensions:',//                        
     1'  Zone   jcube_max  kcube_max  lcube_max',/                      
     2' ------  ---------  ---------  ---------')                       
13    format(i6,i10,2i11)                                               
14    format(//,' > Parameter statement needs to be updated:',/         
     1'       Current values: max_subs = ',i5,/                         
     2' Suggest changing to : max_subs = ',i5,/)                        
21    format(' Zone =',i4,', Meshname = ',a)                            
26    format(i5,4i6,3i5,2i6,3i6)                                        
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'global_indices'


c======================================================================|
      subroutine write_subgrids(nzone,iunit,form,nodeid,numprocs,
     1 ioverlap,idimj,idimk,idiml,jdim,kdim,ldim,
     2 j_subgrid,k_subgrid,l_subgrid,
     3 ippzone,jklmax,
     4 x,y,z,
     5 jindex_global,kindex_global,lindex_global)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c 	This subroutine writes out the subgrids formed in subroutine
c 	SUBGRID_DIMENSIONS and subroutine GLOBAL_INDICES.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c	break
c	breakup_for_fun
c
c Called routines, in order of appearance:
c	get_global_index
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
#include "forttype.h"

      dimension j_subgrid(numprocs),k_subgrid(numprocs),
     1 l_subgrid(numprocs)
      dimension jindex_global(numprocs),kindex_global(numprocs),
     1  lindex_global(numprocs)
c
      __REAL x(jklmax),y(jklmax),z(jklmax)
c
      logical form
      logical DEBUG
c
      DEBUG = .false.
c      DEBUG = .true.
c
      if (DEBUG) then
        print*,''
        print*, ' >>> In preproc.f/write_subgrids <<<'
        print*, '  ioverlap = ',ioverlap
        print*, '  jklmax = ',jklmax
        print*,''
      endif
c
c      jklmax=jdim*kdim*ldim
      jmax_compout=jdim
      kmax_compout=kdim
      lmax_compout=ldim

c User decides whether output is formatted or unformatted
c
      write(6,*) ' '
c
       open(75,file='GRID_DIVIDED_3D.G',form='formatted',
     1   iostat=ierr_unit75,status='unknown')
      if(ierr_unit75.ne.0) then
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine WRITE_SUBGRIDS                   '
       write(6,*) '   File GRID_DIVIDED_3D.G could not be opened.      '
       write(6,*) '                                                    '
       write(6,*) '                                                    '
       write(6,*) '   PROGRAM STOPPED: Problem with i/o file.          '
       write(6,*) '                                                    '
       stop '1: Subroutine WRITE_SUBGRIDS '
      endif
c
      write(6,*) 'Output file GRID_DIVIDED_3D.G opened as formatted.' 
      write(6,*) '                                                  '
c
c
c write out plot file for subgrids in PLOT3D format
c
c write header

c
      if (DEBUG) then
        print*,''
        print*, ' number_of_processors = ',numprocs
        print*, ' ioverlap = ',ioverlap
        print*, ' jklmax = ',jklmax
        print*, ' jmax_compout = ',jmax_compout
        print*, ' kmax_compout = ',kmax_compout
        print*, ' lmax_compout = ',lmax_compout
        print*, ' idimj,idimk,idiml = ',idimj,idimk,idiml
        print*, ' nzone = ',nzone
      endif

      write(75,*) numprocs
      print*,'j_subgrid = ',j_subgrid
      print*,'k_subgrid = ',k_subgrid
      print*,'l_subgrid = ',l_subgrid
c      stop 'stopped'
c
c      do 5 nz=1,nzone
       write(75,8) 
     1(((j_subgrid((j-1)*idimk*idiml+(k-1)*idiml+l),
     2   k_subgrid((j-1)*idimk*idiml+(k-1)*idiml+l),
     3   l_subgrid((j-1)*idimk*idiml+(k-1)*idiml+l),
     4                l=1,idiml),k=1,idimk),j=1,idimj) 
5      continue
c
c
c Extract each subgrid from the original PLOT3D file
c
      write(6,700) iunit
c
      i_total_count=0
c
c      do 30 nz=1,nzone
      nz = 1
c
      jmaxc=jmax_compout
      kmaxc=kmax_compout
      lmaxc=lmax_compout

      write(6,*) ' '
      write(6,*) ' jmaxc and jmax_compout = ',jmaxc
      write(6,*) ' kmaxc and kmax_compout = ',kmaxc
      write(6,*) ' lmaxc and lmax_compout = ',lmaxc
      write(6,*) ' '
c
c read zone grid points from original grid file
c
      write(6,800) nz
      write(6,802)
c
c       read(iunit,*,end=1000,err=1500) 
c     1    (((x((jm-1)*kmaxc*lmaxc+(km-1)*lmaxc+lm),
c     2         jm=1,jmaxc),km=1,kmaxc),lm=1,lmaxc),
c     3    (((y((jm-1)*kmaxc*lmaxc+(km-1)*lmaxc+lm),
c     4         jm=1,jmaxc),km=1,kmaxc),lm=1,lmaxc),
c     5    (((z((jm-1)*kmaxc*lmaxc+(km-1)*lmaxc+lm),
c     6         jm=1,jmaxc),km=1,kmaxc),lm=1,lmaxc)
c
c
c      goto 1200
c
c error check: end of grid file reached prematurely
c
c1000  continue
c      write(6,*) '                                                     '
c      write(6,*) ' ERROR: Subroutine WRITE_SUBGRIDS                    '
c      write(6,*) '   End of input grid file reached prematurely        '
c      write(6,*) '   while trying to read xyz data.                    '
c      write(6,*) '                                                     '
c      write(6,*) '   PROGRAM STOPPED.                                  '
c      write(6,*) '                                                     '
c      stop '4: Subroutine WRITE_SUBGRIDS '
c
c error check: error reading file
c
c1500  continue
c      write(6,*) '                                                     '
c      write(6,*) ' ERROR: Subroutine WRITE_SUBGRIDS                    '
c      write(6,*) '   Not able to read xyz data of input grid file.     '
c      write(6,*) '                                                     '
c      write(6,*) '   PROGRAM STOPPED.                                  '
c      write(6,*) '                                                     '
c      stop '5: Subroutine WRITE_SUBGRIDS '
c
1200  continue
c
c construct subgrids for the zone just read in
c
      itotal=idimj*idimk*idiml
c
      print*,''
      print*, ' Writing to GRID_DIVIDED_3D.G:'
      print*,''
c
c loop over sub-blocks of grids
      num_grid_points = 0

      do 30 j=1,idimj
      do 30 k=1,idimk
      do 30 l=1,idiml
c
       i_total_count=i_total_count+1
c
       index_ref=(j-1)*idimk*idiml + (k-1)*idiml + l
c
       write(6,705) index_ref,itotal,nz,i_total_count,numprocs,ioverlap
c
c
      print*,''
      print*, ' >>> Calling get_global_index'
      print*,''

       call get_global_index
     1 (nz,j,k,l,1,1,1,
     2 jp3d_min,jp3d_max,kp3d_min,kp3d_max,lp3d_min,lp3d_max,
     3 ioverlap,nodeid,numprocs,
     4 j_subgrid(1),k_subgrid(1),l_subgrid(1),
     5 idimj,idimk,idiml,jdim,kdim,ldim,
     6 jindex_global(1),kindex_global(1),lindex_global(1))
c
       write(6,900) nz,j,k,l,jp3d_min,jp3d_max,kp3d_min,
     1 kp3d_max,lp3d_min,lp3d_max
       write(6,901) jp3d_max - jp3d_min + 1, kp3d_max - kp3d_min + 1,
     1   lp3d_max - lp3d_min + 1
c
c
c formatted
       write(75,1900) 
     1      (((x((lp-1)*kmaxc*jmaxc + (kp-1)*jmaxc + jp),
     2                         jp=jp3d_min,jp3d_max),
     3                         kp=kp3d_min,kp3d_max),
     4                         lp=lp3d_min,lp3d_max),
     5      (((y((lp-1)*kmaxc*jmaxc + (kp-1)*jmaxc + jp),
     6                         jp=jp3d_min,jp3d_max),
     7                         kp=kp3d_min,kp3d_max),
     8                         lp=lp3d_min,lp3d_max),
     9      (((z((lp-1)*kmaxc*jmaxc + (kp-1)*jmaxc + jp),
     A                         jp=jp3d_min,jp3d_max),
     B                         kp=kp3d_min,kp3d_max),
     C                         lp=lp3d_min,lp3d_max)
c     1      (((x((jp-1)*kmaxc*lmaxc+(kp-1)*lmaxc+lp),
c     2                         jp=jp3d_min,jp3d_max),
c     3                         kp=kp3d_min,kp3d_max),
c     4                         lp=lp3d_min,lp3d_max),
c     5      (((y((jp-1)*kmaxc*lmaxc+(kp-1)*lmaxc+lp),
c     6                         jp=jp3d_min,jp3d_max),
c     7                         kp=kp3d_min,kp3d_max),
c     8                         lp=lp3d_min,lp3d_max),
c     9      (((z((jp-1)*kmaxc*lmaxc+(kp-1)*lmaxc+lp),
c     A                         jp=jp3d_min,jp3d_max),
c     B                         kp=kp3d_min,kp3d_max),
c     C                         lp=lp3d_min,lp3d_max)
c
c
      num_grid_points = num_grid_points + 
     1  (jp3d_max - jp3d_min + 1)*(kp3d_max - kp3d_min + 1)*
     2   (lp3d_max - lp3d_min + 1)*3

30    continue  


      print*,''
      print*, ' Total number of grid points in GRID_DIVIDED_3D.G: ',
     1 num_grid_points
      print*,''
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
8     format(3i10)                                                      
700   format(//,' Constructing subgrids. reading from unit ',i3,/)      
705   format(                                                           
     1' Constructing subgrid #',i5,' of',i5,' subgrids in zone ',i4,/   
     2'   Total count: subgrid #',i4,' of ',i4,' subgrids',/            
     3'   Grid cells overlap',i2,' cell widths on grid boundaries.')    
800   format(/,                                                         
     1         ' ---------------------------------------------------',/ 
     2         ' Reading original grid file, zone #',i5)                
802   format(  '  This grid does not contain IBLANK data.',/            
     1         ' ---------------------------------------------------',/)
900   format(                                                           
     1'   nz = ',i4,' jcube = ',i4,' kcube = ',i4,' lcube = ',i4,/      
     2'       jp3d_min = ',i4,'  jp3d_max = ',i4,/                      
     3'       kp3d_min = ',i4,'  kp3d_max = ',i4,/                      
     4'       lp3d_min = ',i4,'  lp3d_max = ',i4)                       
901   format(
     1'Local block dimensions: ',i4,' x ',i4,' x ',i4,/
     2)
1900  format(5e19.11)
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'write_subgrids'   


c======================================================================|
      subroutine get_global_index
     1 (nz,j,k,l,jdo,kdo,ldo,j1,j2,k1,k2,l1,l2,ioverlap,
     2 nodeid,numprocs,j_subgrid,k_subgid,l_subgrid,
     3 idimj,idimk,idiml,jdim,kdim,ldim,
     4 jindex_global,kindex_global,lindex_global)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c 	This subroutine is called to give the index range for subgrids
c  	using global indices. It is called from various locations
c  	throughout the code.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c	break
c	link_overset
c	out1planetoplt3d_break
c	out1planetoproc_break
c	outallplanestoplt3d_break
c	outallplanestoproc_break
c	write_subgrids
c
c Called routines, in order of appearance:
c	NONE
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   	Arguments	Description
c 	---------	-----------------------------------------
c Input:
c	nz		zone number
c	j		index for global index array
c	k		index for global index array
c	l		index for global index array
c	jdo		if 1, determine j1 and j2 indices
c	kdo		if 1, determine k1 and k2 indices
c	ldo		if 1, determine l1 and l2 indices
c       ioverlap        extent of grid cell overlap between subgrids
c Output (global values):
c	j1		min j index value for subgrid j,k,l in zone nz
c	j2		max j index value for subgrid j,k,l in zone nz
c	k1		min k index value for subgrid j,k,l in zone nz
c	k2		max k index value for subgrid j,k,l in zone nz
c	l1		min l index value for subgrid j,k,l in zone nz
c	l2		max l index value for subgrid j,k,l in zone nz
c
c
      dimension j_subgrid(numprocs),k_subgrid(numprocs),
     1 l_subgrid(numprocs)

      dimension jindex_global(numprocs),kindex_global(numprocs),
     1 lindex_global(numprocs)

      integer idimj,idimk,idiml,jdim,kdim,ldim,ioverlap

      logical DEBUG
c
      DEBUG = .false.
c      DEBUG = .true.
c
      jmax=jdim
      kmax=kdim
      lmax=ldim
      jklmax=jmax*kmax*lmax

c check value for ioverlap
      if((ioverlap.ne.0).and.(ioverlap.ne.1).and.(ioverlap.ne.2)) then
       write(3,*) '                                                    '
       write(3,*) ' ERROR: Subroutine GET_GLOBAL_INDEX                 '
       write(3,*) '   ioverlap in sub. get_global_index is not valid   '
       write(3,10) ioverlap
       write(3,*) '                                                    '
       write(3,*) '  PROGRAM STOPPED.                                  '
       write(3,*) '                                                    '
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine GET_GLOBAL_INDEX                 '
       write(6,*) '   ioverlap in sub. get_global_index is not valid   '
       write(6,10) ioverlap
       write(6,*) '                                                    '
       write(6,*) '  PROGRAM STOPPED.                                  '
       write(6,*) '                                                    '
       stop '1: Subroutine GET_GLOBAL_INDEX '
      endif
c
c set indices
      index    =(j-1)*idimk*idiml+(k-1)*idiml+l
      index_jm1=(j-2)*idimk*idiml+(k-1)*idiml+l
      index_km1=(j-1)*idimk*idiml+(k-2)*idiml+l
      index_lm1=(j-1)*idimk*idiml+(k-1)*idiml+l-1

c
      if (DEBUG) then
        print*,''
        print*, ' In sub preproc.f/get_global_index:'
        print*, '   numprocs = ',numprocs
        print*, '   index = ',index
        print*, '   index_jm1 = ',index_jm1
        print*, '   index_km1 = ',index_km1
        print*, '   index_lm1 = ',index_lm1
        print*, '   ioverlap = ',ioverlap
        print*, '   j,k,l = ',j,k,l
        print*, '   idimj,idimk,idiml = ',idimj,idimk,idiml
        print*, '   jdim,kdim,ldim = ',jdim,kdim,ldim
        print*, '   jdo,kdo,ldo = ',jdo,kdo,ldo
        print*,''
        print*, 'Max global indices for this subgrid:'
        print*, index,jindex_global(index),kindex_global(index),
     1    lindex_global(index)
c        do 555 i=1,numprocs
c           print*,i,jindex_global(i),kindex_global(i),lindex_global(i)
c555     continue
      endif
c
      if(jdo.ne.0) then
c j values
       if (idimj.eq.1) then
        j2=jindex_global(index)
        j1=1
       else if(j.eq.1) then
        j2=jindex_global(index)+ioverlap
        j1=1
       else if(j.eq.idimj) then
        j2=jindex_global(index)
c        j1=j2-j_subgrid(index)+1
        j1=jindex_global(index_jm1)
       else
        j2=jindex_global(index)+ioverlap
c        j1=j2-j_subgrid(index)+1
        j1=jindex_global(index_jm1)
       endif
      endif
c
      if(kdo.ne.0) then
c k values
       if (idimk.eq.1) then
        k2=kindex_global(index)
        k1=1
       else if(k.eq.1) then
        k2=kindex_global(index)+ioverlap
        k1=1
       else if(k.eq.idimk) then
        k2=kindex_global(index)
c        k1=k2-k_subgrid(index)+1 
        k1=kindex_global(index_km1)
       else
        k2=kindex_global(index)+ioverlap
c        k1=k2-k_subgrid(index)+1
        k1=kindex_global(index_km1)
       endif
      endif
c
      if(ldo.ne.0) then
c l values
       if (idiml.eq.1) then
        l2=lindex_global(index)
        l1=1
       else if(l.eq.1) then
        l2=lindex_global(index)+ioverlap
        l1=1
       else if(l.eq.idiml) then
        l2=lindex_global(index)
c        l1=l2-l_subgrid(index)+1
        l1=lindex_global(index_lm1)
       else
        l2=lindex_global(index)+ioverlap
c        l1=l2-l_subgrid(index)+1
        l1=lindex_global(index_lm1)
       endif
      endif
c
      if (DEBUG) then
        print*, ' j1,k1,l1 = ',j1,k1,l1
        print*, ' j2,k2,l2 = ',j2,k2,l2
        print*,''
      endif

c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
10    format('      ioverlap = ',i10)                                   
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'get_global_index'

c======================================================================|
       subroutine out1planetoplt3d_break( 
     1   nzone,ioverlap,nodeid,numprocs,
     2   j_subgrid,k_subgrid,l_subgrid,
     3   idimj,idimk,idiml,jdim,kdim,ldim,
     4   ippzone,jklmax,
     5   x,y,z,
     6   jindex_global,kindex_global,lindex_global
     7   )
c
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Purpose:
c   	Write one plane of each 
c	subgrid to file GRIDS2D (UNIT 20) in PLOT3D format for graphics.
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Author: 
c	D. W. Barnette
c   	Sandia National Laboratories
c   	Org. 9221, MS 1111
c     	Albuquerque, NM 87112
c  	email: dwbarne@cs.sandia.gov
c++++++++++++++++++++++++++++++++++++++++++++++++++
c Calling routines:
c	break
c
c Called routines, in order of appearance:
c	read_zone_header
c	read_zone
c	get_global_index
c	get_global_index
c	get_global_index
c	
c++++++++++++++++++++++++++++++++++++++++++++++++++
c History:
c	02-13-96	Reference (start) date
c++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   Guess at these; error catching will inform user if not correct
c      parameter(iicntr_max=maxdim*500,
c     1 ibcntr_max=iicntr_max)
c
#include "forttype.h"

      __REAL x(jklmax),y(jklmax),z(jklmax)
      dimension j_subgrid(numprocs),k_subgrid(numprocs),
     1  l_subgrid(numprocs)
      dimension jindex_global(numprocs),kindex_global(numprocs),
     1  lindex_global(numprocs)
c
      character*1 jkl
      character*2 idirect
      logical DEBUG
      character*10 filename

      DEBUG = .false.
c      DEBUG = .true.
c
      jmax_compout=jdim
      kmax_compout=kdim
      lmax_compout=ldim
c
      write(6,*) '                                                     '
      write(6,*) ' >> WRITING 2-D GRIDS TO FILE >GRIDS2D< IN PLOT3D FORM
     1AT'
c
      write(6,*) ' Writing one plane of data for each grid in PLOT3D for
     1mat.'
c
c
c OUTPUT: open the GRIDS2D output file used to store all grids in PLOT3D 
c  format
      filename='GRIDS2D'
      open(20,file=filename,iostat=ierr_unit20,form='formatted',
     1 status='unknown')
      if(ierr_unit20.ne.0) then
       write(6,*) '                                                    '
       write(6,*) ' ERROR: Subroutine OUT1PLANETOPLT3D_BREAK           '
       write(6,*) '   File GRIDS2D could not be opened.                '
       write(6,*) '                                                    '
       write(6,*) '   PROGRAM STOPPED: Problem with i/o file.          '
       write(6,*) '                                                    '
       stop '1: Subroutine OUT1PLANETOPLT3D_BREAK '
      endif
c
c calculate total number of processors
      iproc=0
      do 5 i=1,nzone
       iproc=iproc+ippzone
5     continue
      write(6,*) ' Writing data for a total # of processors = ',iproc
      write(6,*) '                                                     '
c      
c set the error flag to zero
      ierr=0
c
1     continue
c      write(6,10)
c      read(5,'(a1)') jkl
      jkl='k'
c
      write(6,*) ' Output one plane normal to direction = ',jkl
      write(6,*) '                                                     '
c
c If K
       kplane=1
c check on kplane bounds
       do 30 j=1,idimj
       do 30 k=1,idimk
       do 30 l=1,idiml
        index=(j-1)*idimk*idiml+(k-1)*idiml+l
        if((1.le.kplane).and.
     1  (k_subgrid(index).ge.kplane)) then
         goto 30
        else
         ierr=1
        endif
30     continue
       if(ierr.eq.1) then
        write(6,*) '                                                   '
        write(6,*) ' The >k< value is not within bounds. Try again.    '
        write(6,*) '                                                   '
        return
       endif
c
c
c output data
c
50    continue
c
       k1=kplane
       k2=kplane
c
c output xz data
140   continue
      write(6,*)' '
c
      idirect='xz'
c
      write(6,*) ' Output grid data lying in the plane of ',idirect
      write(6,*) '                                                     '
c
      write(20,*) iproc
c
       if (DEBUG) then
        print*, ' idimj,idimk,idiml = ',idimj,idimk,idiml
       endif
       do 310 j=1,idimj
       do 310 k=1,idimk
       do 310 l=1,idiml
        index=(j-1)*idimk*idiml + (k-1)*idiml + l
        if (DEBUG) then
         print*, ' index,j,k,l = ',index,j,k,l
        endif
        write(20,*) j_subgrid(index),l_subgrid(index)
310    continue

c ERROR CHECK
c rewind grid file and read past header
      rewind(20)
      read(20,*) idum_zones
      if (DEBUG) then
       print*, ' reading idum_zones from 20 = ',idum_zones
      endif

      if (idum_zones.ne.numprocs) then
        write(6,*) ' '
        write(6,*) ' idum_zones = ',idum_zones
        write(6,*) ' numprocs = ',numprocs
        write(6,*) ' '
        write(6,*) '     PROGRAM STOPPED: idum_zones.ne.numprocs'
        write(6,*) ' '
        stop '1: error in file header: idum_zones.ne.numprocs'
      endif
      read(20,*) (idum1,idum2,nz=1,numprocs)

c
      irunning_count=0
c
c
      jmaxc=jmax_compout
      kmaxc=kmax_compout
      lmaxc=lmax_compout
c
      do 100 j=1,idimj
      do 100 k=1,idimk
      do 100 l=1,idiml
c
      index=(j-1)*idimk*idiml+(k-1)*idiml+l
c
      irunning_count=irunning_count+1
c
        call get_global_index
     1  (i,j,k,l,1,0,1,j1,j2,k1,k2,l1,l2,ioverlap,
     2  nodeid,numprocs,
     3  j_subgrid(1),k_subgrid(1),l_subgrid(1),
     4  idimj,idimk,idiml,jdim,kdim,ldim,
     5  jindex_global(1),kindex_global(1),lindex_global(1))
c
      write(6,600) irunning_count,iproc,i,index,jmaxc,kmaxc,lmaxc,
     1 j,k,l,ioverlap,j1,j2,k1,k2,l1,l2
c
c
        write(20,*)
c     1        (((x((jd-1)*kmaxc*lmaxc+(kd-1)*lmaxc+ld),
     1        (((x((ld-1)*kmaxc*jmaxc + (kd-1)*jmaxc + jd),
     2                   jd=j1,j2),kd=k1,k2),ld=l1,l2),
c     3        (((z((jd-1)*kmaxc*lmaxc+(kd-1)*lmaxc+ld),
     3        (((z((ld-1)*kmaxc*jmaxc + (kd-1)*jmaxc + jd),
     4                   jd=j1,j2),kd=k1,k2),ld=l1,l2)
c
100   continue
c
      write(6,*) '                                                     '
      write(6,*) ' >> 2-D GRIDS HAVE BEEN WRITTEN TO FILE >GRIDS2D< IN P
     1LOT3D FORMAT'
      write(6,*) '                                                     '
c
      close(20)
c
c
c format statements
c
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
600   format(/,                                                         
     1' ===========================================' ,/                 
     2' Writing the following data to file GRIDS2D',/                   
     3' Subgrid #',i4,' of ',i4,' total subgrids.',/                    
     4' Zone # :',i4,'     Subgrid index = ',i4,/                       
     5' File format: blanked formatted 2-D multiblock PLOT3D',/         
     6'     jmax_compout = ',i4,/                                       
     7'     kmax_compout = ',i4,/                                       
     8'     lmax_compout = ',i4,/                                       
     9'     subgrid location : j = ',i4,', k = ',i4,', l = ',i4,/       
     A'     ioverlap = ',i1,/                                           
     B'     dimensions: j = ',i4,'  to ',i4,/                           
     C'                 k = ',i4,'  to ',i4,/                           
     D'                 l = ',i4,'  to ',i4,/                           
     E' ===========================================')                   
c---x----1----x----2----x----3----x----4----x----5----x----6----x----7--|
c
c
      return
      end
c
c end 'out1planetoplt3d_break'
