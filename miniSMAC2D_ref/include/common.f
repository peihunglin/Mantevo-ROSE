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
c====================================================================
c     
c
#include "precis.h"
c only one zone per processor, so nzne = 1
      parameter(nzne=1,ibcmax=500,jkmax=150,nwallmax=20)
      common/original/jmax_o,kmax_o
      common/airfr/alpha,clift,cdrag,cmom
      dimension dfp(jkmax,3), dfm(jkmax,3)
      dimension f(jkmax,3),fv(jkmax,3)
      common/dflx/dfp, dfm
      common/flxs/f,fv
      common/iparmi/nt,niter,ntmax,ntime,iflxo,ivis,iturb
      common/iparmr/cdis,epscon
      common/impli/impsch,njsp,nksp,nbcimp,ntjsp,ntksp
      common/implr/underr,dcoef2,dcoef4
      common/mass/pin,pout
      common/peri/kend2, kendm
      common/rmssi/jres,kres
      common/rmssr/resmax0,resmax,divmax,turres
      common/timecpu/cputime,cputiter,cpuitpt
      common/timm/beta,dtau,dt,time
      common/visc/reynum,vnu
      dimension ibcval(ibcmax), nzbc(ibcmax), jkbc(ibcmax),
     & jbcb(ibcmax), jbce(ibcmax),
     & kbcb(ibcmax), kbce(ibcmax),
     & nzwall(ibcmax),
     & jkwall(ibcmax), jkinc(ibcmax),
     & jwall1(ibcmax), jwall2(ibcmax),
     & kwall1(ibcmax), kwall2(ibcmax),
     & xwallval(jkmax),ywallval(jkmax),
     & xwallv(jkmax,nwallmax),ywallv(jkmax,nwallmax)
      common/bcsurf/nbcreg, ibcval, nzbc, jkbc,
     & jbcb, jbce,
     & kbcb, kbce,
     & nwall, nzwall,
     & jkwall, jkinc,
     & jwall1, jwall2,
     & kwall1, kwall2,
     & xwallval,ywallval,
     & xwallv,ywallv
      dimension nzbc_t_wake(ibcmax),nzbc_b_wake(ibcmax),
     & jbcb_t_wake(ibcmax),jbce_t_wake(ibcmax),jinc_t_wake(ibcmax),
     & jbcb_b_wake(ibcmax),jbce_b_wake(ibcmax),jinc_b_wake(ibcmax)
      common/bcwake/nbcreg_wake,nzbc_t_wake,nzbc_b_wake,
     & jbcb_t_wake,jbce_t_wake,jinc_t_wake,
     & jbcb_b_wake,jbce_b_wake,jinc_b_wake
      dimension nzt(ibcmax), nzb(ibcmax),
     & jbt(ibcmax), jet(ibcmax), jit(ibcmax),
     & kbt(ibcmax), ket(ibcmax), kit(ibcmax),
     & jbb(ibcmax), jeb(ibcmax), jib(ibcmax),
     & kbb(ibcmax), keb(ibcmax), kib(ibcmax)
      common/pzne/nreg, nzt, nzb,
     & jbt, jet, jit,
     & kbt, ket, kit,
     & jbb, jeb, jib,
     & kbb, keb, kib
      dimension nzint(ibcmax), nintrp(ibcmax),
     & iintbeg(ibcmax), iintend(ibcmax),
     & iint(ibcmax,6)
      common/intpi/nintreg, nzint, nintrp,
     & iintbeg, iintend,
     & iint
      common/unitno/istdout
      common/yaml/title,geometry
      character*60 title,geometry
      double precision cputime,cputiter,cpuitpt
