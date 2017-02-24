cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  rrcov : Scalable Robust Estimators with High Breakdown Point
cc
cc  This program is free software; you can redistribute it and/or modify
cc  it under the terms of the GNU General Public License as published by
cc  the Free Software Foundation; either version 2 of the License, or
cc  (at your option) any later version.
cc
cc  This program is distributed in the hope that it will be useful,
cc  but WITHOUT ANY WARRANTY; without even the implied warranty of
cc  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
cc  GNU General Public License for more details.
cc
cc  You should have received a copy of the GNU General Public License
cc  along with this program; if not, write to the Free Software
cc  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
cc
cc  I would like to thank Peter Rousseeuw and Katrien van Driessen for
cc  providing the initial code of this function.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rfltsreg(dat,n,nvar,nhalff,krep, inbest,objfct,
     *     intercept,intadjust,nvad,datt,
     *     weights,temp,index1,index2,aw2,aw,residu,y,nmahad,ndist,
     *     am,am2,slutn,
     *     jmiss,xmed,xmad,a,da,h,hvec,c,cstock,mstock,c1stock,
     *     m1stock,dath,sd,means,bmeans, i_trace)
c
c    dat   = cbind(x,y)   hence  n x (p+1)
c    nvar  = p
c    nvad  = p+1
c    nhalff: 'quan' = quan.f = h.alpha.n(alpha, n, rk)  which is
c                   = (n + p + 1) %/% 2  when alpha= 1/2
c    krep  = nsamp  (e.g. = 5000 for "best")
c
      implicit none
      integer kmini, nmini, k1,k2,k3,km10,nmaxi,maxmini
      integer i_trace
      integer i_aux(4) ! just for printing when i_trace > 0
c
ccc   parameter (nvmax=115)
ccc   parameter (nmax=57000)
cc
      parameter (kmini=5)
      parameter (nmini=300)
      parameter (k1=2)
      parameter (k2=2)
      parameter (k3=100)

cc  krep := the total number of trial subsamples
cc          to be drawn when n exceeds 2*nmini;
c           krep = 0  :<==>  "exact"  <==>  all possible subsamples

cc
ccc   parameter (nvmax1=nvmax+1)
ccc   parameter (nvmax2=nvmax*nvmax)
ccc   parameter (nvm11=nvmax*(nvmax+1))

      parameter (km10=10*kmini)
      parameter (nmaxi=nmini*kmini)
C--   VT   parameter (maxmini=int((3*nmini-1)/2)+1)
      parameter (maxmini=450)
cc
      integer n,nvar,nvad,nhalff

      integer inbest(nhalff)
      double precision dat(n,nvad)
      double precision datt(n,nvad)

      double precision weights(n)
      integer temp(n)
      integer index1(n)
      integer index2(n)
      double precision aw2(n),aw(n)
      double precision residu(n)
      double precision y(n)
      double precision nmahad(n)
      double precision ndist(n)
      double precision am(n),am2(n),slutn(n)

      integer krep, matz,iseed, seed, tottimes,step
      integer intercept,intadjust
      integer pnsel, replow
      integer i,ii,iii, j,jj,jjj, jndex, k,kk, lll, m,mm, nn
      integer jmin,jmax, jerd,jnc, jreg, kstep,kount
c unused      integer jdefaul, jbreak
      integer minigr
      integer nfac,nerr, ngroup, nhalf,nlen,nmax,nmore, nmore2, nquant
      integer nvmax1, nvm11, nvmax, nsel, nstop, nrep

      double precision bstd, dist2, eps, factor, objfct, object
      double precision fckw, fckwi, fckw1, percen
      double precision MADeps

      logical all,part,fine,final,rfodd,more1,more2
c unused      integer rfnbreak
      integer rfncomb

      integer flag(km10)
      integer mini(kmini)
      integer subdat(2,nmaxi)
      integer subndex(maxmini)
      double precision faclts(11)
      double precision mcdndex(10,2,kmini)

c     Function
      double precision rffindq

ccc   integer jmiss(nvmax1)
ccc   double precision xmed(nvmax1)
ccc   double precision xmad(nvmax1)
ccc   double precision a(nvmax1), da(nvmax1)
ccc   double precision h(nvmax,nvmax1),hvec(nvm11)
ccc   double precision c(nvmax,nvmax1)
ccc   double precision cstock(10,nvmax2)
ccc   double precision mstock(10,nvmax)
ccc   double precision c1stock(km10,nvmax2)
ccc   double precision m1stock(km10,nvmax)
ccc   double precision dath(nmaxi,nvmax1)
ccc   double precision sd(nvmax)
ccc   double precision means(nvmax)
ccc   double precision bmeans(nvmax)

      integer jmiss(nvad)
      double precision xmed(nvad)
      double precision xmad(nvad)
      double precision a(nvad), da(nvad)
      double precision h(nvar,nvad),hvec(nvar*nvad)
      double precision c(nvar,nvad)
      double precision cstock(10,nvar*nvar)
      double precision mstock(10,nvar)
      double precision c1stock(km10,nvar*nvar)
      double precision m1stock(km10,nvar)
      double precision dath(nmaxi,nvad)
      double precision sd(nvar)
      double precision means(nvar)
      double precision bmeans(nvar)

      data faclts/2.6477,2.5092,2.3826,2.2662,2.1587,
     *     2.0589,1.9660,1.879,1.7973,1.7203,1.6473/

      if(i_trace .ge. 2) then
         call intpr('Entering rfltsreg() - krep: ',-1,krep,1)
      endif

      call rndstart
C     -------- == GetRNGstate() in C

CCCC  10.10.2005 - substitute the parameters nmax and nvmax
      nmax = n
      nvmax = nvar
      nvmax1 = nvmax+1
      nvm11 = nvmax*(nvmax+1)

      nrep = krep

      if(nvar .lt.5 ) then
         eps=1.0D-12
      else
         if(nvar .ge. 5 .and. nvar .le. 8) then
            eps=1.0D-14
         else
            eps=1.0D-16
         endif
      endif
c     Tolerance for rfstatis():  |MAD| < MADeps  : <==> "problem"
      MADeps=1.0D-6

cc	nhalff=int((n+nvar+1)/2)

      jmin=(n/2)+1
      jmax = max((3*n/4)+(nvar+1)/4, nhalff)
      nquant=min(11, 1+ nint(40*(dble(nhalff)/n - 0.5)))
      factor=faclts(nquant)
c unused      jbreak=rfnbreak(nhalff,n,nvar)
c unused      jdefaul=(n+nvar+1)/2
      percen = (1.D0*nhalff)/(1.D0*n)
      if(nvad.eq.1) goto 9000
cc
CDDD  CALL INTPR('>>> Enter RFLTSREG ... iseed=',-1,iseed,1)
      seed=iseed
      matz=1
      nsel=nvar
      ngroup=1
      part=.false.
      fine=.false.
      final=.false.
      all=.true.
      do 21,i=1,nmaxi
        subdat(1,i)=1000000
        subdat(2,i)=1000000
 21   continue
cc
      mini(1)=0
      mini(2)=0
      mini(3)=0
      mini(4)=0
      mini(5)=0
      if(krep.gt.0 .and. n.gt.(2*nmini-1)) then
        kstep=k1
        part=.true.
        ngroup=int(n/dble(nmini))
        if(n.ge.(2*nmini) .and. n.le.(3*nmini-1)) then
          if(rfodd(n)) then
            mini(1)=int(n/2)
            mini(2)=int(n/2)+1
          else
            mini(1)=n/2
            mini(2)=n/2
          endif
        else if(n.ge.(3*nmini) .and. n.le.(4*nmini-1)) then
          if(3*(n/3) .eq. n) then
            mini(1)=n/3
            mini(2)=n/3
            mini(3)=n/3
          else
            mini(1)=int(n/3)
            mini(2)=int(n/3)+1
            if(3*(n/3) .eq. n-1) then
              mini(3)=int(n/3)
            else
              mini(3)=int(n/3)+1
            endif
          endif
        else if(n.ge.(4*nmini) .and. n.le.(5*nmini-1)) then
          if(4*(n/4) .eq. n) then
            mini(1)=n/4
            mini(2)=n/4
            mini(3)=n/4
            mini(4)=n/4
          else
            mini(1)=int(n/4)
            mini(2)=int(n/4)+1
            if(4*(n/4) .eq. n-1) then
              mini(3)=int(n/4)
              mini(4)=int(n/4)
            else
              if(4*(n/4) .eq. n-2) then
                mini(3)=int(n/4)+1
                mini(4)=int(n/4)
              else
                mini(3)=int(n/4)+1
                mini(4)=int(n/4)+1
              endif
            endif
          endif
        else
          mini(1)=nmini
          mini(2)=nmini
          mini(3)=nmini
          mini(4)=nmini
          mini(5)=nmini
        endif

        nhalf=int(mini(1)*percen)
        if(ngroup.gt.kmini) ngroup=kmini
        nrep=int(dble(krep)/ngroup)
        minigr=mini(1)+mini(2)+mini(3)+mini(4)+mini(5)
        if(i_trace .ge. 2) call intpr(' rftls.... minigr=',-1,minigr,1)
        call rfrdraw(subdat,n,minigr,mini,ngroup,kmini)
      else
c          krep == 0  or   n <=  2*nmini-1  ( = 599 by default)

        minigr=n
        nhalf=nhalff
        kstep=k1

C VT::25.11.2010 - added krep==0 means "exact" (all combinations)
        if(krep.eq.0 .or. n.le.replow(nsel)) then
c		use all combinations; happens iff  nsel = nvar = p <= 6
           nrep=rfncomb(nsel,n)
           if(i_trace .ge. 2) then
               call intpr('will use *all* combinations: ',-1,nrep,1)
           endif
        else
          nrep = krep
          all=.false.
        endif
      endif

      seed=iseed
cc

CDDD  CALL INTPR('>>> Start initialization ... nrep=',-1,nrep,1)

      do 31, j=1,nvmax
        do k=1,10
          mstock(k,j)=1000000.D0
          do kk=1,kmini
            m1stock((kk-1)*10+k,j)=1000000.D0
          end do
          do i=1,nvmax
            do kk=1,kmini
              c1stock((kk-1)*10+k,(j-1)*nvmax+i)=1000000.D0
            end do
            cstock(k,(j-1)*nvmax+i)=1000000.D0
          end do
        end do
        means(j)=0.D0
        bmeans(j)=0.D0
        sd(j)=0.D0
        do k=1,nvmax1
          c(j,k)=0.D0
          h(j,k)=0.D0
        end do
 31   continue

      do 41, j=1,nmax
        nmahad(j)=0.D0
        ndist(j)=0.D0
        index1(j)=1000000
        index2(j)=1000000
        temp(j)=1000000
        weights(j)=0.D0
        aw(j)=0.D0
        aw2(j)=0.D0
        residu(j)=0.D0
        y(j)=0.D0
        am(j)=0.D0
        am2(j)=0.D0
        slutn(j)=0.D0
 41   continue

      do j=1,km10
        flag(j)=1
      end do
      do 45, j=1,nvmax1
        jmiss(j)=0
        xmed(j)=0.D0
        xmad(j)=0.D0
        a(j)=0.D0
        da(j)=0.D0
        do k=1,nmaxi
          dath(k,j)=0.D0
        end do
 45   continue

      do j=1,maxmini
        subndex(j)=0
      end do
      do j=1,nvm11
        hvec(j)=0.D0
      end do

      if(i_trace .ge. 2)
     +    call intpr(' rftls.... initialization ready',-1,0,1)
 9000 continue

      if(nvad.eq.1) then
        do jj=1,n
          ndist(jj)=dat(jj,1)
        end do
        call rfshsort(ndist,n)
        call rfmcduni(ndist,n,nhalff,slutn,bstd,am,am2,factor,
     *       n-nhalff+1)
        goto 9999
      endif
cc
      if(.not.fine .and. .not.final) then
        call rfstatis(dat,xmed,xmad,aw2,intercept,nvad, nvmax1,nmax,n,
     *       nstop,MADeps,weights,y,nvar,index2)
        if(nstop.eq.1) goto 9999
      endif

cc
      jreg=1
      call rflsreg(nvmax1, nvmax,nvar,n,a,dat, weights, da, h,
     *     fckw,hvec,nvm11,jmiss,nvad,n)
cc             nfac=nvad-1
      nfac=nvar-1
      call rfrtran(nvar,intercept,nfac,nvad,nvmax1,xmed,
     *     xmad,a,nvad,fckw)
      call rftrc(h,da,nvmax,nvmax1,nvar,intercept,nfac,nvad,
     *     xmed,xmad)
      jerd=0

      tottimes=0
c---- - - - - - - - - Outermost loop - - - - - - - - - - - - - - - - - - -
c----
 5555 object=10.D25

      if(.not. part .or. final)             nn=n
      if(part .and. fine .and. .not. final) nn=minigr

      if(fine.or.(.not.part.and.final)) then
        nrep=10
        nsel=nhalf
        kstep=k2
        if (final) then
          nhalf=nhalff
          ngroup=1
          if (n*nvar .le.100000) then
            kstep=k3
          else if (n*nvar .gt.100000 .and. n*nvar .le.200000) then
            kstep=10
          else if (n*nvar .gt.200000 .and. n*nvar .le.300000) then
            kstep=9
          else if (n*nvar .gt.300000 .and. n*nvar .le.400000) then
            kstep=8
          else if (n*nvar .gt.400000 .and. n*nvar .le.500000) then
            kstep=7
          else if (n*nvar .gt.500000 .and. n*nvar .le.600000) then
            kstep=6
          else if (n*nvar .gt.600000 .and. n*nvar .le.700000) then
            kstep=5
          else if (n*nvar .gt.700000 .and. n*nvar .le.800000) then
            kstep=4
          else if (n*nvar .gt.800000 .and. n*nvar .le.900000) then
            kstep=3
          else if (n*nvar .gt.900000 .and. n*nvar .le.1000000) then
            kstep=2
          else
            kstep=1
          endif

          if (n .gt. 5000) nrep=1
        else
          nhalf= int(minigr*percen)
        endif
      endif

      if(i_trace .ge. 2) then
         i_aux(1) = nrep
         i_aux(2) = kstep
         i_aux(3) = nhalf
         call intpr('Main (number of trials nrep, kstep, nhalf):',
     +        -1, i_aux, 3)
      endif
      do i=1,nsel-1
         index1(i)=i
      end do
      index1(nsel)=nsel-1
cc
      if(.not. final) then
        do i=1,10
          do j=1,ngroup
            mcdndex(i,1,j)=10.D25
            mcdndex(i,2,j)=10.D25
          end do
        end do
      endif
      if (fine .and. .not. final) then
        do j=1,minigr
          do k=1,nvad
            dath(j,k)=dat(subdat(1,j),k)
          end do
        end do
      endif
      kount=0

CDDD  CALL INTPR('>>> MAIN LOOP BY GROUPS: NGROUP= ',-1,ngroup,1)

      do 1111 ii = 1,ngroup
         if(i_trace .ge. 3)
     +        call intpr(' rftls.... looping by group ii=',-1,ii,1)
        if(.not.fine) kount=0
        if(part .and. .not. fine) nn=mini(ii)
        do i=1,nn
          index2(i)=i
        end do

        if(part .and. .not. fine) then
          jndex=0
          do j=1,minigr
            if(subdat(2,j).eq.ii) then
              jndex=jndex+1
              subndex(jndex)=subdat(1,j)
            endif
          end do

          do j=1,mini(ii)
            do k=1,nvad
              dath(j,k)=dat(subndex(j),k)
            end do
          end do
        endif

        do 1000 i=1,nrep
          if(i_trace .ge. 4)
     +          call intpr(' rftls.... for(i = 1,nrep): i=',-1,i,1)
          pnsel=nsel
          tottimes=tottimes+1
          fckwi=0.D0
          fckw1=0.D0
          step=0
 132      if((part.and..not.fine).or.(.not.part.and..not.final)) then
            if(part) then
              call rfrangen(mini(ii),nsel,index1)
            else
              if(all) then
                call rfgenpn(n,nsel,index1)
              else
                call rfrangen(n,nsel,index1)
              endif
            endif
          endif

c 9550     continue
          if(.not.fine .and. part) then
             do j=1,pnsel
                do m=1,nvad
                   c(j,m)=dath(index1(j),m)
                end do
             end do
          else if(.not.part .and. .not.final) then
             do j=1,pnsel
                do m=1,nvad
                   c(j,m)=dat(index1(j),m)
                end do
             end do
          endif

          if((.not.part.and..not.final).or.(.not.fine.and.part)) then
            if(nvar.gt.1) then
              call rfequat(c,nvmax,nvmax1,hvec,nvm11,nvar,1,nerr)
              if(nerr.lt.0) then
                 jerd=jerd+1
                 if(.not.all .and. i.gt.2) goto 132
                 goto 1000
              endif
            else
              if(c(1,1).ne.0.D0) c(1,1)=c(1,2)/c(1,1)
            endif

            do jnc=1,nvar
               a(jnc)=c(jnc,1)
            end do
          endif

          if (final) then
            if(mstock(i,1).ne.1000000.D0) then
              do jj=1,nvar
                a(jj)=mstock(i,jj)
              end do
            else
              goto 1111
            endif
          endif
          if (fine.and..not.final) then
            if(m1stock((ii-1)*10+i,1).ne.1000000.D0) then
              do jj=1,nvar
                a(jj)=m1stock((ii-1)*10+i,jj)
              end do
            else
              goto 1111
            endif
          endif

          do jnc=1,nn
            residu(jnc)=0.D0
            do j=1,nvar
              if(part.and..not.final) then
                residu(jnc)=residu(jnc)+a(j)*dath(jnc,j)
              else
                residu(jnc)=residu(jnc)+a(j)*dat(jnc,j)
              endif
            end do
            if(part.and..not.final) then
              residu(jnc)=dath(jnc,nvad)-residu(jnc)
            else
              residu(jnc)=dat(jnc,nvad)-residu(jnc)
            endif
            aw(jnc)=residu(jnc)
          end do

          more1=.false.
          more2=.false.
          nmore=200
          nmore2=nmore/2

          if(intadjust.eq.1) then

CDDD  CALL INTPR('>>> INTERCEPT ADJUSTMENT 1',-1,i,1)
            if(intercept.eq.1.and.((.not.fine.and.part).or.
     *           .not.part.or.((nn-nhalf).le.nmore))) then
              call rfshsort(aw,nn)
              call rfmcduni(aw,nn,nhalf,slutn,bstd,am,am2,
     *             factor,nn-nhalf+1)
              a(nvar)=a(nvar)+slutn(1)
              do jnc=1,nn
                residu(jnc)=residu(jnc)-slutn(1)
              end do
            else if(intercept.eq.1) then
              call rfshsort(aw,nn)
              do jj=1,nn
                am2(jj)=abs(aw(jj))
              end do
              dist2=rffindq(am2,nn,nhalf,index1)
              do jj=1,nhalf
                aw2(jj)=aw(index1(jj))
              end do
              dist2=rffindq(aw2,nhalf,1,index2)
              jnc=index1(index2(1))
              if(jnc+nmore-nmore2+nhalf-1.gt.nn.or.jnc-nmore2.lt.1)
     *             then
                call rfmcduni(aw,nn,nhalf,slutn,bstd,am,am2,
     *               factor,nn-nhalf+1)
                a(nvar)=a(nvar)+slutn(1)
                do jnc=1,nn
                   residu(jnc)=residu(jnc)-slutn(1)
                end do
              else
 555            do jj=0,nhalf-1+nmore
                   aw2(jj+1)=aw(jnc-nmore2+jj)
                end do
                nlen=nmore+1
                call rfmcduni(aw2,nhalf+nmore,nhalf,slutn,
     *               bstd,am,am2,factor,nlen)
                if(nlen.eq.1.and..not.more1) then
                  if(.not.more2) then
                    nmore=nmore2
                    nmore2=nmore2+nmore2
                    more1=.true.
                    if(jnc-nmore2.ge.1) goto 555
                  endif
                else if(nlen.eq.(nmore+1).and..not.more2) then
                  if(.not.more1) then
                    nmore=nmore2
                    nmore2=-nmore2
                    more2=.true.
                    if(jnc+nmore-nmore2+nhalf-1.le.nn)
     *                   goto 555
                  endif
                else if(nlen.eq.1.and.more1) then
                  if(.not.more2) then
                    nmore2=nmore2+100
                    if(jnc-nmore2.ge.1) goto 555
                  endif
                else if(nlen.eq.(nmore+1).and.more2) then
                  if(.not.more1) then
                    nmore2=nmore2+100
                    if(jnc+nmore-nmore2+nhalf-1.le.nn) goto 555
                  endif
                endif
                a(nvar)=a(nvar)+slutn(1)
                do jnc=1,nn
                   residu(jnc)=residu(jnc)-slutn(1)
                end do
              endif
            endif
          endif

          do jnc=1,nn
             residu(jnc)=abs(residu(jnc))
          end do
          dist2=rffindq(residu,nn,nhalf,index2)
c     9555
          do step=1,kstep
            tottimes=tottimes+1
            do j=1,nhalf
              temp(j)=index2(j)
            end do
            call rfishsort(temp,nhalf)
            do j=1,nhalf
              if(.not.part.or.final) then
                do mm=1,nvad
                  datt(j,mm)=dat(temp(j),mm)
                end do
              else
                do mm=1,nvad
                  datt(j,mm)=dath(temp(j),mm)
                end do
              endif
            end do
            call rflsreg(nvmax1, nvmax,nvar,n,a,datt, weights, da, h,
     *           fckw,hvec,nvm11,jmiss,nvad,nn)

            do jnc=1,nn
              residu(jnc)=0.D0
              do j=1,nvar
                if(part.and..not.final) then
                  residu(jnc)=residu(jnc)+a(j)*dath(jnc,j)
                else
                  residu(jnc)=residu(jnc)+a(j)*dat(jnc,j)
                endif
              end do
              if(part.and..not.final) then
                residu(jnc)=dath(jnc,nvad)-residu(jnc)
              else
                residu(jnc)=dat(jnc,nvad)-residu(jnc)
              endif
              aw(jnc)=residu(jnc)
            end do

            more1=.false.
            more2=.false.
            nmore=200
            nmore2=nmore/2

            if(intadjust.eq.1) then

CDDD  CALL INTPR('>>> INTERCEPT ADJUSTMENT 2',-1,step,1)
              if(intercept .eq. 1 .and. ((.not.fine.and.part) .or.
     *             .not.part.or.((nn-nhalf).le.nmore))) then
                call rfshsort(aw,nn)
                call rfmcduni(aw,nn,nhalf,slutn,bstd,am,am2,
     *               factor,nn-nhalf+1)
                a(nvar)=a(nvar)+slutn(1)
                do jnc=1,nn
                  residu(jnc)=residu(jnc)-slutn(1)
                end do
              else if(intercept.eq.1) then
                call rfshsort(aw,nn)
                do jj=1,nn
                  am2(jj)=abs(aw(jj))
                end do
                dist2=rffindq(am2,nn,nhalf,index1)
                do jj=1,nhalf
                  aw2(jj)=aw(index1(jj))
                end do
                dist2=rffindq(aw2,nhalf,1,index2)
                jnc=index1(index2(1))
                if(jnc+nmore-nmore2+nhalf-1.gt.nn.or.jnc-nmore2.lt.1)
     *               then
                  call rfmcduni(aw,nn,nhalf,slutn,bstd,am,am2,
     *                 factor,nn-nhalf+1)
                  a(nvar)=a(nvar)+slutn(1)
                  do jnc=1,nn
                    residu(jnc)=residu(jnc)-slutn(1)
                  end do
                else
c---              repeat { ....
 666              continue
                  do jj=0,nhalf-1+nmore
                    aw2(jj+1)=aw(jnc-nmore2+jj)
                  end do
                  nlen=nmore+1
                  call rfmcduni(aw2,nhalf+nmore,nhalf,slutn,bstd,
     *                 am,am2,factor,nlen)
                  if(nlen.eq.1.and..not.more1) then
                    if(.not.more2) then
                      nmore=nmore2
                      nmore2=nmore2+nmore2
                      more1=.true.
                      if(jnc-nmore2.ge.1) goto 666
                    endif
                  else if(nlen.eq.(nmore+1).and..not.more2) then
                    if(.not.more1) then
                      nmore=nmore2
                      nmore2=-nmore2
                      more2=.true.
                      if(jnc+nmore-nmore2+nhalf-1.le.nn) goto 666
                    endif
                  else if(nlen.eq.1.and.more1) then
                    if(.not.more2) then
                      nmore2=nmore2+100
                      if(jnc-nmore2.ge.1) goto 666
                    endif
                  else if(nlen.eq.(nmore+1).and.more2) then
                    if(.not.more1) then
                      nmore2=nmore2+100
                      if(jnc+nmore-nmore2+nhalf-1.le.nn) goto 666
                    endif
                  endif
c                 } (end repeat)

c                 update a[nvar]
                  a(nvar)=a(nvar)+slutn(1)
                  do jnc=1,nn
                    residu(jnc)=residu(jnc)-slutn(1)
                  end do
                endif
              endif
            endif

            do jnc=1,nn
              residu(jnc)=abs(residu(jnc))
            end do
            dist2=rffindq(residu,nn,nhalf,index2)
            fckw=0.D0
            do jnc=1,nhalf
              fckw=fckw+residu(jnc)**2
            end do

            if(step.ge.2 .and. fckw.eq.fckw1) exit ! break {step loop}
            fckw1=fckwi
            fckwi=fckw
            if(final .and. ((i.eq.1.and.step.eq.1 .and. .not.fine)
     *           .or.fckw.lt.object)) then
              if(i_trace .ge. 3) then
                 call dblepr('Setting objfct= ', -1, fckw,1)
              endif
              object=fckw
              objfct=fckw
              do jjj=1,nhalf
                inbest(jjj)=index2(jjj)
              end do
              call rfcovcopy(a,bmeans,nvar,1)
            endif

          end do ! step in 1..kstep

c 5000
          if(.not. final) then
            if(part .and. .not. fine) then
              iii=ii
            else
              iii=1
cc		   At the end of the algorithm, only the ten
cc		   best solutions need to be stored.
            endif

            if( flag((iii-1)*10+1).eq.1) then
              lll=1
            else
              lll=2
            endif

            do j = lll,10
              if (fckw .le. mcdndex(j,2,iii)) then
                if(fckw.ne.mcdndex(j,2,iii)) then
                  if(.not.fine.and.part) goto 203
                  goto 205
                else
                  do k=j,10
                    if(fckw.eq.mcdndex(k,2,iii)) then
                      do jjj=1,nvar
                        if(part.and..not.fine) then
                          if(a(jjj).ne.m1stock((iii-1)*10+
     *                         k,jjj)) then
                            goto 203
                          endif
                        else
                          if(a(jjj).ne.mstock(k,jjj)) goto 205
                        endif
                      end do
                    endif
                  end do
                endif
                exit ! j-loop
c               ----

cvvv            using m1stock[,]
 203            do k=10,j+1,-1
                  do kk=1,nvar
                    m1stock((iii-1)*10+k,kk)=
     *                   m1stock((iii-1)*10+k-1,kk)
                  end do
                  mcdndex(k,1,iii)=mcdndex(k-1,1,iii)
                  mcdndex(k,2,iii)=mcdndex(k-1,2,iii)
                end do
                do kk=1,nvar
                  m1stock((iii-1)*10+j,kk)=a(kk)
                end do
                mcdndex(j,1,iii)=i
                mcdndex(j,2,iii)=fckw
                exit ! j-loop
c               ----

cvvv            using mstock[,]
 205            do k=10,j+1,-1
                  do kk=1,nvar
                    mstock(k,kk)= mstock(k-1,kk)
                  end do
                  mcdndex(k,1,iii)=mcdndex(k-1,1,iii)
                  mcdndex(k,2,iii)=mcdndex(k-1,2,iii)
                end do
                do kk=1,nvar
                  mstock(j,kk)=a(kk)
                end do
                mcdndex(j,1,iii)=i
                mcdndex(j,2,iii)=fckw
                exit ! j-loop
c               ----
              endif
            end do ! j-loop
          endif

 1000   continue
c.....           end for( i = 1..nrep )


 1111 continue
c---- -------- end for( ii = 1..ngroup )

cc
      if(part .and. .not. fine) then
        fine=.true.
        goto 5555
      endif
      if(.not. final .and. (.not.part .or. fine)) then
        final=.true.
        goto 5555
      endif

 9999 continue
      call rndend
C     ------ == PutRNGstate() in C
      return
      end
ccccc end {rfltsreg} ===================================================

      subroutine rfstatis(x,xmed,xmad,aw2,intercept,nvad,nvmax1,
     *     nmax,n,nstop,MADeps,weights,y,nvar,index2)
cc
      implicit none
      integer intercept, nvad,nvmax1, nmax, n, nstop, nvar
      double precision xmed(nvmax1), x(n,nvad), xmad(nvmax1)
      double precision aw2(nmax)
      double precision weights(nmax)
      double precision y(nmax)
      double precision MADeps
      double precision rfamdan
      integer index2(nmax)
c
      integer j,jnc
cc
      nstop=0
c     nstop=0: success;  =1 : "problem": mad ~= 0

      if (intercept.eq.0) then
c       regression without intercept
        do 50 j=1,nvad
          xmed(j)=0.0
          do jnc=1,n
            aw2(jnc)=abs(x(jnc,j))
          end do

          xmad(j)=rfamdan(nmax,aw2,n,index2)*1.4826
          if(abs(xmad(j)) .le. MADeps) then
            xmad(j)=0.0
            do jnc=1,n
              xmad(j)=xmad(j)+aw2(jnc)
            end do
            xmad(j)=(xmad(j)/n)*1.2533
            if(abs(xmad(j)) .le. MADeps) then
              nstop=1
              return
            endif
          endif
          do jnc=1,n
            x(jnc,j)=x(jnc,j)/xmad(j)
          end do
 50     continue

      else
c     	regression with intercept
        xmed(nvar)=0.D0
        xmad(nvar)=1.D0
        do 120 j=1,nvad
          if(j.eq.nvar) goto 120
          do jnc=1,n
            aw2(jnc)=x(jnc,j)
          end do
          xmed(j)=rfamdan(nmax,aw2,n,index2)
          do jnc=1,n
            aw2(jnc)=abs(aw2(jnc)-xmed(j))
          end do
          xmad(j)=rfamdan(nmax,aw2,n,index2)*1.4826
          if(abs(xmad(j)) .le. MADeps) then
            xmad(j)=0.0
            do jnc=1,n
              xmad(j)=xmad(j)+aw2(jnc)
            end do
            xmad(j)=(xmad(j)/n)*1.2533
            if(dabs(xmad(j)) .le. MADeps) then
              nstop=1
              return
            endif
          endif

          do jnc=1,n
            x(jnc,j)=(x(jnc,j)-xmed(j))/xmad(j)
          end do
 120    continue

      endif

      do jnc=1,n
        weights(jnc)=1.0
        y(jnc)=x(jnc,nvad)
      end do
      return
      end
cc
      function rfamdan(nmax,aa,n,index2)
cc
      implicit none
      integer nmax, n, index2(nmax)
      double precision aa(n)

      double precision rffindq, rfamdan
      integer jndl
cc
      jndl=int(n/2.0)
      if(mod(n,2).eq.0) then
        rfamdan=(rffindq(aa,n,jndl,  index2)+
     *           rffindq(aa,n,jndl+1,index2))/2.0
      else
        rfamdan=rffindq(aa,n,jndl+1,index2)
      endif
      return
      end
cc
      subroutine rflsreg(nvmax1, nvmax,k, n, f, x, w, da, h,fckw,
     *  hvec,nvm11,jmiss,nvad,nnn)
cc Arguments
      implicit none
      integer nvmax1, nvmax, k, n, nvm11, nvad, nnn
      integer jmiss(nvmax1)
      double precision x(n,nvad), f(k), w(n), da(k)
      double precision hvec(nvm11), h(nvmax,nvmax1), fckw
cc External Functions
      double precision rfqlsrg
cc Var
      double precision dfckw,dfact, dwjnc,dyj,dfka
      double precision ank,anul, hda
      integer j,l, jnc,ka,kplus, mm
cc
      kplus=k+1
      do jnc=1,k
        do j=1,kplus
          h(jnc,j)=0.D0
        end do
      end do
      anul=0.0
      do 30 jnc=1,nnn
        call rffcn(k,f,x,jnc,n,nvad)
        dwjnc=dble(w(jnc))
        anul=anul+w(jnc)
        dyj=dble(x(jnc,kplus))
        do ka=1,k
          dfka=dble(f(ka))
          h(ka,k+1)=h(ka,k+1)+dwjnc*dfka*dyj
          do l=1,ka
            h(ka,l)=h(ka,l)+dwjnc*dfka*dble(f(l))
          end do
        end do
 30   continue
      do j=1,k
        do jnc=1,j
          h(jnc,j)=h(j,jnc)
        end do
      end do
      call rfmatnv(h,nvmax,nvmax1,hvec,nvm11,k,1,jmiss)
      mm=k+1
      fckw = rfqlsrg(k,n,nvmax1,nvmax,f,x, w,h,mm,nvad,nnn)
      do jnc=1,k
        f(jnc)=h(jnc,k+1)
      end do
      dfckw=dble(fckw)
      ank=anul-k
      dfact=dble(ank)
      dfact=dfckw/dfact
      do jnc=1,k
        do j=1,k
          h(jnc,j)=h(jnc,j)*dfact
        end do
      end do
      do jnc=1,k
        hda=h(jnc,jnc)
        da(jnc)=sqrt(hda)
      end do
      return
      end
ccccc
ccccc
      subroutine rffcn(k,f,x,jnc,n,nvad)
cc
      implicit none
      integer k, jnc,n,nvad, j
      double precision f(k), x(n,nvad)
cc
      do 10,j=1,k
        f(j)=x(jnc,j)
 10   continue
      return
      end
ccccc
ccccc
      subroutine rfmatnv(an,nvmax,nvmax1,hvec,nvm11,na,nb, jmiss)
cc
      implicit none
      integer nvmax,nvmax1,nvm11, na,nb
      integer jmiss(nvmax1)
      double precision an(nvmax,nvmax1), hvec(nvm11)

      double precision deter,turn,swap
      integer j,n, nc, jcl,jdelc,jdla,jdlb,jdm,jhfd,
     *     jnc,jncb,jncc,jncd,jnk,jpaal, nma,npnb, ldel

      deter=1.0D0
      n=na
      npnb=n+nb
      jnk=0
      do j=1,npnb
        jnk=(j-1)*nvmax
        do nc=1,nvmax
          jnk=jnk+1
          hvec(jnk)=an(nc,j)
        end do
      end do
      ldel=0
      jdm=nvmax
      nma=n-1
      jdelc=1-jdm
      do 130 jhfd=1,n
        turn=0.0D0
        jdelc=jdelc+jdm
        jdla=jdelc+jhfd-1
        jdlb=jdelc+nma
        do 40 jncb=jdla,jdlb
          if(dabs(hvec(jncb)) .gt. dabs(turn)) then
            turn=hvec(jncb)
            ldel=jncb
          endif
 40     continue
        if (turn .eq. 0) goto 180

        jpaal=ldel-jdelc+1
        jmiss(jhfd)=jpaal
        if(jpaal .gt. jhfd) then
          deter=-deter
          jpaal=jpaal-jdm
          jncd=jhfd-jdm
          do 70 jnc=1,npnb
            jpaal=jpaal+jdm
            jncd=jncd+jdm
            swap=hvec(jncd)
            hvec(jncd)=hvec(jpaal)
            hvec(jpaal)=swap
 70       continue
        endif
        deter=deter*turn
        turn=1.0D0/turn
        jncd=jdelc+nma
        do jnc=jdelc,jncd
          hvec(jnc)=-hvec(jnc)*turn
        end do
        hvec(jdla)=turn
        jncb=jhfd-jdm
        jpaal=1-jdm
        do 120 jnc=1,npnb
          jpaal=jpaal+jdm
          jncb=jncb+jdm
          if(jnc .ne. jhfd) then
            jcl=jpaal+nma
            swap=hvec(jncb)
            jncd=jdelc-1
            do jncc=jpaal,jcl
              jncd=jncd+1
              hvec(jncc)=hvec(jncc)+swap*hvec(jncd)
            end do
            hvec(jncb)=swap*turn
          endif
 120    continue
 130  continue
      do 160 jncb=1,n
        jhfd=n+1-jncb
        ldel=jmiss(jhfd)
        if(ldel .ne. jhfd) then
          jpaal=(ldel-1)*jdm+1
          jcl=jpaal+nma
          jdelc=(jhfd-1)*jdm+1-jpaal
          do jncc=jpaal,jcl
            jncd=jncc+jdelc
            swap=hvec(jncc)
            hvec(jncc)=hvec(jncd)
            hvec(jncd)=swap
          end do
        endif
 160  continue
c---
 180  jnk=0
      do j=1,npnb
        do nc=1,nvmax
          jnk=jnk+1
          an(nc,j)=hvec(jnk)
        end do
      end do
      return
      end
ccccc
ccccc
      double precision
     * function rfqlsrg(k,n,nvmax1,nvmax,f,x, w,h,mm,nvad,nnn)
cc
      implicit none
      integer k,n,nvmax1,nvmax, mm,nvad,nnn
      double precision f(k), x(n,nvad), w(n), h(nvmax,nvmax1)

      double precision q,hsum
      integer jnc,jncb
cc
      q=0.D0
      do 30 jnc=1,nnn
        call rffcn(k,f,x,jnc,n,nvad)
        hsum=0.D0
        do jncb=1,k
           hsum=h(jncb,mm)*f(jncb)+hsum
        enddo
        q=(hsum-x(jnc,mm))*(hsum-x(jnc,mm))*w(jnc)+q
 30   continue
      rfqlsrg=q
      return
      end
ccccc
ccccc
      subroutine rfrtran(nvar,jcst,nfac,nvad,nvmax1,
     *  xmed,xmad, aa,jal,fckw)
cc
      implicit none
      integer nvar,jcst,nfac,nvad,nvmax1, jal
      double precision aa(jal), xmed(nvmax1), xmad(nvmax1), fckw
c  Var
      integer j

      if(nvar.le.1) then
        aa(1)=aa(1)*xmad(nvad)/xmad(1)
      else
         do j=1,nfac
            aa(j)=aa(j)*xmad(nvad)/xmad(j)
         end do
         if(jcst.eq.0) then
            aa(nvar)=aa(nvar)*xmad(nvad)/xmad(nvar)
         else
            aa(nvar)=aa(nvar)*xmad(nvad)
            do j=1,nfac
               aa(nvar)=aa(nvar)-aa(j)*xmed(j)
            end do
            aa(nvar)=aa(nvar)+xmed(nvad)
         endif
      endif
      fckw=fckw*(xmad(nvad)*xmad(nvad))
      return
      end
ccccc
ccccc
      subroutine rftrc(h,da,nvmax,nvmax1,nvar,jcst,nfac,nvad,
     *  xmed,xmad)
cc
      implicit none
      integer nvmax,nvmax1,nvar,jcst,nfac,nvad
      double precision h(nvmax,nvmax1), da(nvmax)
      double precision xmed(nvmax1),xmad(nvmax1)

      double precision xmp2,hnn
      integer j,k, k2
cc
      xmp2=dble(xmad(nvad))*dble(xmad(nvad))
      if(jcst.eq.0) then
        do 10 j=1,nvar
           do k=1,j
              h(j,k)=h(j,k)*(xmp2/(dble(xmad(j))*dble(xmad(k))))
           end do
           da(j)=dsqrt(h(j,j))
 10     continue
      else
        do j=1,nvar
           h(j,nvad)=h(j,j)
        end do
        do 30, j=1,nvar
          do k=1,j
            h(j,k)=h(j,k)*xmp2/(dble(xmad(j))*dble(xmad(k)))
          end do
          da(j)=dsqrt(h(j,j))
 30     continue
        do 50 k=1,nfac
          h(nvar,k)=h(k,nvar)*xmp2/dble(xmad(k))
          do 60 k2=1,nvar
            if(k.eq.k2) then
               h(nvar,k)=h(nvar,k)-dble(xmed(k))*xmp2/
     *              (dble(xmad(k2))*dble(xmad(k)))*h(k2,nvad)
            else if(k.lt.k2) then
               h(nvar,k)=h(nvar,k)-(dble(xmed(k2))*xmp2)/
     *              (dble(xmad(k2))*dble(xmad(k)))*h(k,k2)
            else ! k > k2
               h(nvar,k)=h(nvar,k)-dble(xmed(k2))*xmp2/
     *              (dble(xmad(k2))*dble(xmad(k)))*h(k2,k)
            endif
 60       continue
 50     continue
        h(nvar,nvar)=h(nvar,nvad)*xmp2
        do 70 k=1,nvar
          h(nvar,nvar)=h(nvar,nvar)+
     *      (dble(xmed(k))*dble(xmed(k)))*xmp2/
     *      (dble(xmad(k))*dble(xmad(k)))*h(k,nvad)
 70     continue
        do 80 k=1,nvar
           if(k.ne.nvar) then
             h(nvar,nvar)=h(nvar,nvar)-2.0D0*xmp2*dble(xmed(k))/
     *         (dble(xmad(k)))*h(k,nvar)
           else
             h(nvar,nvar)=h(nvar,nvar)-2.0D0*xmp2*dble(xmed(k))/
     *         (dble(xmad(k)))*h(nvar,nvad)
           endif
 80     continue
        do j=1,nfac
          do k=j+1,nvar
            hnn=2.0D0*dble(xmed(j))*dble(xmed(k))*xmp2
            h(nvar,nvar)=h(nvar,nvar)+hnn/
     *           (dble(xmad(j))*dble(xmad(k)))*h(j,k)
          end do
        end do
        da(nvar)=dsqrt(h(nvar,nvar))
      endif
      return
      end
ccccc
ccccc
      subroutine rfequat(am,nvmax,nvmax1, hvec,nvm11,na,nb,nerr)

      implicit none
      integer nvmax,nvmax1, nvm11, na,nb,nerr
      double precision am(nvmax,nvmax1), hvec(nvm11)

      double precision turn,swap,deter
      integer j,n, ldel, jbegc,jbegx,jdel,jdm,jendc,jendx,jhfd,jmat,
     *     jnc,jncb,jncc,jncd,jnce,jncf,jnk,jrow, lclpl, nc,neqa,nznde

      ldel=0
      jdm=nvmax
      deter=1.0D0
      n=na
      jmat=n+nb
      jnk=0
      do j=1,jmat
        jnk=(j-1)*nvmax
        do nc=1,nvmax
          jnk=jnk+1
          hvec(jnk)=am(nc,j)
        end do
      end do
      nznde=n-1
      lclpl=-jdm

      do 120 jhfd=1,n
         turn=0.D0
         lclpl=lclpl+jdm+1
         jdel=lclpl+n-jhfd
         do jncb=lclpl,jdel
           if(dabs(hvec(jncb)) .gt. dabs(turn)) then
             turn=hvec(jncb)
             ldel=jncb
           endif
         end do
         if(dabs(turn) .le. 1D-8) then
            nerr=-1
            goto 180
         endif

         if(ldel .ne. lclpl) then
            deter=-deter
            ldel=ldel-jdm
            jncb=lclpl-jdm
            do jncc=jhfd,jmat
              ldel=ldel+jdm
              jncb=jncb+jdm
              swap=hvec(jncb)
              hvec(jncb)=hvec(ldel)
              hvec(ldel)=swap
            end do
         end if
         deter=deter*turn

         if(jhfd.eq.n) goto 120
         turn=1./turn
         jncb=lclpl+1
         do jncc=jncb,jdel
           hvec(jncc)=hvec(jncc)*turn
         end do
         jncd=lclpl
         jrow=jhfd+1
         do jncb=jrow,n
            jncd=jncd+1
            jnce=lclpl
            jncf=jncd
            do jncc=jrow,jmat
              jnce=jnce+jdm
              jncf=jncf+jdm
              hvec(jncf)=hvec(jncf)-hvec(jnce)*hvec(jncd)
            end do
         end do
 120  continue

      nerr=0
      neqa=n+1
      jbegx=nznde*jdm+1
      do 150 jnc=neqa,jmat
         jbegx=jbegx+jdm
         jendx=jbegx+n
         jbegc=n*jdm+1
         jendc=jbegc+nznde
         do 140 jncb=1,nznde
            jendx=jendx-1
            jbegc=jbegc-jdm
            jendc=jendc-jdm-1
            hvec(jendx)=hvec(jendx)/hvec(jendc+1)
            swap=hvec(jendx)
            jncd=jbegx-1
            do jncc=jbegc,jendc
               jncd=jncd+1
               hvec(jncd)=hvec(jncd)-hvec(jncc)*swap
            end do
 140     continue
         hvec(jbegx)=hvec(jbegx)/hvec(1)
 150  continue
      jnc=-jdm
      jbegx=nznde*jdm+1
      jendx=jbegx+nznde
      do 160 jncb=neqa,jmat
         jbegx=jbegx+jdm
         jendx=jendx+jdm
         jnc=jnc+jdm
         jncd=jnc
         do jncc=jbegx,jendx
            jncd=jncd+1
            hvec(jncd)=hvec(jncc)
         end do
 160  continue

 180  jnk=0
      do j=1,jmat
         do nc=1,nvmax
            jnk=jnk+1
            am(nc,j)=hvec(jnk)
         end do
      end do
      return
      end
ccccc
C--     VT-- The following functions were added
C--
C--  MM: moved to ./rf-common.f - since they are used from ./rffastmcd.f too

