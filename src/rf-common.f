c
c--   Routines common to
c--   fastLTS ( ./rfltsreg.f )  and
c--   fastMCD ( ./rffastmcd.f )
c
c
      subroutine rfrangen(n, nsel, index)
c
c     Randomly draws nsel cases out of n cases.
c     Here, index is the index set.
c
      implicit none
      integer n, nsel, index(nsel)
c     unifrnd() == R C API's  unif_rand()  --> see ./R-rng4ftn.c
      double precision unifrnd
      integer i,j, num
c
      do i=1,nsel
cOLD 10  num=int(uniran(seed)*n)+1
 10      num=int(unifrnd()*n)+1
C        if(num .gt. n) then
C           call intpr('** rfrangen(): num > n; num=', -1, num, 1)
C           num=n
C        endif
         if(i.gt.1) then
            do j=1,i-1
               if(index(j).eq.num) goto 10
            end do
         endif
         index(i)=num
      end do
      return
      end
c     ---------------------------------------------------------

cOLD    function uniran(seed)
cOLD cc
cOLD cc  Draws a random number from the uniform distribution on [0,1].
cOLD cc
cOLD    real uniran
cOLD    integer seed
cOLD    integer quot
cOLD cc
cOLD    seed=seed*5761+999
cOLD    quot=seed/65536
cOLD    seed=seed-quot*65536
cOLD    uniran=float(seed)/65536.D0
cOLD    return
cOLD    end
c     ---------------------------------------------------------

      subroutine rfgenpn(n,nsel,index)
cc
cc    Constructs all subsets of nsel cases out of n cases.
cc
      implicit none
      integer n,nsel,index(nsel)
cc
      integer k,i

      k=nsel
      index(k)=index(k)+1
c while
 10   if(k.eq.1 .or. index(k).le.(n-(nsel-k))) goto 100
      k=k-1
      index(k)=index(k)+1
      do i=k+1,nsel
         index(i)=index(i-1)+1
      end do
      goto 10
c end{while}
 100  return
      end
c     ---------------------------------------------------------

      subroutine rfshsort(a,n)
cc
cc  Sorts the array a of length n.
cc
      implicit none
      integer n
      double precision a(n)
c
      double precision t
      integer gap, i,j, nextj

      gap=n
c --- repeat
 100  gap=gap/2
      if(gap.eq.0) goto 200
      do 180 i=1,n-gap
         j=i
 120     if(j.lt.1) goto 180
         nextj=j+gap
         if(a(j).gt.a(nextj)) then
            t=a(j)
            a(j)=a(nextj)
            a(nextj)=t
         else
            j=0
         endif
         j=j-gap
         goto 120
 180  continue
      goto 100
c     ---- --- end repeat
 200  return
      end
c     ---------------------------------------------------------

      subroutine rfishsort(a,n)
cc
cc  Sorts the integer array a of length n.
cc
      implicit none
      integer n, a(n)
c
      integer t, gap, i,j, nextj

      gap=n
c --- repeat
 100  gap=gap/2
      if(gap.eq.0) goto 200
      do 180 i=1,n-gap
         j=i
 120     if(j.lt.1) goto 180
         nextj=j+gap
         if(a(j).gt.a(nextj)) then
            t=a(j)
            a(j)=a(nextj)
            a(nextj)=t
         else
            j=0
         endif
         j=j-gap
         goto 120
 180  continue
      goto 100
c     ---- --- end repeat
 200  return
      end
c     ---------------------------------------------------------

      integer function replow(k)
cc
cc    Find out which combinations of n and p are
cc    small enough in order to perform exaustive search
cc    Returns the maximal n for a given p, for which
cc    exhaustive search is to be done
cc
cc    k is the number of variables (p)
cc
      implicit none
      integer k
c
      integer irep(6)
      data irep/500,50,22,17,15,14/
c
      if(k .le. 6) then
         replow = irep(k)
      else
         replow = 0
      endif
      return
      end
c     ---------------------------------------------------------

      integer function rfncomb(k,n)
cc
cc  Computes the number of combinations of k out of n.
cc  (To avoid integer overflow during the computation,
cc  ratios of reals are multiplied sequentially.)
cc  For comb > 1E+009 the resulting 'comb' may be too large
cc  to be put in the integer 'rfncomb', but the main program
cc  only calls this function for small enough n and k.
cc
      implicit none
      integer k,n
c
      double precision comb,fact
      integer j
c
      comb=dble(1.0)
      do j=1,k
         fact=(dble(n-j+1.0))/(dble(k-j+1.0))
         comb=comb*fact
      end do
c     Should give error now instead of integer overflow!
c     Don't know how to get .Machine$integer.max in Fortran, portably
      if(comb .gt. 2147483647) then
         comb=2147483647.
         call
     +    dblepr('** too many combinations; using max.integer instead:',
     +           -1,comb,1)
      endif
      rfncomb=int(comb+0.5D0)
      return
      end
c     ---------------------------------------------------------

      subroutine rfcovcopy(a,b,n1,n2)
cc
cc  Copies matrix a to matrix b.
cc
      double precision a(n1,n2)
      double precision b(n1,n2)
c
      do i=1,n1
         do j=1,n2
            b(i,j)=a(i,j)
         end do
      end do
      return
      end
c     ---------------------------------------------------------

      double precision function rffindq(aw, ncas, k, index)
c
c     Finds the k-th order statistic of the array aw[1..ncas],
c     sorting the array aw[.] until aw[k] is sure to contain the k-th value
c
c     MM{FIXME}: "rather" use R's C API   rPsort (double* X, int N, int K)

      implicit none
      integer ncas,k,index(ncas)
      double precision aw(ncas)
c
      double precision ax,wa
      integer i,j,l,lr,jnc
c
      do j=1,ncas
        index(j)=j
      end do
c     lower (= l) and upper ( =lr ) bounds:
      l=1
      lr=ncas

c--- while(l < lr)
 20   if(l .lt. lr) then
         ax=aw(k)
         jnc=l
         j=lr

c---   while(jnc < j)
 30      if(jnc .le. j) then
 40         if(aw(jnc).ge.ax) goto 50
            jnc=jnc+1
            goto 40
 50         if(aw(j).le.ax) goto 60
            j=j-1
            goto 50
 60         if(jnc .le. j) then ! swap jnc <--> j
               i=index(jnc)
               index(jnc)=index(j)
               index(j)=i
               wa=aw(jnc)
               aw(jnc)=aw(j)
               aw(j)=wa
               jnc=jnc+1
               j=j-1
            endif
            goto 30
         end if

         if(j.lt.k) l=jnc
         if(k.lt.jnc) lr=j
         goto 20
      end if

      rffindq=aw(k)
      return
      end
c     ---------------------------------------------------------

      subroutine rfrdraw(a,n,ntot,mini,ngroup,kmini)
cc
cc  Draws ngroup nonoverlapping subdatasets out of a dataset of size n,
cc  such that the selected case numbers are uniformly distributed from 1 to n.
cc
      implicit none
      integer n, ntot, kmini, a(2,ntot), mini(kmini), ngroup
c     unifrnd() == R C API's  unif_rand()  --> see ./R-rng4ftn.c
      double precision unifrnd
c
      integer jndex, nrand, k,m,i,j
cc
      jndex=0
      do k=1,ngroup
         do 20 m=1,mini(k)
cOLD        nrand=int(uniran(seed)*(n-jndex))+1
            nrand=int(unifrnd()*(n-jndex))+1
C           if(nrand .gt. n-jndex) then
C              call intpr(
C      1         '** rfrdraw(): need to correct nrand > n-jndex; nrand=',
C      2                          -1, nrand, 1)
C              nrand=n-jndex
C           endif

            jndex=jndex+1
            if(jndex.eq.1) then
               a(1,jndex)=nrand
               a(2,jndex)=k
            else
               a(1,jndex)=nrand+jndex-1
               a(2,jndex)=k
               do i=1,jndex-1
                  if(a(1,i).gt.nrand+i-1) then
                     do j=jndex,i+1,-1
                        a(1,j)=a(1,j-1)
                        a(2,j)=a(2,j-1)
                     end do
                     a(1,i)=nrand+i-1
                     a(2,i)=k
                     goto 20
c                    ------- break
                  endif
               end do
            endif
 20      continue
      end do
      return
      end
c     ---------------------------------------------------------

      logical function rfodd(n)

      rfodd=.true.
      if(2*(n/2).eq.n) rfodd=.false.
      return
      end
c     ---------------------------------------------------------

c unused        function rfnbreak(nhalf,n,nvar)
c unused cc
c unused cc  Computes the breakdown value - in percent! - of the MCD estimator
c unused cc
c unused         implicit none
c unused        integer rfnbreak, nhalf, n, nvar
c unused
c unused        if (nhalf.le.(n+nvar+1)/2) then
c unused          rfnbreak=(nhalf-nvar)*100/n
c unused        else
c unused          rfnbreak=(n-nhalf+1)*100/n
c unused        endif
c unused        return
c unused        end
c     ---------------------------------------------------------

      subroutine rfmcduni(w,ncas,jqu,slutn,bstd,aw,aw2,factor,len)
cc
cc  rfmcduni : calculates the MCD in the univariate case.
cc           w contains the ordered observations
cc
c This version returns the index (jint) in 'len'
c which is used in rfltreg.f

      implicit double precision (a-h,o-z), integer(i-n)
      integer ncas, jqu, len
      double precision w(ncas), aw(ncas), aw2(ncas)
      double precision slutn(len)
cc
      sq=0.D0
      sqmin=0.D0
      ndup=1
      do j=1,ncas-jqu+1
         slutn(j)=0.D0
      end do

      do jint=1,ncas-jqu+1
         aw(jint)=0.D0
         do j=1,jqu
            aw(jint)=aw(jint)+w(j+jint-1)
            if (jint.eq.1) sq=sq+w(j)*w(j)
         end do
         aw2(jint)=aw(jint)*aw(jint)/jqu
         if (jint.eq.1) then
            sq=sq-aw2(jint)
            sqmin=sq
            slutn(ndup)=aw(jint)
            len=jint
         else
            sq=sq - w(jint-1)*w(jint-1) + w(jint+jqu-1)*w(jint+jqu-1)
     *           - aw2(jint) + aw2(jint-1)
            if(sq.lt.sqmin) then
               ndup=1
               sqmin=sq
               slutn(ndup)=aw(jint)
               len=jint
            else
               if(sq.eq.sqmin) then
                  ndup=ndup+1
                  slutn(ndup)=aw(jint)
               endif
            endif
         endif
      end do
      slutn(1)=slutn(int((ndup+1)/2))/jqu
      bstd=factor*sqrt(sqmin/jqu)
      return
      end
c     ---------------------------------------------------------
