cc -*- mode: fortran; kept-new-versions: 25; kept-old-versions: 20 -*-
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
cc
cc   Computes the MCD estimator of multivariate location and scatter.
cc   This estimator is given by the subset of h observations for which
cc   the determinant of their covariance matrix is minimal. The MCD
cc   location estimate is then the mean of those h points, and the MCD
cc   scatter estimate is their covariance matrix. This value of h may be
cc   chosen by the user; its default value is roughly n/2.
cc
cc   The MCD estimator was first introduced in:
cc
cc      Rousseeuw, P.J. (1984), "Least Median of Squares Regression,"
cc      Journal of the American Statistical Association, Vol. 79,
cc      pp. 871-881. [See page 877.]
cc
cc   The MCD is a robust estimator in the sense that the estimates are
cc   not unduly influenced by outliers in the data, even if there
cc   are many outliers. Its robustness was proved in:
cc
cc      Rousseeuw, P.J. (1985), "Multivariate Estimation with High
cc      Breakdown Point," in Mathematical Statistics and Applications,
cc      edited by  W. Grossmann, G. Pflug, I. Vincze, and W. Wertz.
cc      Dordrecht: Reidel Publishing Company, pp. 283-297.
cc
cc      Rousseeuw, P.J. and Leroy, A.M. (1987), Robust Regression and
cc      Outlier Detection, Wiley-Interscience, New York. [Chapter 7]
cc
cc   The program also computes the distance of each observation
cc   from the center (location) of the data, relative to the shape
cc   (scatter) of the data:
cc
cc   * Using the classical estimates yields the Mahalanobis distance
cc     MD(i). Often, outlying points fail to have a large Mahalanobis
cc     distance because of the masking effect.
cc
cc   * Using the MCD estimates yields a robust distance RD(i).
cc     These distances allow us to easily identify the outliers.
cc
cc   For applications of robust distances in a regression context see:
cc
cc      Rousseeuw, P.J. and van Zomeren, B.C. (1990), "Unmasking
cc      Multivariate Outliers and Leverage Points," Journal of the
cc      American Statistical Association, Vol. 85, 633-639.
cc
cc   There also a diagnostic plot is given to distinguish between
cc   regular observations, vertical outliers, good leverage points,
cc   and bad leverage points.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc   The new FAST_MCD algorithm introduced here is due to
cc
cc      Rousseeuw, P.J. and Van Driessen, K. (1997), "A Fast
cc      Algorithm for the Minimum Covariance Determinant
cc      Estimator," in preparation.
cc
cc   The algorithm works as follows:
cc
cc      The dataset contains n cases, and nvar variables are used.
cc      Let n_0 := 2 * nmini (== 600 by default).
cc      When n <  n_0, the algorithm will analyze the dataset as a whole,
cc      when n >= n_0, the algorithm will use several subdatasets.
cc
cc   1. n < n_0 : When the dataset is analyzed as a whole, a trial
cc      subsample of nvar+1 cases is taken, of which the mean and
cc      covariance matrix is calculated. The h cases with smallest
cc      relative distances are used to calculate the next mean and
cc      covariance matrix, and this cycle is repeated k1 times.
cc      [For small n we can consider all subsets of nvar+1 out of n, else
cc      the algorithm draws 500 random subsets.]
cc      Afterwards, the best 10 solutions (covariance matrices and
cc      corresponding means) are used as starting values for the final
cc      iterations. These iterations stop when two subsequent determinants
cc      become equal. (At most k3 iteration steps are taken.)
cc      The solution with smallest determinant is retained.
cc
cc   2. n > n_0 --- more than n_0 = 2*nmini cases: The algorithm
cc      does part of the calculations on (at most) kmini nonoverlapping
cc      subdatasets, of (roughly) nmini cases.
cc
cc      Stage 1: For each trial subsample in each subdataset,
cc      k1 iterations are carried out in that subdataset.
cc      For each subdataset, the 10 best solutions are stored.
cc
cc      Stage 2 considers the union of the subdatasets, called the
cc      merged set. (If n is large, the merged set is a proper subset of
cc      the entire dataset.) In this merged set, each of the 'best
cc      solutions' of stage 1 are used as starting values for k2
cc      iterations. Also here, the 10 best solutions are stored.
cc
cc      Stage 3 depends on n, the total number of cases in the
cc      dataset. If n <= 5000, all 10 preliminary solutions are iterated
cc      k3 times. If n > 5000, only the best preliminary
cc      solution is iterated, and the number of iterations decreases to 1
cc      according to n*nvar. (If n*nvar <= 100,000 we iterate k3 times,
cc      whereas for n*nvar > 1,000,000 we take only one iteration step.)
cc
cc   An important advantage of the algorithm FAST_MCD is that it allows
cc   for exact fit situations, where more than h observations lie on
cc   a hyperplane. Then the program still yields the MCD location and
cc   scatter matrix, the latter being singular (as it should be), as
cc   well as the equation of the hyperplane.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rffastmcd(dat, n,nvar, nhalff, krep, nmini,kmini,
     *     initcov, initmean,
     *     inbest, det, weight, fit, coeff, kount, adcov,
     *     temp, index1, index2, indexx, nmahad, ndist, am, am2, slutn,
     *     med, mad, sd, means, bmeans, w, fv1, fv2,
     *     rec, sscp1, cova1, corr1, cinv1, cova2, cinv2, z,
     *     cstock, mstock, c1stock, m1stock, dath,
     *     cutoff, chimed, i_trace)

cc      VT::10.10.2005 - a DATA operator was used for computing the
cc              median and the 0.975 quantile of the chisq distribution
cc              with nvar degrees of freedom. Since now we have no
cc              restriction on the number of variables, these will be
cc              passed as parameters - cutoff and chimed

      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  ALGORITHM PARAMETERS:
c
c      The number of iteration steps in stages 1,2 and 3 can be changed
c      by adapting the parameters k1, k2, and k3.

      integer k1,k2,k3, int_max
      parameter (k1=2)
      parameter (k2=2)
      parameter (k3=100)
c int_max: easily recognized, slightly smaller than 2147483647 = .Machine$integer.max
      parameter (int_max = 2146666666)
c Arguments
      integer n,nvar ! (n, p)
      integer nhalff ! == quan := h(alpha) >= n/2 "n half"
      integer krep   ! krep == nsamp
c  krep := the total number of trial subsamples
c          to be drawn when n exceeds 2*nmini;
c          krep = 0  :<==>  "exact"  <==>  all possible subsamples
c      was hardcoded krep := 500; now an *argument*
      integer kmini ! the maximal number of subdatasets   and
      integer nmini ! their minimal size

      double precision dat(n,nvar)
      double precision initcov(nvar*nvar), initmean(nvar)
      integer inbest(nhalff)
      double precision det
      integer weight(n), fit
      double precision coeff(kmini,nvar)
      integer kount
      double precision adcov(nvar*nvar)

      integer temp(n)
      integer index1(n), index2(n), indexx(n)
      double precision nmahad(n), ndist(n)
      double precision am(n), am2(n), slutn(n)

      double precision med(nvar), mad(nvar), sd(nvar), means(nvar),
     *     bmeans(nvar), w(nvar), fv1(nvar), fv2(nvar)

      double precision rec(nvar+1),
     *     sscp1((nvar+1)*(nvar+1)), corr1(nvar*nvar),
     *     cova1(nvar*nvar), cinv1(nvar*nvar),
     *     cova2(nvar*nvar), cinv2(nvar*nvar),
     *     z(nvar*nvar)

      double precision cstock(10,nvar*nvar), mstock(10,nvar),
     *     c1stock(10*kmini, nvar*nvar),
     *     m1stock(10*kmini, nvar*nvar),
     *     dath(nmini*kmini, nvar)

      double precision cutoff, chimed
      integer i_trace

c Functions from ./rf-common.f :
      integer replow
      integer rfncomb
      double precision rffindq

c ------------------------------------------------------------------

c Variables
      integer i,ii,iii, ix, j,jj,jjj, k,kk,kkk,kstep,
     *     l,lll, m,mm,minigr,
     *     nn, ngroup,nhalf,nrep,nsel, nv_2
      double precision bstd, deti,detimin1,dist,dist2, eps,
     *     object, qorder, t

c km10, nmaxi: now *variable* as nmini
      integer km10, nmaxi,
     *     ierr,matz,pnsel, tottimes, step,
     *     flag(10*kmini), mini(kmini),
     *     subdat(2, nmini*kmini)
      double precision mcdndex(10,2,kmini)
c     subndex: vector of indices;
c     length(subndex) = maximal value of n_j := mini(j) {j in 1:ngroup} below; n0 := nmini
c     mini(j) = n1 or n1+1,  where  n0 <= n1 < n_max := max_j n_j <= n1+1 <= 1+ (3 n0 - 1)/2 = (3 n0 + 1)/2
c     ==> see vignette ../vignettes/fastMcd-kmini.Rnw
      integer subndex((3*nmini + 1)/ 2)
      double precision med1,med2, percen, pivot,rfmahad,medi2
      logical all,part,fine,final,class
c     -Wall (false alarm):
      all = .true.
      part= .false.
c     Consistency correction now happens in R

      if(i_trace .ge. 2) then
         call pr1mcd(i_trace, n, nvar, nhalff, krep, nmini, kmini)
      endif

      call rndstart
C     -------- == GetRNGstate() in C

      nrep = krep
      kstep = k1
      medi2 = 0

cc  From here on, the sample size n is known.
cc  Some initializations can be made. First of all, h (= the number of
cc  observations on which the MCD is based) is given by the integer variable
cc  nhalff.
cc  If nhalff equals n, the MCD is the classical covariance matrix.
cc  The logical value class indicates this situation.
cc  The variable jbreak is the breakdown point of the MCD estimator
cc  based on nhalff observations, whereas jdefaul = (n+nvar+1)/2
cc  would be the optimal value of nhalff, with maximal breakdown point.
cc  The variable percen is the corresponding percentage (MM: rather "fraction").
cc
c     unused    jbreak=rfnbreak(nhalff,n,nvar)

      percen = dble(nhalff) / n ! the fraction, also called  'alpha'

      if(nvar.lt.5) then
         eps=1.0D-12
      else
         if(nvar.ge.5.and.nvar.le.8) then
            eps=1.0D-14
         else
            eps=1.0D-16
         endif
      endif

      class = (nhalff .ge. n)
      if(class) goto 9500 ! compute *only* the classical estimate

      if(nvar.eq.1) then
         do jj=1,n
            ndist(jj)=dat(jj,1)
         end do
         call rfshsort(ndist,n)
cc. consistency correction now happens in R code
cc.       nquant=min(int(real(((nhalff*1.D0/n)-0.5D0)*40))+1,11)
cc.       factor=faclts(nquant)
cc.       call rfmcduni(ndist,n,nhalff,slutn,bstd,am,am2, factor,
         call rfmcduni(ndist,n,nhalff,slutn,bstd,am,am2, 1.d0,
     *        n-nhalff+1)
         initmean(1)=slutn(1)
         adcov(1)=bstd
         initcov(1)=bstd
         goto 9999
      endif

cc  p >= 2   in the following
cc  ------
c These are "constants" given the arguments:
      nmaxi = nmini*kmini
      km10 = 10*kmini
      nv_2 = nvar*nvar

cc  Some initializations:
cc    matz = auxiliary variable for the subroutine rs, indicating whether
cc           or not eigenvectors are calculated
cc    nsel = number of variables + 1
cc    ngroup = number of subdatasets, is in {1,2,.., kmini}
cc    part = logical value, true if the dataset is split up
cc    fine = logical value, becomes true when the subsets are merged
cc    final = logical value, to indicate the final stage of the algorithm
cc    all = logical value, true if all (p+1)-subsets out n of should be drawn;
cc          always true for (very) small n, but also when krep=0 (special value)
cc    subdat = matrix with a first row containing indices of observations
cc             and a second row indicating the corresponding subdataset
cc
      matz=1
      nsel=nvar+1
      ngroup=1
      fine=.false.
      final=.false.
      do i=1,nmaxi
         subdat(1,i)=int_max
         subdat(2,i)=int_max
      end do

cc  Determine whether the dataset needs to be divided into subdatasets
cc  or can be treated as a whole. The subroutine rfrdraw constructs
cc  nonoverlapping subdatasets, with uniform distribution of the case numbers.
cc  For small n, the number of trial subsamples is determined.

c     part := Shall we partition the data into sub-datasets / "groups"?
      part = (krep.gt.0 .and. n .ge. (2*nmini))
      all = .not. part
      if(part) then
         do i=1,kmini
            mini(i)=0
         end do
         kstep=k1
         ngroup = n / nmini ! =: k = n % nmini (integer division)
         if(ngroup .lt. kmini) then
c          we distribute n evenly into ngroup subdatasets, of size
            mm = n / ngroup ! =: n_0 =  n % k ==> rest r = n - k*N = n-k*n_0
c          The rest r in {0,..,k-1} gives one extra obs. in the last r groups, i.e.,
c          group numbers j > jj := k - r :
            ii = n - ngroup*mm ! =: r
            jj = ngroup - ii   ! = k - r
            do j =  1,jj
               mini(j) = mm
            end do
            do j =  jj+1,ngroup
               mini(j) = mm +1
            end do
            minigr = ngroup*mm + ii
         else !  ngroup = k := floor(n/nmini) >= kmini =: k_0 :
            ngroup = kmini
            do j=1,kmini
               mini(j)=nmini
            end do
            minigr = kmini*nmini
         end if

         nhalf = int(mini(1)*percen)
         nrep = krep / ngroup ! integer division
         if(i_trace .ge. 2)
     +        call prp1mcd (n,ngroup,minigr,nhalf,nrep, mini)
         call rfrdraw(subdat,n,minigr,mini,ngroup,kmini)
      else
c  "not part" : not partitioning; either  krep == 0  or   n <= 2*nmini-1 ( = 599 by default)
         minigr=n
         nhalf=nhalff
         kstep=k1
         if(krep.eq.0 .or. n.le.replow(nsel)) then
c             use all combinations; happens iff  nsel = nvar+1 = p+1 <= 6
            nrep = rfncomb(nsel,n)
            if(i_trace .ge. 2) call intpr('*all* combinations ',-1,0,0)
         else
            nrep=krep
            all = .false.
         endif
      endif
c     seed=iseed

c     above: pr1mcd(i_trace, n, nvar, nhalff, krep, nmini, kmini)
      if(i_trace .ge. 2) then
         call pr2mcd(part, all, kstep, ngroup, minigr, nhalf, nrep)
      endif

cc
cc  Some more initializations:
cc    m1stock = matrix containing the means of the ngroup*10 best estimates
cc              obtained in the subdatasets.
cc    c1stock = matrix containing the covariance matrices of the ngroup*10
cc              best estimates obtained in the subdatasets.
cc    mstock = matrix containing the means of the ten best estimates
cc             obtained after merging the subdatasets and iterating from
cc             their best estimates.
cc    cstock = matrix containing the covariance matrices of the ten best
cc             estimates obtained after merging the subdatasets
cc             and iterating from their best estimates.
cc    means = mean vector
cc    bmeans = initial MCD location estimate
cc    sd = standard deviation vector
cc    nmahad = vector of mahalanobis distances
cc    ndist = vector of general (possibly robust) distances
cc    inbest = best solution vector
cc    index1 = index vector of subsample observations
cc    index2 = index vector of ordered mahalanobis distances
cc    indexx = temporary index vector, parallel to index1, used when
cc              generating all possible subsamples
cc    temp  = auxiliary vector
cc    flag = vector with components indicating the occurrence of a
cc           singular intermediate MCD estimate.
cc
      do j=1,nvar
         do k=1,10
            mstock(k,j)=1234567.D0
            do kk=1,kmini
               m1stock((kk-1)*10+k,j)=1234567.D0
            end do
            do i=1,nvar
               do kk=1,kmini
                  c1stock((kk-1)*10+k,(j-1)*nvar+i)=1234567.D0
               end do
               cstock(k,(j-1)*nvar+i)=1234567.D0
            end do
         end do
         means(j)=0.D0
         bmeans(j)=0.D0
         sd(j)=0.D0
      end do

      do j=1,n
         nmahad(j)=0.D0
         ndist(j)=0.D0
         index1(j)=int_max
         index2(j)=int_max
         indexx(j)=int_max
         temp(j)=int_max
      end do
      do j=1,km10
         flag(j)=1
      end do


 9500 continue
c==== ********* Compute the classical estimates **************
c
      call rfcovinit(sscp1,nvar+1,nvar+1)
      do i=1,n
         do j=1,nvar
            rec(j)=dat(i,j)
         end do
         call rfadmit(rec,nvar,sscp1)
      end do
      call rfcovar(n,nvar,sscp1,cova1,means,sd)
      do j=1,nvar
         if(sd(j).eq.0.D0)      goto 5001
      end do

      call rfcovcopy(cova1,cinv1,nvar,nvar)
      det= 0.
      do j=1,nvar
         pivot=cinv1((j-1)*nvar+j)
         det=det + log(pivot)
         if(pivot.lt.eps)       goto 5001

         call rfcovsweep(cinv1,nvar,j)
      end do
      call rfcorrel(nvar,cova1,corr1,sd)

c     if just classical estimate, we are done
      if(class)         goto 9999

      goto 5002

c     singularity '1' (exact fit := 1) :
 5001 continue
      call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
      call rfdis(dat,z,ndist,n,nvar,n,nvar,means)
      call rfexact(kount,n,ndist, nvar,
     *     sscp1,rec,dat, cova1,means,sd,weight)
      call rfcovcopy(cova1,initcov,nvar,nvar)
      call rfcovcopy(means,initmean,nvar,1)
      do j=1,nvar
         coeff(1,j)=z(j)
      end do
      fit=1
      goto 9999


 5002 continue
cc
cc Compute and store classical Mahalanobis distances.
cc
      do j=1,n
         do i=1,nvar
            rec(i)=dat(j,i)
         end do
         nmahad(j)=rfmahad(rec,nvar,means,cinv1)
      end do


cc ******* Compute the MCD estimates ************** ----------------------------

cc Main loop: inspects the subsamples.
cc   Every time the sscp of the subsample is placed in sscp1,
cc   its covariance matrix in cova1, and its inverse in cinv1 .
cc   The minimum covariance determinant matrix is placed in cova2,
cc   and its inverse in cinv2.
cc   The robust distances are placed in ndist.
cc
c    tottimes := counting the total number of iteration steps in the main loop
cc
cc   The algorithm returns here twice when the dataset is divided
cc   at the beginning of the program. According to the situation,
cc   new initializations are made.
c  fine  == TRUE : <==> We are in the second stage, where the subdatasets are merged,
c  final == TRUE : <==> We are in the last stage, when the whole dataset is considered
c     			In the last stage, the number of iterations 'nrep'
c                       is determined according to the total number of observations and the dimension.
      tottimes=0

 5555 object=10.D25

      if(.not. part .or. final) then
         nn=n
      else if (fine) then !->  part  &  fine  &  .not. final
         nn=minigr
      else !->  part - "phase 1" (.not. fine  & .not. final)
         nn=-1
      endif
      if(i_trace .ge. 2)  ! " Main loop, phase[%s]: ... "
     1     call pr3mcd(part, fine, final, nrep, nn,
     2                 nsel, nhalf, kstep, nmini, kmini)

      if(fine .or.(.not.part.and.final)) then
         nrep = 10
c        ----   == hardcoded
         nsel = nhalf
         kstep = k2
         if (final) then ! "final": stage 3 --
            nhalf=nhalff
            ngroup=1
c           ksteps := k3 (= 100) unless n*p is "large" where
c	    ksteps jumps down to at most 10 <<- "discontinuous!" FIXME
            if (n*nvar .le.100000) then
               kstep=k3  ! = 100 ("hardcoded default")
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
            else ! n*p > 1e6
               kstep=1
            endif
            if (n.gt.5000) then
               nrep=1
            endif
         else
            nhalf=int(minigr*percen)
         endif
      endif

      do i=1,nsel-1
         index1(i)=i
         indexx(i)=i
      end do
      index1(nsel)=nsel-1
      indexx(nsel)=nsel-1
cc
cc  Initialization of the matrices to store partial results. For the
cc  first stage of the algorithm, the currently best covariance matrices and
cc  means are stored in the matrices c1stock and m1stock initialized earlier.
cc  The corresponding objective values and the number of the trial subset
cc  are stored in the matrix mcdndex.
cc  For the second stage of the algorithm or for small datasets, only the
cc  currently best objective values are stored in the same matrix mcdndex
cc  and the corresponding covariance matrices and mean vectors are stored in
cc  the matrices cstock and mstock initialized earlier.
cc
      if(.not. final) then
         do i=1,10
            do j=1,ngroup
               mcdndex(i,1,j)=10.D25
               mcdndex(i,2,j)=10.D25
            end do
         end do
      endif
      if(.not.fine .and. .not.final) then !-- first phase
         do j=1,nvar
            do i=1,n
               am (i)=dat(i,j)
               am2(i)=dat(i,j)
            end do
            if(2*n/2 .eq. n) then
               med1=rffindq(am, n, n/2,   index2)
               med2=rffindq(am2,n,(n+2)/2,index2)
               med(j)=(med1+med2)/2
            else
               med(j)=rffindq(am,n,(n+1)/2,index2)
            endif
            do i=1,n
               ndist(i)=dabs(dat(i,j)-med(j))
            end do
            mad(j)=rffindq(ndist,n,nhalff,index2)
            if(mad(j)-0.D0 .lt. eps) then
               do k=1,j-1
                  do i=1,n
                     dat(i,k)=dat(i,k)*mad(k)+med(k)
                  end do
               end do
               call rfcovinit(sscp1,nvar+1,nvar+1)
               do k=1,nsel
                  do m=1,nvar
                     rec(m)=dat(index2(k),m)
                  end do
                  call rfadmit(rec,nvar,sscp1)
               end do
               call rfcovar(nsel,nvar,sscp1,cova1,means,sd)
               call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)

C       VT::15.11.2014, fixing array overrun, found by MM
C       The following code expects that z (the plane coefficients)
C       are all zeros with 1 in the position of the variable with MAD=0
C       If not, tries to find it.
C
            if(.FALSE.) then
               if(z(j).ne.1) then
                  do kk=1,nvar
                     if(z(kk*nvar+j).eq.1) then
                        do l=1,nvar
                           z(l)=z(kk*nvar+l)
                        end do
                        goto 76 ! break
                     endif
                  end do
               endif
 76            continue
            else
C       Instead of this, we set all coefficients to 0 and the one of
C       variable j to 1. The exactfit code will be set 3 and will be
C       handled respectively by the R code.
               do kk=1,nvar
                 z(kk) = 0
               end do
               z(j) = 1
            end if

               call rfdis(dat,z,ndist,n,nvar,n,nvar,means)
               call rfexact(kount,n,ndist, nvar,
     *              sscp1,rec,dat, cova1,means,sd,weight)
               call rfcovcopy(cova1,initcov,nvar,nvar)
               call rfcovcopy(means,initmean,nvar,1)
               do jjj=1,nvar
                  coeff(1,jjj)=z(jjj)
               end do
               fit=3
               goto 9999
            endif
            do i=1,n
               dat(i,j)=(dat(i,j)-med(j))/mad(j)
            end do
         end do
      endif
cc
cc  The matrix dath contains the observations to be used in the
cc  algorithm. In the first stage of the split-up procedure dath contains
cc  nmini objects, corresponding to the original observations, with the index
cc  of the processed group in the array subdat. For the second stage, the
cc  data points of all the subdatasets are merged in dath.
cc  The variable kount indicates the occurrence of a singular subsample leading
cc  to the corresponding plane. In some situations the variable kount counts
cc  the number of observations on that plane.
cc
      if (fine .and. .not. final) then
         do j=1,minigr
            do k=1,nvar
               dath(j,k)=dat(subdat(1,j),k)
            end do
         end do
      endif
      kount=0

c---- For-Loop over groups  - - - - - - - - - - - - - - - - - - - - -
      do 1111 ii= 1,ngroup
         if(.not.fine) kount=0
         if(part .and. .not. fine) then
            nn=mini(ii)
            kk=0
            do j=1,minigr
               if(subdat(2,j).eq.ii) then
                  kk=kk+1
                  subndex(kk)=subdat(1,j)
               endif
            end do
            do j=1,mini(ii)
               do k=1,nvar
                  dath(j,k)=dat(subndex(j),k)
               end do
            end do
         endif

         if(i_trace .ge. 3) call prgrmcd(ii, nn, i_trace)
         do i=1,nn
            index2(i)=i
         end do

cc  The number of trial subsamples is represented by nrep, which depends
cc  on the data situation.
cc  When all (p+1)-subsets out of n can be drawn, the subroutine rfgenpn
cc  is used. Otherwise, random subsamples are drawn by the routine
cc  rfrangen. The trial subsamples are put in the array index1. The
cc  same thing happens for large datasets, except that the number of
cc  observations is nmini instead of n.
cc
cc  When a trial subsample is singular, the algorithm counts the number of
cc  observations that lie on the hyperplane corresponding to this sample.
cc  If, for small datasets, this number is larger than nhalff, the program
cc  stops (exact fit) and gives the mean and the covariance matrix
cc  of the observations on the hyperplane, together with the equation
cc  of the hyperplane.
cc  For large datasets, the algorithm first checks whether there are more
cc  than nhalff observations on the hyperplane. If this is the case, the
cc  program stops for the same reason of exact fit and gives the covariance
cc  matrix and mean of the observations on the hyperplane. If not, the
cc  algorithm counts the number of observations that lie on the hyperplane.
cc  When this number is smaller than the current nhalf in the subdataset, these
cc  observations are extended to nhalf observations by adding those
cc  observations that have smallest orthogonal distances to the hyperplane
cc  and the algorithm continues.
cc  When larger, the coefficients of the hyperplane are stored in the matrix
cc  m1stock for use as starting value in the next stage, and the flag of this
cc  estimate gets the value zero.
cc
cc  In the second stage of the algorithm, when the subdatasets are merged,
cc  the array index2 contains the indices of the observations
cc  corresponding to the nhalf observations with minimal relative distances
cc  with respect to the best estimates of the first stage.
cc  When the estimate of the first stage is a hyperplane, the algorithm
cc  investigates whether there are more than the current nhalf observations of
cc  the merged subdataset on that hyperplane. If so, the coefficients of the
cc  hyperplane are again stored, now in the matrix mstock, for the final
cc  stage of the algorithm.
cc  If not, the observations on the hyperplane are extended to nhalf
cc  observations by adding the observations in the merged dataset with
cc  smallest orthogonal distances to that hyperplane.
cc  For small datasets or for larger datasets with n <= nmaxi := nmini*kmini,
cc  the algorithm already stops when one solution becomes singular,
cc  since we then have an exact fit.
cc
cc  In the third stage, the covariance matrices and means of the best
cc  solutions of the second stage are used as starting values.
cc  Again, when a solution becomes singular, the subroutine 'exact'
cc  determines the hyperplane through at least nhalff observations and stops
cc  because of the exact fit.
cc
cc  When the program stops because of an exact fit, the covariance matrix and
cc  mean of the observations on the hyperplane will always be given.
cc
C        VT::27.10.2014 - an issue with nsamp="exact" fixed:
         do ix=1,n
             indexx(ix)=index1(ix)
         end do

         do 1000 i=1,nrep
            pnsel=nsel
            tottimes=tottimes+1
            if(i_trace .ge. 4) call pr4mcd(i)
            call rchkusr() ! <- allow user interrupt
            deti= -1.d300
            detimin1=deti
            step=0
            call rfcovinit(sscp1,nvar+1,nvar+1)
            if((part.and..not.fine).or.(.not.part.and..not.final)) then
               if(part) then
                  call rfrangen(mini(ii),nsel,index1)
               else if(all) then
                  call rfgenpn(n,nsel,indexx)
                  do ix=1,n
                     index1(ix)=indexx(ix)
                  end do
               else
                  call rfrangen(n,nsel,index1)
               endif
            endif
cc
cc  The covariance matrix and mean of the initial subsamples are
cc  calculated with the subroutine covar and represented by
cc  the variables cova1 and means.
cc
cc  In the following stages of the algorithm, the covariance matrices and means
cc  used as starting values are already stored in the matrices c1stock
cc  and m1stock (for the second stage), and in the matrices cstock and mstock
cc  (for the third stage).
cc
cc  The inverse cinv1 of the covariance matrix is calculated by the
cc  subroutine rfcovsweep, together with its determinant det.
c
c Repeat
 9550       call rfcovinit(sscp1,nvar+1,nvar+1)
            if(.not.fine.and.part) then
               do j=1,pnsel
                  do m=1,nvar
                     rec(m)=dath(index1(j),m)
                  end do
                  call rfadmit(rec,nvar,sscp1)
               end do
               call rfcovar(pnsel,nvar,sscp1,cova1,means,sd)
            endif
            if(.not.part.and..not.final) then
               do j=1,pnsel
                  do m=1,nvar
                     rec(m)=dat(index1(j),m)
                  end do
                  call rfadmit(rec,nvar,sscp1)
               end do
               call rfcovar(pnsel,nvar,sscp1,cova1,means,sd)
            endif
            if (final) then
               if(mstock(i,1) .ne. 1234567.D0) then
                  do jj=1,nvar
                     means(jj)=mstock(i,jj)
                     do kk=1,nvar
                        cova1((jj-1)*nvar+kk)=cstock(i,(jj-1)*nvar+kk)
                     end do
                  end do
               else
                  goto 1111
               endif
               if(flag(i).eq.0) then
                  qorder=1.D0
                  do jjj=1,nvar
                     z(jjj)=coeff(1,jjj)
                  end do
                  call rfdis(dat,z,ndist,n,nvar,nn,nvar, means)
                  dist2=rffindq(ndist,nn,nhalf,index2)
                  goto 9555
               endif
            endif
            if (fine .and. .not.final) then
               if(m1stock((ii-1)*10+i,1) .ne. 1234567.D0) then
                  do jj=1,nvar
                     means(jj)=m1stock((ii-1)*10+i,jj)
                     do kk=1,nvar
                        cova1((jj-1)*nvar+kk)=c1stock((ii-1)*10+i,
     *                       (jj-1)*nvar+kk)
                     end do
                  end do
               else
                  goto 1111
               endif
               if(flag((ii-1)*10+i).eq.0) then
                  qorder=1.D0
                  do jjj=1,nvar
                     z(jjj)=coeff(ii,jjj)
                  end do
                  call rfdis(dath,z,ndist,nmaxi,nvar,nn,nvar, means)
                  call rfshsort(ndist,nn)
                  qorder=ndist(nhalf)
                  if(dabs(qorder-0.D0).lt.10.D-8 .and. kount.eq.0
     *                 .and. n.gt.nmaxi) then
                     kount=nhalf
                     do kkk=nhalf+1,nn
                        if(dabs(ndist(kkk)-0.D0).lt.10.D-8) then
                           kount=kount+1
                        endif
                     end do
                     flag(1)=0
                     do kkk=1,nvar
                        coeff(1,kkk)=z(kkk)
                     end do
                     call rfstore2(nvar,cstock,mstock,nv_2,
     *                    kmini,cova1,means,i,mcdndex,kount)
                     kount=1
                     goto 1000
                  else if(dabs(qorder-0.D0).lt.10.D-8 .and.
     *                    kount.ne.0 .and. n.gt.nmaxi) then
                     goto 1000
                  else
                     flag(1)=1
                     dist2=rffindq(ndist,nn,nhalf,index2)
                     goto 9555
                  endif
               endif
            endif
            call rfcovcopy(cova1,cinv1,nvar,nvar)
            det=0.
            do 200 j=1,nvar
               pivot=cinv1((j-1)*nvar+j)
               det=det + log(pivot)
               if(pivot.lt.eps) then
                  call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
                  qorder=1.D0
                  if(.not.part.or.final) then
                     call rfdis(dat,z,ndist,n,nvar,nn,nvar,means)
                  else
                     call rfdis(dath,z,ndist,nmaxi,nvar,nn,nvar,means)
                  endif
                  call rfshsort(ndist,nn)
                  qorder=ndist(nhalf)
                  if(dabs(qorder-0.D0).lt. 10.D-8 .and. .not.part) then
                     call transfo(cova1,means,dat,med,mad,nvar,n)
                     call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
                     call rfdis(dat,z,ndist,n,nvar,nn,nvar,means)
                     call rfexact(kount,n,ndist, nvar,
     *                    sscp1,rec,dat, cova1,means,sd,weight)
                     call rfcovcopy(cova1,initcov,nvar,nvar)
                     call rfcovcopy(means,initmean,nvar,1)
                     do jjj=1,nvar
                        coeff(1,jjj)=z(jjj)
                     end do
                     fit=2
                     goto 9999
                  else if(dabs(qorder-0.D0).lt. 10.D-8 .and. part .and.
     *                    kount.eq.0) then
                     call rfdis(dat,z,ndist,n,nvar,n,nvar, means)
                     call rfshsort(ndist,n)
                     if(dabs(ndist(nhalff)-0.D0).lt.10.D-8) then
                        call transfo(cova1,means,dat,med,mad,nvar,n)
                        call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
                        call rfdis(dat,z,ndist,n,nvar,nn,nvar,means)
                        call rfexact(kount,n,ndist, nvar,sscp1,
     *                       rec,dat, cova1,means,sd,weight)
                        call rfcovcopy(cova1,initcov,nvar,nvar)
                        call rfcovcopy(means,initmean,nvar,1)
                        do jjj=1,nvar
                           coeff(1,jjj)=z(jjj)
                        end do
                        fit=2
                        goto 9999
                     endif
                     call rfdis(dath,z,ndist,nmaxi,nvar,nn,nvar, means)
                     call rfshsort(ndist,nn)
                     kount=nhalf
                     do kkk=nhalf+1,nn
                        if(dabs(ndist(kkk)-0.D0) .lt. 10.D-8) then
                           kount=kount+1
                        endif
                     end do
                     flag((ii-1)*10+1)=0
                     do kkk=1,nvar
                        coeff(ii,kkk)=z(kkk)
                     end do
                     call rfstore1(nvar,c1stock,m1stock,nv_2,
     *                    kmini,cova1,means,i,km10,ii,mcdndex, kount)
                     kount=1
                     goto 1000
                  else if(dabs(qorder-0.D0).lt. 10.D-8 .and. part .and.
     *                    kount.ne.0) then
                     goto 1000
                  else
C
C                 VT::27.10.2014 - an issue with nsamp="exact" fixed:
C
C                 Add one more observation and return to recompute the
C                 covariance. In case of complete enumeration, when all
C                 p+1 subsamples are generated, the array 'index1' must
C                 be preserved 8around label 9550).
C
                     if(i_trace .ge. 2)
     *                call intpr('Singularity -> extended subsample: ',
     *                    -1,index1,nsel)

                     call rfishsort(index1,pnsel)
                     call prdraw(index1,pnsel, nn)
                     pnsel=pnsel+1
                     goto 9550
c                    --------- until
                  endif
               endif
               call rfcovsweep(cinv1,nvar,j)
 200        continue
cc
cc  Mahalanobis distances are computed with the subroutine rfmahad
cc  and stored in the array ndist.
cc  The k-th order statistic of the mahalanobis distances is stored
cc  in dist2. The array index2 containes the indices of the
cc  corresponding observations.
cc
            do j=1,nn
               if(.not.part.or.final) then
                  do mm=1,nvar
                     rec(mm)=dat(j,mm)
                  end do
               else
                  do mm=1,nvar
                     rec(mm)=dath(j,mm)
                  end do
               endif
               t=rfmahad(rec,nvar,means,cinv1)
               ndist(j)=t
            end do
            dist2=rffindq(ndist,nn,nhalf,index2)
cc
cc  The variable kstep represents the number of iterations of the current stage (1,2, or 3),
cc  i.e., the situation of the program, kstep = k1, k2, or k3.  Within each
cc  iteration the mean and covariance matrix of nhalf observations are
cc  calculated. The nhalf smallest corresponding mahalanobis distances
cc  determine the subset for the next iteration.
cc  The best subset for the whole data is stored in the array inbest.
cc  The iteration stops when two subsequent determinants become equal.
cc
 9555       do 400 step=1,kstep
               tottimes=tottimes+1
               if(i_trace .ge. 4) call pr5mcd(step, tottimes)
               call rchkusr() ! <- allow user interrupt
               call rfcovinit(sscp1,nvar+1,nvar+1)
               do j=1,nhalf
                  temp(j)=index2(j)
               end do
               call rfishsort(temp,nhalf)
               do j=1,nhalf
                  if(.not.part.or.final) then
                     do mm=1,nvar
                        rec(mm)=dat(temp(j),mm)
                     end do
                  else
                     do mm=1,nvar
                        rec(mm)=dath(temp(j),mm)
                     end do
                  endif
                  call rfadmit(rec,nvar,sscp1)
               end do
               call rfcovar(nhalf,nvar,sscp1,cova1,means,sd)
               call rfcovcopy(cova1,cinv1,nvar,nvar)
               det= 0.
               do 600 j=1,nvar
                  pivot=cinv1((j-1)*nvar+j)
                  det=det + log(pivot)
                  if(pivot.lt.eps) then
                     if(final .or. .not.part .or.
     *                  (fine.and. .not.final .and. n .le. nmaxi)) then
                        call transfo(cova1,means,dat,med,mad,nvar,n)
                        call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
                        if(final.or..not.part) then
                           call rfdis(dath,z,ndist,n, nvar,nn,
     *                          nvar,means)
                        else
                           call rfdis(dath,z,ndist,nmaxi,nvar,nn,
     *                          nvar,means)
                        endif
                        call rfexact(kount,n,ndist,nvar,sscp1,
     *                       rec,dat, cova1,means,sd,weight)
                        call rfcovcopy(cova1,initcov,nvar,nvar)
                        call rfcovcopy(means,initmean,nvar,1)
                        do jjj=1,nvar
                           coeff(1,jjj)=z(jjj)
                        end do
                        fit=2
                        goto 9999
                     endif
                     if(part.and..not.fine.and.kount.eq.0) then
                        call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
                        call rfdis(dat,z,ndist,n,nvar,n,nvar, means)
                        call rfshsort(ndist,n)
                        if(dabs(ndist(nhalff)-0.D0).lt.10.D-8) then
                           call transfo(cova1,means,dat,med,mad,nvar,n)
                           call rs(nvar,nvar,cova1,w,matz,z,
     *                          fv1,fv2,ierr)
                           call rfdis(dat,z,ndist,n,nvar,n,nvar,means)
                           call rfexact(kount,n,ndist,nvar,sscp1,
     *                          rec,dat, cova1,means,sd,weight)
                           call rfcovcopy(cova1,initcov,nvar,nvar)
                           call rfcovcopy(means,initmean,nvar,1)
                           do jjj=1,nvar
                              coeff(1,jjj)=z(jjj)
                           end do
                           fit=2
                           goto 9999
                        endif
                        call rfdis(dath,z,ndist,nmaxi,nvar,nn,
     *                       nvar,means)
                        call rfshsort(ndist,nn)
                        kount=nhalf
                        do,kkk=nhalf+1,nn
                           if(dabs(ndist(kkk)-0.D0).lt.10.D-8) then
                              kount=kount+1
                           endif
                        end do
                        flag((ii-1)*10+1)=0
                        do kkk=1,nvar
                           coeff(ii,kkk)=z(kkk)
                        end do
                        call rfstore1(nvar,c1stock,m1stock,nv_2,
     *                       kmini,cova1,means,i,km10,ii,mcdndex, kount)
                        kount=1
                        goto 1000
                     else
                        if(part.and..not.fine.and.kount.ne.0) then
                           goto 1000
                        endif
                     endif
                     if(fine.and..not.final.and.kount.eq.0) then
                        call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
                        call rfdis(dat,z,ndist,n,nvar,n,nvar, means)
                        call rfshsort(ndist,n)
                        if(dabs(ndist(nhalff)-0.D0).lt.10.D-8) then
                           call transfo(cova1,means,dat,med,mad,nvar,n)
                           call rs(nvar,nvar,cova1,w,matz,z,
     *                          fv1,fv2,ierr)
                           call rfdis(dat,z,ndist,n,nvar,n,nvar,means)
                           call rfexact(kount,n,ndist,nvar,sscp1,
     *                          rec,dat, cova1,means,sd,weight)
                           call rfcovcopy(cova1,initcov,nvar,nvar)
                           call rfcovcopy(means,initmean,nvar,1)
                           do jjj=1,nvar
                              coeff(1,jjj)=z(jjj)
                           end do
                           fit=2
                           goto 9999
                        endif
                        call rfdis(dath,z,ndist,nmaxi,nvar,nn,
     *                       nvar,means)

                        call rfshsort(ndist,nn)
                        kount=nhalf
                        do kkk=nhalf+1,nn
                           if(dabs(ndist(kkk)-0.D0).lt.10.D-8) then
                              kount=kount+1
                           endif
                        end do
                        flag(1)=0
                        do kkk=1,nvar
                           coeff(1,kkk)=z(kkk)
                        end do
                        call rfstore2(nvar,cstock,mstock,nv_2,
     *                       kmini,cova1,means,i,mcdndex,kount)
                        kount=1
                        goto 1000
                     else
                        if(fine.and..not.final.and.kount.ne.0) then
                           goto 1000
                        endif
                     endif
                  endif
                  call rfcovsweep(cinv1,nvar,j)
 600           continue

               if(step.ge.2 .and. det.eq.detimin1) then
                  goto 5000
               endif
               detimin1=deti
               deti=det
               do j=1,nn
                  if(.not.part.or.final) then
                     do mm=1,nvar
                        rec(mm)=dat(j,mm)
                     end do
                  else
                     do mm=1,nvar
                        rec(mm)=dath(j,mm)
                     end do
                  endif
                  t=rfmahad(rec,nvar,means,cinv1)
                  ndist(j)=t
               end do
               dist2=rffindq(ndist,nn,nhalf,index2)
               dist=dsqrt(dist2)
               if(final .and. ((i.eq.1 .and. step.eq.1 .and. .not.fine)
     *             .or. det .lt. object)) then
                  medi2=rffindq(ndist,nn,int(n/2),index1)
                  object=det
                  do jjj=1,nhalf
                     inbest(jjj)=index2(jjj)
                  end do
                  call rfcovcopy(cova1,cova2,nvar,nvar)
                  call rfcovcopy(cinv1,cinv2,nvar,nvar)
                  call rfcovcopy(means,bmeans,nvar,1)
               endif
 400        continue
            if(i_trace .ge. 4) call intpr("", -1,1,0)

cc  After each iteration, it has to be checked whether the new solution
cc  is better than some previous one and therefore needs to be stored. This
cc  isn't necessary in the third stage of the algorithm, where only the best
cc  solution is kept.

 5000       if(.not. final) then
               if(part .and. .not. fine) then
                  iii=ii
               else
                  iii=1
               endif
c     At the end of the algorithm, only the ten
c     best solutions need to be stored.

cc  For each data group :
cc    If the objective function is lower than the largest value in the
cc    matrix mcdndex :
cc    A distinction is made between different stages of the algorithm:
cc      * At the first stage of the split-up situation:
cc        -If the new objective value did not yet occur in mcdndex
cc         its value and corresponding covariance matrix and mean are
cc         stored at the right place in the matrices mcdndex, c1stock and
cc         m1stock, and other values are shifted to their new position
cc         in these arrays.
cc        -If the new objective value already occurs in mcdndex, a
cc         comparison is made between the new mean vector and covariance matrix
cc         and those estimates with the same determinant.
cc         When for an equal determinant, the mean vector or covariance matrix
cc         do not correspond, both of them are kept in the matrices mcdndex
cc         and nbest.
cc      * In the second stage of the algorithm, the covariances and means
cc        are stored :
cc        - If the new objective value did not yet occur
cc          in the matrix mcdndex, it is inserted by shifting the greater
cc          determinants upwards and doing the same in the arrays mstock
cc          and cstock.
cc        - If the new objective value already occurs in the array mcdndex,
cc          it is compared with all solutions with the same determinant.
cc          In the case of an equality, the means and covariances
cc          are compared to determine whether or not to insert the
cc          new solution.
cc    Otherwise nothing happens. When a singularity occurs,
cc    the determinant in the matrix mcdndex is zero and the
cc    corresponding flag is zero too, so the search in the arrays mcdndex,
cc    m1stock, c1stock, mstock and cstock is done on the rows with flag one.
cc

               if( flag((iii-1)*10+1).eq.1) then
                  lll=1
               else
                  lll=2
               endif
               do j=lll,10
                  if (det .le. mcdndex(j,2,iii)) then
                     if(det.ne.mcdndex(j,2,iii)) then
                        if(.not.fine.and.part)                  goto 203
                        goto 205
                     else
                        do kkk=j,10
                           if(det.eq.mcdndex(kkk,2,iii)) then
                              do jjj=1,nvar
                                 if(part.and..not.fine) then
                                    if(means(jjj) .ne.
     *                                 m1stock((iii-1)*10+ kkk,jjj))
     *                                   goto 203
                                 else
                                    if(means(jjj).ne.mstock(kkk,jjj))
     *                                   goto 205
                                 endif
                              end do
                              do jjj=1,nvar*nvar
                                 if(part.and..not.fine) then
                                    if(cova1(jjj) .ne.
     *                                   c1stock((iii-1)*10+ kkk,jjj))
     *                                   goto 203
                                 else
                                    if(cova1(jjj).ne.cstock(kkk,jjj))
     *                                   goto 205
                                 endif
                              end do
                           endif
                        end do  ! kkk
                     endif
                     goto 1000

c---
 203                 do k=10,j+1,-1
                        do kk=1,nvar*nvar
                           c1stock((iii-1)*10+k,kk)=
     *                          c1stock((iii-1)*10+k-1,kk)
                        end do
                        do kk=1,nvar
                           m1stock((iii-1)*10+k,kk)=
     *                          m1stock((iii-1)*10+k-1,kk)
                        end do
                        mcdndex(k,1,iii)=mcdndex(k-1,1,iii)
                        mcdndex(k,2,iii)=mcdndex(k-1,2,iii)
                     end do

                     do kk=1,nvar
                        do kkk=1,nvar
                           c1stock((iii-1)*10+j,(kk-1)*nvar+kkk)=
     *                          cova1((kk-1)*nvar+kkk)
                           m1stock((iii-1)*10+j,kk)=means(kk)
                        end do
                     end do
                     mcdndex(j,1,iii)=i
                     mcdndex(j,2,iii)=det
                     goto 1000
c---
 205                 do k=10,j+1,-1
                        do kk=1,nvar*nvar
                           cstock(k,kk)= cstock(k-1,kk)
                        end do

                        do kk=1,nvar
                           mstock(k,kk)= mstock(k-1,kk)
                        end do

                        mcdndex(k,1,iii)=mcdndex(k-1,1,iii)
                        mcdndex(k,2,iii)=mcdndex(k-1,2,iii)
                     end do
                     do kk=1,nvar
                        do kkk=1,nvar
                           cstock(j,(kk-1)*nvar+kkk)=
     *                          cova1((kk-1)*nvar+kkk)
                           mstock(j,kk)=means(kk)
                        end do
                     end do
                     mcdndex(j,1,iii)=i
                     mcdndex(j,2,iii)=det
                     goto 1000

                  endif
               end do  ! j

            endif
c     (not final)

 1000    continue !end{ i = 1..nrep }
 1111 continue
c---- - - - - - end [ For (ii = 1 .. ngroup) ]  - - - - - - - - -

cc  Determine whether the algorithm needs to be run again or not.
cc
      if(part .and. .not. fine) then
         fine= .true.
         goto 5555
      else if(.not. final .and. ((part.and.fine).or. .not.part)) then
         final= .true.
         goto 5555
      endif

cc******** end { Main Loop } ************** --------------------------------


c MM: 'temp' is thrown away in calling R code:
c      do j=1,nhalf
c         temp(j)=inbest(j)
c      end do
c      call rfishsort(temp,nhalf)

      do j=1,nvar
         means(j)=bmeans(j)*mad(j)+med(j)
      end do
      call rfcovcopy(means,initmean,nvar,1)

      do i=1,nvar
         do j=1,nvar
            cova1((i-1)*nvar+j)=cova2((i-1)*nvar+j)*mad(i)*mad(j)
         end do
      end do
      call rfcovcopy(cova1,initcov,nvar,nvar)
      det=object
      do j=1,nvar
         det=det + 2*log(mad(j))
      end do


cc      VT::chimed is passed now as a parameter
cc        call rfcovmult(cova1,nvar,nvar,medi2/chimed(nvar))
cc        call rfcovmult(cova2,nvar,nvar,medi2/chimed(nvar))
cc        call rfcovmult(cinv2,nvar,nvar,1.D0/(medi2/chimed(nvar)))

      medi2 = medi2/chimed
      call rfcovmult(cova1, nvar,nvar, medi2)
      call rfcovmult(cova2, nvar,nvar, medi2)
      call rfcovmult(cinv2, nvar,nvar, 1.D0/medi2)
      call rfcovcopy(cova1, adcov,nvar,nvar)
cc
cc      The MCD location is in bmeans.
cc      The MCD scatter matrix is in cova2,
cc      and its inverse in cinv2.
cc
cc      For every observation we compute its MCD distance
cc      and compare it to a cutoff value.
cc
      call rfcovinit(sscp1,nvar+1,nvar+1)

      do i=1,n
         do mm=1,nvar
            rec(mm)=dat(i,mm)
         end do
         dist2=rfmahad(rec,nvar,bmeans,cinv2)
         if(dist2.le.cutoff) then
            weight(i)=1
         else
            weight(i)=0
         endif
      end do

      call transfo(cova2,bmeans,dat,med,mad,nvar,n)
      goto 9999
cc ******************************************************************

 9999 continue
      if(i_trace .ge. 2) call pr9mcd(tottimes)
      call rndend
C     ------ == PutRNGstate() in C
      return
      end
ccccc end {rffastmcd}

ccccc
ccccc
ccccc
ccccc
      subroutine rfexact(kount,nn,ndist, nvar,sscp1,
     *     rec,dat, cova1,means,sd,weight)
cc
cc Determines how many objects lie on the hyperplane with equation
cc z(1,1)*(x_i1 - means_1)+ ... + z(p,1)* (x_ip - means_p) = 0
cc and computes their mean and their covariance matrix.
cc
      double precision ndist(nn)
      double precision sscp1(nvar+1,nvar+1)
      double precision rec(nvar+1)
      double precision dat(nn,nvar)
      double precision cova1(nvar,nvar)
      double precision means(nvar), sd(nvar)
      integer weight(nn)

      call rfcovinit(sscp1,nvar+1,nvar+1)
      kount=0
      do kk=1,nn
         if(dabs(ndist(kk)-0.D0).lt.10.D-8) then
            kount=kount+1
            weight(kk)=1
            do j=1,nvar
               rec(j)=dat(kk,j)
            end do
            call rfadmit(rec,nvar,sscp1)
         else
            weight(kk)=0
         endif
      end do
      call rfcovar(kount,nvar,sscp1,cova1,means,sd)
      return
      end
ccccc
ccccc
      subroutine transfo(cova,means,dat,med,mad,nvar,n)
cc
      implicit none
      integer n, nvar
      double precision dat(n,nvar), cova(nvar,nvar)
      double precision means(nvar), med(nvar), mad(nvar)
      integer i,j,k
      do j=1,nvar
        means(j)=means(j)*mad(j)+med(j)
        do k=1,nvar
           cova(j,k)=cova(j,k)*mad(j)*mad(k)
        end do
        do i=1,n
           dat(i,j)=dat(i,j)*mad(j)+med(j)
        end do
      end do

      return
      end
ccccc
ccccc
      subroutine rfcovmult(a,n1,n2,fac)
cc
cc  Multiplies the matrix a by the real factor fac.
cc
      double precision a(n1,n2)
      double precision fac
cc
      do i=1,n1
        do j=1,n2
          a(i,j)=a(i,j)*fac
        end do
      end do
      return
      end
ccccc
ccccc
      subroutine rfadmit(rec,nvar,sscp)
cc
cc  Updates the sscp matrix with the additional case rec.
cc
      double precision rec(nvar)
      double precision sscp(nvar+1,nvar+1)
cc
      sscp(1,1)=sscp(1,1)+1.D0
      do j=1,nvar
        sscp(1,j+1)=sscp(1,j+1)+rec(j)
        sscp(j+1,1)=sscp(1,j+1)
      end do
      do i=1,nvar
        do j=1,nvar
          sscp(i+1,j+1)=sscp(i+1,j+1)+rec(i)*rec(j)
        end do
      end do
      return
      end
ccccc
ccccc
      subroutine rfcovar(n,nvar, sscp,cova, means,sd)
cc
cc  Computes the classical mean and covariance matrix.
cc
      implicit none
      integer n,nvar, i,j
      double precision sscp(nvar+1,nvar+1), cova(nvar,nvar)
      double precision means(nvar), sd(nvar), f

      do i=1,nvar
        means(i)=sscp(1,i+1)
        sd(i)=sscp(i+1,i+1)
        f=(sd(i)-means(i)*means(i)/n)/(n-1)
        if(f.gt.0.D0) then
          sd(i)=dsqrt(f)
        else
          sd(i)=0.D0
        endif
        means(i)=means(i)/n
      end do
      do i=1,nvar
        do j=1,nvar
          cova(i,j)=sscp(i+1,j+1)
        end do
      end do
      do i=1,nvar
        do j=1,nvar
          cova(i,j)=cova(i,j)-n*means(i)*means(j)
          cova(i,j)=cova(i,j)/(n-1)
        end do
      end do
      return
      end
ccccc
ccccc
      subroutine rfcorrel(nvar,a,b,sd)
cc
cc  Transforms the scatter matrix a to the correlation matrix b: <==> R's  cov2cor(.)
cc
      implicit none
      integer nvar
      double precision a(nvar,nvar), b(nvar,nvar), sd(nvar)
      integer j,i

      do j=1,nvar
         sd(j)=1/sqrt(a(j,j))
      end do
      do i=1,nvar
         do j=1,nvar
            if(i.eq.j) then
               b(i,j)=1.0
            else
               b(i,j)=a(i,j)*sd(i)*sd(j)
            endif
         end do
      end do
      return
      end

      subroutine prdraw(a,pnsel, nn)

      implicit none
      integer nn, a(nn), pnsel
c
      double precision unifrnd
      integer jndex, nrand, i,j

      jndex=pnsel
c     OLD       nrand=int(uniran(seed)*(nn-jndex))+1
      nrand=int(unifrnd() * (nn-jndex))+1
C     if(nrand .gt. nn-jndex) then
C     call intpr(
C     1          '** prdraw(): correcting nrand > nn-jndex; nrand=',
C     2          -1, nrand, 1)
C     nrand=nn-jndex
C     endif

      jndex=jndex+1
      a(jndex)=nrand+jndex-1
      do i=1,jndex-1
         if(a(i).gt.nrand+i-1) then
            do j=jndex,i+1,-1
               a(j)=a(j-1)
            end do
            a(i)=nrand+i-1
            goto 10
c           ------- break
         endif
      end do
 10   continue
      return
      end
ccccc
ccccc
      double precision function rfmahad(rec,nvar,means,sigma)
cc
cc  Computes a Mahalanobis-type distance.
cc
      double precision rec(nvar), means(nvar), sigma(nvar,nvar), t

      t = 0.
      do j=1,nvar
         do k=1,nvar
            t = t + (rec(j)-means(j))*(rec(k)-means(k))*sigma(j,k)
         end do
      end do
      rfmahad=t
      return
      end
ccccc
ccccc
      subroutine rfdis(da,z,ndist,nm,nv,nn,nvar, means)
cc
cc Computes the distance between the objects of da and a hyperplane with
cc equation z(1,1)*(x_i1 - means_1) + ... + z(p,1)*(x_ip - means_p) = 0
cc
      double precision da(nm,nv)
      double precision z(nvar,nvar)
      double precision ndist(nn)
      double precision means(nvar)

      do i=1,nn
         ndist(i)=0
         do j=1,nvar
            ndist(i)=z(j,1)*(da(i,j)-means(j))+ndist(i)
         end do
         ndist(i)=dabs(ndist(i))
      end do
      return
      end
ccccc
ccccc
      subroutine rfstore2(nvar,cstock,mstock,nv_2,
     *     kmini,cova1,means,i,mcdndex,kount)
cc
cc  Stores the coefficients of a hyperplane
cc  z(1,1)*(x_i1 - means_1) + ... +  z(p,1)*(x_ip - means_p) = 0
cc  into the first row of the matrix mstock, and shifts the other
cc  elements of the arrays mstock and cstock.
cc
      double precision cstock(10, nv_2), mstock(10, nvar)
      double precision mcdndex(10, 2, kmini)
      double precision cova1(nvar,nvar), means(nvar)

      do k=10,2,-1
         do kk=1,nvar*nvar
            cstock(k,kk)= cstock(k-1,kk)
         end do
         do kk=1,nvar
            mstock(k,kk)= mstock(k-1,kk)
         end do
         mcdndex(k,1,1)=mcdndex(k-1,1,1)
         mcdndex(k,2,1)=mcdndex(k-1,2,1)
      end do
      do kk=1,nvar
         mstock(1,kk)=means(kk)
         do jj=1,nvar
            cstock(1,(kk-1)*nvar+jj)=cova1(kk,jj)
         end do
      end do
      mcdndex(1,1,1)=i
      mcdndex(1,2,1)=kount
      return
      end
ccccc
ccccc
      subroutine rfstore1(nvar,c1stock,m1stock,nv_2,
     *     kmini,cova1,means,i,km10,ii,mcdndex,kount)

      double precision c1stock(km10, nv_2), m1stock(km10, nvar)
      double precision mcdndex(10,2,kmini)
      double precision cova1(nvar,nvar), means(nvar)

      do k=10,2,-1
         do kk=1,nvar*nvar
            c1stock((ii-1)*10+k,kk)=
     *           c1stock((ii-1)*10+k-1,kk)
         end do
         do kk=1,nvar
            m1stock((ii-1)*10+k,kk)=
     *           m1stock((ii-1)*10+k-1,kk)
         end do
         mcdndex(k,1,ii)=mcdndex(k-1,1,ii)
         mcdndex(k,2,ii)=mcdndex(k-1,2,ii)
      end do
      do kk=1,nvar
         m1stock((ii-1)*10+1,kk)=means(kk)
         do jj=1,nvar
            c1stock((ii-1)*10+1,(kk-1)*nvar+jj)=
     *           cova1(kk,jj)
         end do
      end do
      mcdndex(1,1,ii)=i
      mcdndex(1,2,ii)=kount
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ccccc
ccccc
      subroutine rfcovinit(a,n1,n2)
cc
cc  Initializes the matrix a by filling it with zeroes.
cc
      double precision a(n1,n2)
cc
      do i=1,n1
        do j=1,n2
           a(i,j)=0.D0
        end do
      end do
      return
      end
ccccc
ccccc
      subroutine rfcovsweep(a,nvar,k)
cc
      double precision a(nvar,nvar)
      double precision b, d
cc
      d=a(k,k)
      do j=1,nvar
         a(k,j)=a(k,j)/d
      end do
      do i=1,nvar
         if(i.ne.k) then
            b=a(i,k)
            do j=1,nvar
               a(i,j)=a(i,j)-b*a(k,j)
            end do
            a(i,k) = -b/d
         endif
      end do
      a(k,k)=1/d
      return
      end
ccccc
