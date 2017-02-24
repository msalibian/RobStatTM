c--- For lmrob.lar()  in  ../R/lmrob.M.S.R
c---     ~~~~~~~~~~~
C=======================================================================
      SUBROUTINE rlSTORm2(Y,N,J,YJ)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N)
C-----------------------------------------------------------------------
C     rlSTORm2 SEARCHES THE J-TH VALUE IN ORDER OF MAGNITUDE IN
C     A VECTOR OF LENGTH N.
C-----------------------------------------------------------------------
C--- copied from robust package: src/lmrobmm.f -------------------------
      L=1
      LR=N
 20   IF (L.GE.LR) GOTO 90
      AX=Y(J)
      JNC=L
      JJ=LR
 30   IF(JNC.GT.JJ) GOTO 80
 40   IF (Y(JNC).GE.AX) GOTO 50
      JNC=JNC+1
      GOTO 40
 50   IF(Y(JJ).LE.AX) GOTO 60
      JJ=JJ-1
      GOTO 50
 60   IF(JNC.GT.JJ) GOTO 70
      WA=Y(JNC)
      Y(JNC)=Y(JJ)
      Y(JJ)=WA
      JNC=JNC+1
      JJ=JJ-1
 70   GOTO 30
 80   IF(JJ.LT.J) L=JNC
      IF(J.LT.JNC) LR=JJ
      GOTO 20
 90   YJ=Y(J)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlCOLbi(V1,V2,MLT,M,IOUT)
C.......................................................................
      DOUBLE PRECISION V1(M),V2(M),MLT
C-----------------------------------------------------------------------
C     AUXILIARY ROUTINE FOR rlLARSbi
C-----------------------------------------------------------------------
C--- copied from robust package: src/lmrobbi.f -------------------------
      DO 220 I=1,M
         IF (I .EQ. IOUT) GOTO 220
         V1(I)=V1(I)-V2(I)*MLT
 220  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlICHGbi(A,B)
C.......................................................................
C     AUXILIARY ROUTINE FOR rlLARSbi
C-----------------------------------------------------------------------
C--- copied from robust package: src/lmrobbi.f -------------------------
      DOUBLE PRECISION A,B,C
      C=A
      A=B
      B=C
      RETURN
      END
C=======================================================================
      SUBROUTINE rlLARSbi(X,Y,N,NP,MDX,MDT,TOL,NIT,K,
     +     KODE,SIGMA,THETA,RS,SC1,SC2,SC3,SC4,BET0)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(MDT),RS(N),SC1(N),SC2(NP),
     +     SC3(NP),SC4(NP)
      INTEGER OUT
      LOGICAL STAGE,TEST
      DATA ZERO,TWO,EPS,BIG/0.D0,2.D0,1.0D-10,3.401D38/
cMM   would think rather this: double precision --- but it breaks our checks
C     DATA ZERO,TWO,EPS,BIG/0.D0,2.D0,2.22D-16,1.796D308/

C-----------------------------------------------------------------------
C     LEAST ABSOLUTE RESIDUALS -- aka  L_1 - Regression
C     --> Result in THETA[1:NP]
C-----------------------------------------------------------------------
C---  copied from robust package: src/lmrobbi.f -------------------------
      DO J=1,NP
         SC4(J)=DBLE(J)
         SC2(J)=ZERO
      end do
      SUM=ZERO
      DO I=1,N
         SC1(I)=DBLE(NP+I)
         THETA(I)=Y(I)
         IF (Y(I) .lt. ZERO) then
            DO J=1,NP
               X(I,J)=-X(I,J)
            end do
            THETA(I)=-THETA(I)
            SC1(I)=-SC1(I)
         endif
         SUM=SUM+THETA(I)
      end do
C-----------------------------------------------------------------------
C     COMPUTE THE MARGINAL COSTS.
C-----------------------------------------------------------------------
      SUMIN=SUM
      DO J=1,NP
         SUM=ZERO
         DO I=1,N
            SUM=SUM+X(I,J)
         end do
         SC3(J)=SUM
      end do
C-----------------------------------------------------------------------
C     STAGE I. DETERMINE THE VECTOR TO ENTER THE BASIS.
C-----------------------------------------------------------------------
      TEST=.FALSE.              ! -Wall
      STAGE=.TRUE.
      KOUNT=0
      KR=1
      KL=1
      IN=1                      ! -Wall

c--   ---------------- LOOP (Stage I) ------------------------------------
 70   VMAX=-1.D0
      DNP=DBLE(NP)
      DO J=KR,NP
         IF (DABS(SC4(J)) .GT. DNP) cycle ! = continue
         D=DABS(SC3(J))
         IF (D-VMAX .LE. ZERO) cycle
         IF (D-VMAX .LE. EPS)  cycle
         VMAX=D
         IN=J
      end do
      IF (SC3(IN) .lt. ZERO) then ! swap signs
         do I=1,N
            X(I,IN)=-X(I,IN)
         end do
         SC3(IN)=-SC3(IN)
         SC4(IN)=-SC4(IN)
      endif

C-----------------------------------------------------------------------
C     DETERMINE THE VECTOR TO LEAVE THE BASIS.
C-----------------------------------------------------------------------

cvvv  ------------ 2nd-level loop ---------------------------------
 100  K=0
      DO I=KL,N
         D=X(I,IN)
         IF (D .LE. TOL) cycle
         K=K+1
         Y(K)=THETA(I)/D
         RS(K)=DBLE(I)
         TEST=.TRUE.
      end do

C---  -------------- 3rd-level loop ------------------
 120  IF (K .le. 0) then
         TEST=.FALSE.           ! and GOTO 150
      else                      ! 130
         VMIN=BIG
         DO I=1,K
            IF (Y(I)-VMIN .GE. ZERO) cycle
            IF (VMIN-Y(I) .LE. EPS)  cycle
            J=I
            VMIN=Y(I)
            OUT=INT(RS(I))
         end do
         Y(J)=Y(K)
         RS(J)=RS(K)
         K=K-1
      endif

C-----------------------------------------------------------------------
C     CHECK FOR LINEAR DEPENDENCE IN STAGE I.
C-----------------------------------------------------------------------
c     150
      IF (.not.TEST .and. STAGE) then
         DO I=1,N
            CALL rlICHGbi(X(I,KR),X(I,IN))
         end do
         CALL rlICHGbi(SC3(KR),SC3(IN))
         CALL rlICHGbi(SC4(KR),SC4(IN))
         KR=KR+1
c     GOTO 260
      else
c     170
         IF (.not. TEST) then
            KODE=2
            GOTO 350
         endif
c     180
         PIVOT=X(OUT,IN)
         IF (SC3(IN)-PIVOT-PIVOT .gt. TOL) then ! not converged
            DO J=KR,NP
               D=X(OUT,J)
               SC3(J)=SC3(J)-D-D
               X(OUT,J)=-D
            end do
            D=THETA(OUT)
            SUMIN=SUMIN-D-D
            THETA(OUT)=-D
            SC1(OUT)=-SC1(OUT)
            GOTO 120
c     -----------end{ 3rd-level loop } -----------------
         endif

C-----------------------------------------------------------------------
C 200   PIVOT ON X(OUT,IN).
C-----------------------------------------------------------------------
         DO J=KR,NP
            IF (J.EQ.IN) cycle  ! = continue
            X(OUT,J)=X(OUT,J)/PIVOT
         end do
         THETA(OUT)=THETA(OUT)/PIVOT
         DO J=KR,NP
            IF (J .EQ. IN) cycle
            D=X(OUT,J)
            SC3(J)=SC3(J)-D*SC3(IN)
            CALL rlCOLbi(X(1,J),X(1,IN),D,N,OUT)
         end do
         SUMIN=SUMIN-SC3(IN)*THETA(OUT)
         DO I=1,N
            IF (I .EQ. OUT) cycle
            D=X(I,IN)
            THETA(I)=THETA(I)-D*THETA(OUT)
            X(I,IN)=-D/PIVOT
         end do
         SC3(IN)=-SC3(IN)/PIVOT
         X(OUT,IN)=1.D0/PIVOT
         CALL rlICHGbi(SC1(OUT),SC4(IN))
         KOUNT=KOUNT+1
         IF (.NOT. STAGE) GOTO 270
C-----------------------------------------------------------------------
C     INTERCHANGE ROWS IN STAGE I.
C-----------------------------------------------------------------------
         KL=KL+1
         DO J=KR,NP
            CALL rlICHGbi(X(OUT,J),X(KOUNT,J))
         enddo
         CALL rlICHGbi(THETA(OUT),THETA(KOUNT))
         CALL rlICHGbi(SC1(OUT),SC1(KOUNT))
      endif

      IF (KOUNT+KR .NE. NP+1) GOTO 70
c                             =======

C-----------------------------------------------------------------------
C     STAGE II. DETERMINE THE VECTOR TO ENTER THE BASIS.
C-----------------------------------------------------------------------
      STAGE=.FALSE.
cvvv
 270  VMAX=-BIG
      DO J=KR,NP
         D=SC3(J)
         IF (D .lt. ZERO) then
            IF (D+TWO .GT. ZERO) cycle
            D=-D-TWO
         endif
         IF (D-VMAX .LE. ZERO) cycle
         IF (D-VMAX .LE. EPS)  cycle
         VMAX=D
         IN=J
      end do
      IF (VMAX .gt. TOL) then   ! not converged
         IF (SC3(IN) .le. ZERO) then
            DO I=1,N
               X(I,IN)=-X(I,IN)
            end do
            SC3(IN)=-SC3(IN)-2.D0
            SC4(IN)=-SC4(IN)
         endif
         GOTO 100
c        ========
      endif
C-----------------------------------------------------------------------
C 310  PREPARE OUTPUT
C-----------------------------------------------------------------------
      L=KL-1
      DO I=1,N
         RS(I)=ZERO
         IF (I .GT. L .OR. THETA(I) .GE. ZERO) cycle
         do J=KR,NP
            X(I,J)=-X(I,J)
         end do
         THETA(I)=-THETA(I)
         SC1(I)=-SC1(I)
      end do
      KODE=0
      IF (KR .eq. 1) then       ! first time only
         do J=1,NP
            D=DABS(SC3(J))
            IF (D .LE. TOL .OR. TWO-D .LE. TOL) GOTO 350
         end do
         KODE=1
      endif
c---
 350  DO I=1,N
         K=INT(SC1(I))
         D=THETA(I)
         IF (K .le. 0) then
            K=-K
            D=-D
         endif
         IF (I .lt. KL) then
            SC2(K)=D
         else
            K=K-NP
            RS(K)=D
         endif
      end do
      K=NP+1-KR
      SUM=ZERO
      DO I=KL,N
         SUM=SUM+THETA(I)
      end do
      SUMIN=SUM
      NIT=KOUNT
      DO J=1,NP
         THETA(J)=SC2(J)
      end do
      DO I=1,N
         Y(I)=DABS(RS(I))
      end do
      N2=N/2+1
      CALL RLSTORM2(Y,N,N2,SIGMA)
      SIGMA=SIGMA/BET0
      RETURN
      END
