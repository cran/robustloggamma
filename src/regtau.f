
      SUBROUTINE RICLLS(XT,Y,N,NP,MDX,MDT,K,IX,IY,SIGMA,THETA,RS1,RS2,
     1                SE,SF,SG,SH,IP)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR: A. MARAZZI
C.......................................................................
C
      DIMENSION XT(MDX,NP),Y(N),THETA(MDT),RS1(N),RS2(N),
     1          SE(NP),SF(NP),SG(NP),SH(NP)
      INTEGER IP(NP)
      LOGICAL NPRCHK
C
C  PARAMETER CHECK AND INITIALIZATION
C
      KP1=K+1
      MDXP1=MDX+1
      KK=MDX*(K-1)+K
      LDIAG=MIN0(N,NP)
      MX=MAX0(N,NP)
      NPRCHK=LDIAG.GT.0.AND.N.LE.MDX.AND.MDT.GE.MX
     1       .AND.((IX.EQ.0.AND.K.GT.0.AND.K.LE.LDIAG).OR.
     2       IX.EQ.1).AND.(IY.EQ.0.OR.IY.EQ.1)
      IF (.NOT.NPRCHK) CALL MESSGE(500,'RICLLS',1)
      IF (IX.EQ.1) CALL RIMTRF(XT,N,NP,MDX,0,1.E-6,K,SF,SG,SH,IP)
C
C  HOUSEHOLDER TRANSFORMATIONS OF THE OBSERVATION VECTOR
C
      IF (IY.EQ.0) GOTO 25
      IF (K.NE.NP) CALL SWAP(XT,SF,K,MDXP1,1,KK,K)
      DO 20 JJ=1,LDIAG
      J=JJ
   20 CALL H12(2,J,J+1,N,XT(1,J),1,SH(J),Y,1,N,1,N)
      IF (K.NE.NP) CALL SWAP(XT,SF,K,MDXP1,1,KK,K)
C
C  SOLVE THE TRANSFORMED LS-PROBLEM
C
   25 DO 30 I=1,N
   30 THETA(I)=Y(I)
      CALL SOLV(XT,THETA,NP,K,MDX,MDT)
      IF (K.EQ.NP) GOTO 50
      DO 40 J=KP1,NP
   40 THETA(J)=0.0
   50 CONTINUE
C
C  COMPUTE THE TRANSFORMED RESIDUAL VECTORS Q*(Y-X1*THETA) (RS1)
C  AND (IF K.LT.NP) Q*(Y-X*THETA) (RS2)
C
      IF (K.EQ.NP) GOTO 60
      CALL RES(3,XT,Y,THETA,RS2,SE,SG,N,NP,K,NP,MDX,MDT)
   60 CALL RES(1,XT,Y,THETA,RS1,SE,SG,N,NP,K,NP,MDX,MDT)
C
C  COMPUTE SIGMA
C
      SIGMA=0.
      IF (K.EQ.N) GOTO 80
      CALL NRM2(RS1,N,1,N,SIGMA)
      SIGMA=SIGMA/SQRT(FLOAT(N-K))
   80 IF (SIGMA.GT.0.) GOTO 90
      CALL MESSGE(101,'RICLLS',0)
   90 CONTINUE
C
C  TRANSFORM RESIDUAL VECTORS FOR OUTPUT
C
      IF (K.NE.NP) CALL SWAP(XT,SF,K,MDXP1,1,KK,K)
      DO 100 J1=1,LDIAG
      J=LDIAG-J1+1
      IF (K.EQ.NP) GOTO 100
      CALL H12(2,J,J+1,N,XT(1,J),1,SH(J),RS2,1,N,1,N)
  100 CALL H12(2,J,J+1,N,XT(1,J),1,SH(J),RS1,1,N,1,N)
      IF (K.NE.NP) CALL SWAP(XT,SF,K,MDXP1,1,KK,K)
C
C  TRANSFORM SOLUTION VECTOR FOR OUTPUT
C
      IF (K.EQ.NP) GOTO 120
      DO 110 J=1,K
      I=J
  110 CALL H12(2,I,KP1,NP,XT(I,1),MDX,SG(I),THETA,1,N,1,NP)
  120 CALL PERM(THETA,IP,LDIAG,NP)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RIMTRF(X,N,NP,MDX,INTCH,TAU,K,SF,SG,SH,IP)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR: A. MARAZZI
C.......................................................................
C
      DIMENSION X(MDX,NP),SF(NP),SG(NP),SH(NP)
      INTEGER IP(NP)
      LOGICAL NPRCHK
      EXTERNAL DIFF
C
C  PARAMETER CHECK AND INITIALIZATION
C
      FACTOR=0.001
      LDIAG=MIN0(N,NP)
      NPRCHK=LDIAG.GT.0.AND.N.LE.MDX.AND.TAU.GE.0.
      IF (.NOT.NPRCHK) CALL MESSGE(500,'RIMTRF',1)
C
      DO 80 JJ=1,LDIAG
      J=JJ
      IF (INTCH.EQ.0) IP(J)=J
      IF (INTCH.EQ.0) GOTO 70
      IF (J.EQ.1) GOTO 20
C
C  UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
C
      LMAX=J
      DO 10 L=J,NP
      SH(L)=SH(L)-X(J-1,L)**2
      IF(SH(L).GT.SH(LMAX)) LMAX=L
   10 CONTINUE
CC Modified by C. Agostinelli 2012-09-04
CC      IF (DIFF(HMAX+FACTOR*SH(LMAX),HMAX)) 20,20,50
      IF (DIFF(HMAX+FACTOR*SH(LMAX),HMAX).GT.0) GOTO 50
C
C  COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
C
   20 LMAX=J
      DO 40 L=J,NP
      SH(L)=0.
      DO 30 I=J,N
   30 SH(L)=SH(L)+X(I,L)**2
      IF (SH(L).GT.SH(LMAX)) LMAX=L
   40 CONTINUE
      HMAX=SH(LMAX)
C
C  LMAX HAS BEEN DETERMINED: INTERCHANGE COLUMNS IF NEEDED
C
   50 CONTINUE
      IP(J)=LMAX
      IF (IP(J).EQ.J) GOTO 70
      DO 60 I=1,N
      TMP=X(I,J)
      X(I,J)=X(I,LMAX)
   60 X(I,LMAX)=TMP
      SH(LMAX)=SH(J)
C
C  COMPUTE THE HOUSEHOLDER TRANSF. Q AND APPLY IT TO X
C
   70 MDC=NP-J
      IF (MDC.GT.0)
     1CALL H12(1,J,J+1,N,X(1,J),1,SH(J),X(1,J+1),1,MDX,MDC,MDX*MDC)
      IF (MDC.EQ.0)
     1CALL H12(1,J,J+1,N,X(1,J),1,SH(J),SF,1,1,0,1)
   80 CONTINUE
C
C  X CONTAINS NOW THE TRANSFORMED DESIGN MATRIX Q*X.
C  DETERMINE THE PSEUDORANK K USING THE TOLERANCE TAU
C
      DO 100 J=1,LDIAG
      IF (ABS(X(J,J)).LE.TAU) GOTO 110
  100 CONTINUE
      K=LDIAG
      GOTO 120
  110 K=J-1
  120 KP1=K+1
C
C  IF THE PSEUDORANK IS LESS THAN NP STORE THE FIRST K
C  DIAG.ELEMENTS OF X FOR FURTHER APPLICATIONS OF Q
C
      IF (K.EQ.NP) GOTO 130
      DO 125 I=1,K
  125 SF(I)=X(I,I)
  130 CONTINUE
C
C  SPECIAL FOR PSEUDORANK=0
C
      IF (K.GT.0) GOTO 140
      CALL MESSGE(401,'RIMTRF',0)
      RETURN
C
C  IF THE PSEUDORANK IS LESS THAN NP COMPUTE Q*X*V
C
  140 IF (K.EQ.NP) GOTO 160
      MDC=MDX*(NP-1)
      DO 150 II=1,K
      I=KP1-II
  150 CALL H12(1,I,KP1,NP,X(I,1),MDX,SG(I),X,MDX,1,I-1,MDC+I-1)
  160 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SWAP(X,Y,N,INCX,INCY,MDX,MDY)
C.......................................................................
C
C   COPYRIGHT 1979 SOCIETY FOR INDUSTRIAL AND APPLIED MATHEMATICS.
C   ALL RIGHTS RESERVED.
C
C   AUTHOR :     LINPACK (SUBROUTINE SSWAP)
C                REPRINTED WITH PERMISSION FROM 
C                LINPACK USER'S GUIDE.
C                ADAPTED FOR ROBETH BY A. MARAZZI
C.......................................................................
C
      REAL X(MDX),Y(MDY)
      LOGICAL NPRCHK
C
C  PARAMETER CHECK
C
      NPRCHK=N.GE.0.AND.INCX.NE.0.AND.IABS(INCX)*(N-1)+1.LE.MDX
     1       .AND.INCY.NE.0.AND.IABS(INCY)*(N-1)+1.LE.MDY
      IF (.NOT.NPRCHK) CALL MESSGE(500,'SWAP  ',1)
C
      IF (N.EQ.0) RETURN
      IF (INCX.EQ.1.AND.INCY.EQ.1) GOTO 20
C
C  CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT
C  EQUAL TO 1
C
      IX=1
      IY=1
      IF (INCX.LT.0) IX=(-N+1)*INCX+1
      IF (INCY.LT.0) IY=(-N+1)*INCY+1
      DO 10 I=1,N
      TEMP=X(IX)
      X(IX)=Y(IY)
      Y(IY)=TEMP
      IX=IX+INCX
      IY=IY+INCY
   10 CONTINUE
      RETURN
C
C  CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 M=MOD(N,3)
      IF (M.EQ.0) GOTO 40
      DO 30 I=1,M
      TEMP=X(I)
      X(I)=Y(I)
      Y(I)=TEMP
   30 CONTINUE
      IF (N.LT.3) RETURN
   40 MP1=M+1
      DO 50 I=MP1,N,3
      TEMP=X(I)
      X(I)=Y(I)
      Y(I)=TEMP
      TEMP=X(I+1)
      X(I+1)=Y(I+1)
      Y(I+1)=TEMP
      TEMP=X(I+2)
      X(I+2)=Y(I+2)
      Y(I+2)=TEMP
   50 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SOLV(X,THETA,NP,K,MDX,MDT)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
C  PURPOSE
C  -------
C  LET U BE THE K BY K UPPER TRIANGULAR MATRIX WITH ELEMENTS
C  X(1,1)...X(1,K),X(2,2)...X(2,K),...X(K,K).
C  SOLV SOLVES THE TRIANGULAR SYSTEM U*THETA=Y (BY BACK SUBSTI-
C  TUTION). ON INPUT Y IS CONTAINED IN THETA.  ON OUTPUT
C  THETA(1)...THETA(K) CONTAIN THE DESIRED SOLUTION.
C
C  ERRORS
C  ------
C  1   AN ELEMENT OF THE PRINCIPAL DIAGONAL OF X IS =0.
C
      DIMENSION X(MDX,NP),THETA(MDT)
      DOUBLE PRECISION SM,DZERO
      DZERO=0.D0
      KP1=K+1
      DO 80 L=1,K
      SM=DZERO
      I=KP1-L
      IF (I.EQ.K) GOTO 60
      IP1=I+1
      DO 50 J=IP1,K
   50 SM=SM+X(I,J)*DBLE(THETA(J))
   60 SM1=SM
CC Modified by C. Agostinelli 2012-09-04
CC    IF (X(I,I)) 80,70,80
      IF (X(I,I).EQ.0) THEN
   70 CALL MESSGE(501,'SOLV  ',1)
      ENDIF
   80 THETA(I)=(THETA(I)-SM1)/X(I,I)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RES(MODE,X,Y,S,R,COV,SG,N,NP,K,NCOV,MDX,MDS)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
C  PURPOSE
C  -------
C  THE ARRAY X CONTAINS THE TRANSFORMED DESIGN MATRIX.
C  S(1)...S(K) CONTAIN A SOLUTION OF THE TRANSFORMED
C  LS-PROBLEM AND Y, THE TRANSFORMED RIGHT SIDE.
C  LET U BE THE K BY K UPPER TRIANGULAR MATRIX WITH ELEMENTS
C  X(1,1)...X(1,K),X(2,2)...X(2,K),...,X(K,K).
C  IF MODE=1 OR 2 RES SETS (R(1)...R(K)=(Y(1)...Y(K))
C  -U*(S(1)...S(K)),(R(K+1)...R(N))=(Y(K+1)...Y(N))
C  (SET MODE=1 IF (S(1)...S(K))=U**(-1)*(Y(1)...Y(K)))
C  IF MODE=3 RES COMPUTES THE RESIDUALS OF THE ORIGINAL
C  LS-PROBLEM AFTER THE APPLICATION OF Q.
C  SG CONTAINS PIVOT SCALARS FOR HOUSEHOLDER TRANSFORMATIONS.
C  COV IS USED AS SCRATCH (ONLY THE LAST NP LOCATIONS AND
C  ONLY IF K.LT.NP)
C
      DIMENSION X(MDX,NP),Y(N),R(N),S(MDS),COV(NCOV),SG(NP)
      DOUBLE PRECISION SM,DZERO
      LDIAG=MIN0(N,NP)
      DZERO=0.D0
C
C  COMPUTE THE RESIDUALS R(1)...R(K).
C
      IF (MODE.EQ.2.OR.MODE.EQ.3) GOTO 20
      DO 10 I=1,K
   10 R(I)=0.
      GOTO 50
   20 CONTINUE
      IF (K.LT.N) GOTO 25
      DO 26 I=1,N
   26 R(I)=0.
      GOTO 130
   25 CONTINUE
      DO 40 I=1,K
      SM=DZERO
      DO 30 J=I,K
   30 SM=SM+X(I,J)*DBLE(S(J))
      SM1=SM
   40 R(I)=Y(I)-SM1
   50 CONTINUE
C
C  COMPUTE THE RESIDUALS R(K+1)...R(NP)
C
      IF (K.EQ.NP) GOTO 110
      IF (K.EQ.N) GOTO 130
      KP1=K+1
      IF (MODE.EQ.3) GOTO 55
      DO 52 I=KP1,LDIAG
   52 R(I)=Y(I)
      GOTO 110
   55 INZ=NCOV-NP
      DO 100 I=KP1,LDIAG
      IM1=I-1
      DO 60 J=1,IM1
      J1=J+INZ
   60 COV(J1)=0.
      DO 70 J=I,NP
      J1=J+INZ
   70 COV(J1)=X(I,J)
      DO 80 II=1,K
      I1=KP1-II
   80 CALL R3V(I1,KP1,NP,X(I1,1),MDX,SG(I1),COV,NCOV,INZ)
      SM=DZERO
      DO 90 J=1,K
      J1=J+INZ
   90 SM=SM+COV(J1)*DBLE(S(J))
      SM1=SM
      R(I)=Y(I)-SM1
  100 CONTINUE
C
C  COMPUTE THE RESIDUALS R(NP+1)...R(N)
C
  110 NPP1=NP+1
      IF (NP.GE.N) GOTO 130
      DO 120 I=NPP1,N
  120 R(I)=Y(I)
  130 RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE R3V(LPIVOT,L1,M,U,IUE,UP,C,NCOV,INZ)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
C  PURPOSE
C  -------
C  R3V APPLIES THE HOUSEHOLDER TRANSFORMATION DEFINED
C  BY THE VECTOR U TO A VECTOR E=(0,...,0,E(I+1),...,E(N))
C  WHERE I.GE.LPIVOT.
C  E IS STORED IN THE LOCATIONS INZ+1,...,INZ+N OF C.
C  PARAMETERS LPIVOT,L1,M,U,IUE,UP AS IN H12.
C
C
      DIMENSION U(IUE,M),C(NCOV)
      DOUBLE PRECISION SM,B,DZERO
      ONE=1.
      DZERO=0.D0
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN
      CL=ABS(U(1,LPIVOT))
CC Modified by C. Agostinelli 2012-09-04
CC      IF (CL) 130,130,70
      IF (CL.LE.0) GOTO 130
   70 CONTINUE
      B=DBLE(UP)*U(1,LPIVOT)
C
C  B MUST BE NONPOSITIVE HERE. IF B=0., RETURN.
C
CC Modified by C. Agostinelli 2012-09-04
CC      IF (B) 80,130,130
      IF (B.GE.0) GOTO 130
   80 B=ONE/B
      I2=LPIVOT-1+INZ
      INCR=L1-LPIVOT
      I2=I2+1
      I3=I2+INCR
      I4=I3
      SM=DZERO
      DO 90 I=L1,M
      SM=SM+C(I3)*DBLE(U(1,I))
   90 I3=I3+1
CC Modified by C. Agostinelli 2012-09-04
CC      IF (SM) 100,120,100
      IF (SM.EQ.0) GOTO 120
  100 SM=SM*B
      C(I2)=C(I2)+SM*DBLE(UP)
      DO 110 I=L1,M
      C(I4)=C(I4)+SM*DBLE(U(1,I))
  110 I4=I4+1
  120 CONTINUE
  130 RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION DIFF(X,Y)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
C  PURPOSE
C  -------
C  SEE (L6) IN H12, P.278.
C
      DIFF=X-Y
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NRM2(X,N,INCX,MDX,XNRM)
C.......................................................................
C
C   COPYRIGHT 1979 SOCIETY FOR INDUSTRIAL AND APPLIED MATHEMATICS.
C   ALL RIGHTS RESERVED.
C
C   AUTHOR :     LINPACK (SUBROUTINE SNRM2)
C                REPRINTED WITH PERMISSION FROM 
C                LINPACK USER'S GUIDE.
C                ADAPTED FOR ROBETH BY A. MARAZZI
C.......................................................................
C
      LOGICAL NPRCHK
      REAL X(MDX)
      DOUBLE PRECISION SUM,ZERO,ONE,XMAX,DXI
C     DATA ZERO,ONE/0.0D0,1.0D0/,CUTLO,CUTHI/4.441E-16,1.304E19/
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,CUTLO=4.441E-16,CUTHI=1.304E19)
C
C  PARAMETER CHECK
C
      NPRCHK=INCX.GT.0.AND.INCX*(N-1)+1.LE.MDX
      IF (.NOT.NPRCHK) CALL MESSGE(500,'NRM2  ',1)
C
      IF (N.GT.0) GOTO 10
      XNRM=0.
      GOTO 300
C
   10 NEXT=1  !ASSIGN 30 TO NEXT
      SUM=ZERO
      NN=N*INCX
C
C  BEGIN MAIN LOOP
C
      I=1
   20 DXI=DBLE(X(I))
      GOTO (30,50,70,110) NEXT
   30 IF (ABS(X(I)).GT.CUTLO) GOTO 85
      NEXT=2  !ASSIGN 50 TO NEXT
      XMAX=ZERO
C
C  PHASE1.  SUM IS ZERO
C
   50 IF (X(I).EQ.0.) GOTO 200
      IF (ABS(X(I)).GT.CUTLO) GOTO 85
C
C  PREPARE FOR PHASE 2.
C
      NEXT=3   !ASSIGN 70 TO NEXT
      GOTO 105
C
C  PREPARE FOR PHASE 4.
C
  100 I=J
      DXI=DBLE(X(I))
      NEXT=4   !ASSIGN 110 TO NEXT
      SUM=(SUM/DXI)/DXI
  105 XMAX=DABS(DXI)
      GOTO 115
C
C  PHASE 2.  SUM IS SMALL. SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF (ABS(X(I)).GT.CUTLO) GOTO 75
C
C  COMMON CODE FOR PHASE 2 AND 4.
C  IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF (DABS(DXI).LE.XMAX) GOTO 115
      SUM=ONE+SUM*(XMAX/DXI)**2
      XMAX=DABS(DXI)
      GOTO 200
C
  115 SUM=SUM+(DXI/XMAX)**2
      GOTO 200
C
C  PREPARE FOR PHASE 3.
C
   75 SUM=(SUM*XMAX)*XMAX
C
C  SET HITEST=CUTHI/N
C
   85 HITEST=CUTHI/FLOAT(N)
C
C  PHASE3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J=I,NN,INCX
      IF (ABS(X(J)).GE.HITEST) GOTO 100
   95 SUM=SUM+X(J)*DBLE(X(J))
      XNRM=DSQRT(SUM)
      GOTO 300
C
  200 CONTINUE
      I=I+INCX
      IF (I.LE.NN) GOTO 20
C
C  END MAIN LOOP
C
      XNRM=XMAX*DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PERM(X,SP,N,NDIM)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
C  PURPOSE
C  -------
C  PERMUTE COMPONENTS OF X TO COMPENSATE COLUMN INTERCH.
C
      DIMENSION X(NDIM)
      INTEGER SP(NDIM)
      DO 10 JJ=1,N
      J=N-JJ+1
      IF (SP(J).EQ.J) GOTO 10
      L=SP(J)
      TMP=X(L)
      X(L)=X(J)
      X(J)=TMP
   10 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE H12(MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV,
     1               MDC)
C.......................................................................
C
C   AUTHORS :     CH.L. LAWSON & R.J. HANSON (1974)
C                 SOLVING LEAST SQUARES PROBLEMS 
C                 REPRINT FROM PP.290-291,308 BY PERMISSION OF 
C                 PRENTICE HALL, ENGLEWOOD CLIFFS, NEW JERSEY.
C                 ADAPTED FOR ROBETH BY A. MARAZZI
C.......................................................................
C
      REAL U(IUE,M),C(MDC)
      DOUBLE PRECISION SM,B
      ONE=1.
C
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN
      CL=ABS(U(1,LPIVOT))
      IF (MODE.EQ.2) GOTO 60
C
C  CONSTRUCT THE TRANSFORMATION
C
      DO 10 J=L1,M
   10 CL=AMAX1(ABS(U(1,J)),CL)
CC Modified by C. Agostinelli 2012-09-04
CC      IF (CL) 130,130,20
      IF (CL.LE.0) GOTO 130
   20 CLINV=ONE/CL
      SM=(DBLE(U(1,LPIVOT))*CLINV)**2
      DO 30 J=L1,M
   30 SM=SM+(DBLE(U(1,J))*CLINV)**2
C
C  CONVERT DBLE. PRE. SM TO SNGL. PREC. SM1
C
      SM1=SM
      CL=CL*SQRT(SM1)
CC Modified by C. Agostinelli 2012-09-04
CC      IF (U(1,LPIVOT)) 50,50,40
      IF (U(1,LPIVOT).GT.0) THEN
   40 CL=-CL
      ENDIF
   50 UP=U(1,LPIVOT)-CL
      U(1,LPIVOT)=CL
      GOTO 70
C
C  APPLY THE TRANSFORMATION I+U*(U**T)/B TO C
C
CC Modified by C. Agostinelli 2012-09-04
CC   60 IF (CL) 130,130,70
   60 IF (CL.LE.0) GOTO 130
   70 IF (NCV.LE.0) RETURN
      B=DBLE(UP)*U(1,LPIVOT)
C
C  B MUST BE NONPOSITIVE HERE. IF B=0., RETURN.
C
CC Modified by C. Agostinelli 2012-09-04
CC      IF (B) 80,130,130
      IF (B.GE.0) GOTO 130
   80 B=ONE/B
      I2=1-ICV+ICE*(LPIVOT-1)
      INCR=ICE*(L1-LPIVOT)
      DO 120 J=1,NCV
      I2=I2+ICV
      I3=I2+INCR
      I4=I3
      SM=C(I2)*DBLE(UP)
      DO 90 I=L1,M
      SM=SM+C(I3)*DBLE(U(1,I))
   90 I3=I3+ICE
CC Modified by C. Agostinelli 2012-09-04
CC      IF (SM) 100,120,100
      IF (SM.EQ.0) GOTO 120
  100 SM=SM*B
      C(I2)=C(I2)+SM*DBLE(UP)
      DO 110 I=L1,M
      C(I4)=C(I4)+SM*DBLE(U(1,I))
  110 I4=I4+ICE
  120 CONTINUE
  130 RETURN
      END
C
C=======================================================================
C
      SUBROUTINE REGTAU(X,Y,N,NN,B1,C1,B2,C2,TOL,ISEED,AO,BO,TO,
     +           RS,SA,SB,TA,TMP1,TMP2)
C.......................................................................
C
C   COPYRIGHT 2009 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION X(N),Y(N),RS(N),SA(NN),SB(NN),TA(NN),TAI
      DOUBLE PRECISION AO,BO,TO,SRJ,SX(2),SY(2),SUMXY,SUMXX,SUMX,SUMY
      REAL TMP1(N),TMP2(N)
      INTEGER IND(2)
      LOGICAL NPRCHK
      EXTERNAL RHO4
      COMMON/RHOTAU/IPSI,C,H1,H2,H3,XK,D

C  ADDED CLAUDIO 16/08/2012
C  TO USE THE UNIFORM RANDOM NUMBER GENERATOR FROM R 
C
      external rndstart
      external rndend
      external rndunif

      DOUBLE PRECISION rndunif

      call rndstart()

C
C  PARAMETER CHECK AND INITIALIZATION
C
      NPRCHK=N.GT.0.AND.TOL.GT.0.AND.ISEED.NE.0
      IF (.NOT.NPRCHK) CALL MESSGE(500,'REGTAU',1)
      N1=N/2
      IPS0=IPSI
      XK0=XK
      IPSI=4
      IO=N
      TAI=1.D6
      KSEED=ISEED
      DO 500 I=1,NN
c       call intpr('i',1,i,1)
        DO 30 K=1,2

C  MODIFIED CLAUDIO 16/08/2012
C  TO USE THE UNIFORM RANDOM NUMBER GENERATOR FROM R 
C

CCCCC   10     CALL RANDOW(KSEED,RND)
 10       RND = rndunif()
          IK=INT(RND*FLOAT(N))+1
          IF (IK.GT.N) IK=N
          IF (K.EQ.1) THEN
            IND(1)=IK
            IND1=IK
            GOTO 30
          ENDIF
          IF (IK.EQ.IND1) GOTO 10
          IF (DABS(X(IND1)-X(IK)).LE.1.D-5) GOTO 10
          IND(2)=IK
   30   CONTINUE
      DO 100 K=1,2
      IK=IND(K)
      SX(K)=X(IK)
  100 SY(K)=Y(IK)
      SB(I)=(SY(2)-SY(1))/(SX(2)-SX(1))
      SA(I)=SY(1)-SB(I)*SX(1)
      DO 120 J=1,N
      RS(J)=Y(J)-SB(I)*X(J)-SA(I)
      SRJ=DABS(RS(J))
      TMP1(J)=SNGL(SRJ)
      TMP2(J)=FLOAT(J)
  120 CONTINUE
      CALL SRT2(TMP1,TMP2,N,1,N)
      SUMXX=0.D0
      SUMXY=0.D0
      SUMX=0.D0
      SUMY=0.D0
      DO 140 J=1,N1
      K=TMP2(J)
      SUMXX=SUMXX+X(K)*X(K)
      SUMXY=SUMXY+X(K)*Y(K)
      SUMX=SUMX+X(K)
      SUMY=SUMY+Y(K)
  140 CONTINUE 
      SB(I)=(SUMXY-SUMX*SUMY/DFLOAT(N1))/(SUMXX-SUMX*SUMX/DFLOAT(N1))
      SA(I)=(SUMY-SB(I)*SUMX)/DFLOAT(N1)
      DO 160 J=1,N
      RS(J)=Y(J)-SB(I)*X(J)-SA(I)
      SRJ=DABS(RS(J))
      TMP1(J)=SNGL(SRJ)
  160 CONTINUE
      CALL SRT1(TMP1,N,1,N)
      S0=TMP1(N1+1)
      IF (2*N1.EQ.N) S0=(TMP1(N1)+S0)/2.0
      S0=S0/0.6745
      TAU=TOL
      IF (S0.LE.TOL) GOTO 300
      IT=0
      H=1.0
      XK=C1
c      call realpr('s0',2,s0,1)
  180 IT=IT+1
      SUM=0.0
C     RHO4(S)=RHO4(-S) !!
      DO 200 J=1,N
      U=TMP1(J)
  200 SUM=SUM+RHO4(U/S0)
      S1=SUM/(FLOAT(N)*B1)
      S1=S0*SQRT(S1)
      H=ABS(S1-S0)/S0
      S0=S1
      IF (H.GT.TOL.AND.IT.LT.50) GOTO 180
c      call realpr('h',1,h,1)
c      call intpr('it',2,it,1)
  300 CONTINUE
      IF (S0.LE.TOL) GOTO 400
      SUM=0.0
      XK=C2
      DO 320 J=1,N
      U=TMP1(J)
  320 SUM=SUM+RHO4(U/S0)
      S1=SUM/(FLOAT(N)*B2)
      TAU=S0*SQRT(S1)      
  400 TA(I)=DBLE(TAU)
      IF (TA(I).GE.TAI) GOTO 500
      IO=I
      TAI=TA(I)
  500 CONTINUE
      IPSI=IPS0
      XK=XK0
      AO=SA(IO)
      BO=SB(IO)
      TO=TA(IO)

C  ADDED CLAUDIO 16/08/2012
C  TO USE THE UNIFORM RANDOM NUMBER GENERATOR FROM R 
C
      call rndend()
      RETURN
      END
C
      SUBROUTINE REGTAUW(X,Y,W,N,NN,B1,C1,B2,C2,TOL,ISEED,AO,BO,TO,
     +           RS,SA,SB,TA,TMP1,TMP2)
C.......................................................................
C
C   COPYRIGHT 2009 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION X(N),Y(N),W(N),RS(N),SA(NN),SB(NN),TA(NN),TAI
      DOUBLE PRECISION AO,BO,TO,SRJ,SX(2),SY(2),SUMXY,SUMXX,SUMX,SUMY
      REAL TMP1(N),TMP2(N)
      INTEGER IND(2)
      LOGICAL NPRCHK
      EXTERNAL RHO4
      COMMON/RHOTAU/IPSI,C,H1,H2,H3,XK,D

C  ADDED CLAUDIO 16/08/2012
C  TO USE THE UNIFORM RANDOM NUMBER GENERATOR FROM R 
C
      external rndstart
      external rndend
      external rndunif

      DOUBLE PRECISION rndunif

      call rndstart()

C
C  PARAMETER CHECK AND INITIALIZATION
C
      NPRCHK=N.GT.0.AND.TOL.GT.0.AND.ISEED.NE.0
      IF (.NOT.NPRCHK) CALL MESSGE(500,'REGTAU',1)
      N1=N/2
      IPS0=IPSI
      XK0=XK
      IPSI=4
      IO=N
      TAI=1.D6
      KSEED=ISEED
      DO 500 I=1,NN
c       call intpr('i',1,i,1)
        DO 30 K=1,2

C  MODIFIED CLAUDIO 16/08/2012
C  TO USE THE UNIFORM RANDOM NUMBER GENERATOR FROM R 
C
CCCCC   10     CALL RANDOW(KSEED,RND)
 10       RND = rndunif()
          IK=INT(RND*FLOAT(N))+1
          IF (IK.GT.N) IK=N
          IF (K.EQ.1) THEN
            IND(1)=IK
            IND1=IK
            GOTO 30
          ENDIF
          IF (IK.EQ.IND1) GOTO 10
          IF (DABS(X(IND1)-X(IK)).LE.1.D-5) GOTO 10
          IND(2)=IK
   30   CONTINUE
      DO 100 K=1,2
      IK=IND(K)
      SX(K)=X(IK)
  100 SY(K)=Y(IK)
      SB(I)=(SY(2)-SY(1))/(SX(2)-SX(1))
      SA(I)=SY(1)-SB(I)*SX(1)
      DO 120 J=1,N
      RS(J)=Y(J)-SB(I)*X(J)-SA(I)
      SRJ=DABS(RS(J))
      TMP1(J)=SNGL(SRJ)
      TMP2(J)=FLOAT(J)
  120 CONTINUE
      CALL SRT2(TMP1,TMP2,N,1,N)
      SUMXX=0.D0
      SUMXY=0.D0
      SUMX=0.D0
      SUMY=0.D0
      DO 140 J=1,N1
      K=TMP2(J)
      SUMXX=SUMXX+X(K)*X(K)
      SUMXY=SUMXY+X(K)*Y(K)
      SUMX=SUMX+X(K)
      SUMY=SUMY+Y(K)
  140 CONTINUE 
      SB(I)=(SUMXY-SUMX*SUMY/DFLOAT(N1))/(SUMXX-SUMX*SUMX/DFLOAT(N1))
      SA(I)=(SUMY-SB(I)*SUMX)/DFLOAT(N1)
      DO 160 J=1,N
      RS(J)=Y(J)-SB(I)*X(J)-SA(I)
      SRJ=DABS(RS(J)*W(J))
      TMP1(J)=SNGL(SRJ)
  160 CONTINUE
      CALL SRT1(TMP1,N,1,N)
      S0=TMP1(N1+1)
      IF (2*N1.EQ.N) S0=(TMP1(N1)+S0)/2.0
      S0=S0/0.6745
      TAU=TOL
      IF (S0.LE.TOL) GOTO 300
      IT=0
      H=1.0
      XK=C1
c      call realpr('s0',2,s0,1)
  180 IT=IT+1
      SUM=0.0
C     RHO4(S)=RHO4(-S) !!
      DO 200 J=1,N
      U=TMP1(J)
  200 SUM=SUM+RHO4(U/S0)
      S1=SUM/(FLOAT(N)*B1)
      S1=S0*SQRT(S1)
      H=ABS(S1-S0)/S0
      S0=S1
      IF (H.GT.TOL.AND.IT.LT.50) GOTO 180
c      call realpr('h',1,h,1)
c      call intpr('it',2,it,1)
  300 CONTINUE

      IF (S0.LE.TOL) GOTO 400
      SUM=0.0
      XK=C2
      DO 320 J=1,N
      U=TMP1(J)
  320 SUM=SUM+RHO4(U/S0)
      S1=SUM/(FLOAT(N)*B2)
      TAU=S0*SQRT(S1)      
  400 TA(I)=DBLE(TAU)
      IF (TA(I).GE.TAI) GOTO 500
      IO=I
      TAI=TA(I)
  500 CONTINUE
      IPSI=IPS0
      XK=XK0
      AO=SA(IO)
      BO=SB(IO)
      TO=TA(IO)

C  ADDED CLAUDIO 16/08/2012
C  TO USE THE UNIFORM RANDOM NUMBER GENERATOR FROM R 
C
      call rndend()

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION RHO4(S)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  GIVES THE VALUE OF THE INTEGRAL FROM 0 TO S OF PSI(T)
C
      COMMON/RHOTAU/IPSI,C,H1,H2,H3,XK,D
      IPS=IABS(IPSI)
      ABST=ABS(S)
      S2=S*S
      IF (IPS.EQ.0) GOTO 100
      IF (IPS.EQ.1) GOTO 200
      IF (IPS.EQ.2) GOTO 300
      IF (IPS.EQ.3) GOTO 400
      IF (IPS.EQ.4) GOTO 500
      IF (IPS.EQ.10) GOTO 700
C      RHO4=URHO(S)
C      RETURN
  100 RHO4=S2/2.
      RETURN
  200 TMP=S2/2.
      IF (ABST.GT.C) TMP=C*(ABST-C/2.)
      GOTO 600
  300 IF (ABST.GT.H2) GOTO 350
      TMP=S2/2.
      IF (ABST.GT.H1) TMP=H1*(ABST-H1/2.)
      GOTO 600
  350 TMP=0.5*H1*(H3+H2-H1)
      IF (ABST.LT.H3) TMP=TMP-.5*H1*(H3-ABST)**2/(H3-H2)
      GOTO 600
  400 TMP=1./6.
      IF (ABST.GE.1.) GOTO 600
      TMP=(S2*(S2-3)+3)*S2/6.
      GOTO 600
  500 TMP=1.
      IF (ABST.GE.XK) GOTO 600
      S2=S2/(XK**2)
      TMP=(S2*(S2-3)+3)*S2
      GOTO 600
  700 TMP=S2/2.
      IF (S.LT.H1) TMP=H1*(S-H1/2.)
      IF (S.GT.H2) TMP=H2*(S-H2/2.)
  600 RHO4=TMP
      RETURN
      END
C
      SUBROUTINE MESSGE(NUMBER,ITEXT,ISTOP)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      CHARACTER *6 ITEXT, CC*34
      IF (ISTOP.EQ.1) THEN
C
C Error exit from R
C
       CC='Input parameter error(s) in '//ITEXT
       CALL REXIT(CC)
      ELSE
       CC='Warning message in '//ITEXT
       CALL INTPR(CC,25,NUMBER,1)
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SRT1(A,N,K1,K2)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
      REAL A(N)
      LOGICAL NPRCHK
C
      NPRCHK=K1.GE.1.AND.K2.GT.K1.AND.K2.LE.N
      IF (.NOT.NPRCHK) CALL MESSGE(500,'SRT1  ',1)
      N1=K2-K1+1
c      I=1
c   10 I=I+I
c      IF (I.LE.N1) GOTO 10
      M=N1
   20 M=M/2
      IF (M.EQ.0) GOTO 90
      K=N1-M
      DO 40 J=1,K
      L=J
   50 IF (L.LT.1) GOTO 40
      LPM=L+M
      LPM1=LPM+K1-1
      L1=L+K1-1
      IF (A(LPM1).GE.A(L1)) GOTO 40
      X=A(LPM1)
      A(LPM1)=A(L1)
      A(L1)=X
      L=L-M
      GOTO 50
   40 CONTINUE
      GOTO 20
   90 CONTINUE
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SRT2(A,B,N,K1,K2)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
      REAL A(N),B(N)
      LOGICAL NPRCHK
C
      NPRCHK=N.GT.0.AND.K1.GE.1.AND.K2.GE.K1.AND.K2.LE.N
      IF (.NOT.NPRCHK) CALL MESSGE(500,'SRT2  ',1)
      N1=K2-K1+1
c      I=1
c   10 I=I+I
c      IF (I.LE.N1) GOTO 10
      M=N1
   20 M=M/2
      IF (M.EQ.0) GOTO 90
      K=N1-M
      DO 40 J=1,K
      L=J
   50 IF (L.LT.1) GOTO 40
      LPM=L+M
      LPM1=LPM+K1-1
      L1=L+K1-1
      IF (A(LPM1).GE.A(L1)) GOTO 40
      X=A(LPM1)
      Y=B(LPM1)
      A(LPM1)=A(L1)
      B(LPM1)=B(L1)
      A(L1)=X
      B(L1)=Y
      L=L-M
      GOTO 50
   40 CONTINUE
      GOTO 20
   90 CONTINUE
      END

