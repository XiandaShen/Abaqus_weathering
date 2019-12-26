C     WEATHERING  

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

C FOR CONSTANT
      PARAMETER(FTOL=1.0D-10,PI=3.1415926D0)
C
      COMMON/PARAMETERS/BETA,ALPHA,BMO,SMO
C
C----------------------------------------------------------------
      DIMENSION SYMIF(3,3,3,3),DTSTRESS(3,3),STIFFM(3,3,3,3),
     1 TSTRAIN(3,3),DTSTRAIN(3,3),TSTRESS(3,3),HCMATRIX(6,6)
C	 
	  DOUBLE PRECISION SMM,BMM,DAMAGE,DAMAGEN,DALAMBDA
C      
C
C   FOR DAMAGE
      DIMENSION OSTRAIN(3,3)
      DATA ZERO,HALF,ONE,TWO,THREE /0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/
C
C---------------------------------------------------------------------------
C 	                      MATERIAL PARAMETERS
C----------------------------------------------------------------------------
C
C STATE VARIABLE: DAMAGE 
C 7 WEAHTERING TIME
	  IF (TIME(2).EQ.0) THEN	
		STATEV(1)=0.0D0
	  END IF
	  DAMAGE=STATEV(1)
	  BMO=PROPS(1)
	  SMO=PROPS(2)
	  BETA=PROPS(3)
	  ALPHA=PROPS(4)	
C
C-----------------------------------------------------------------------
C                          INPUT INITIATION
C-----------------------------------------------------------------------
C    INPUT FOR EACH TIME STEP, FROM XY TO XZY----Y CHANGED TO Z
		CALL MAT1_MAT2(NTENS,STRAN,TSTRAIN,HALF)
		CALL MAT1_MAT2(NTENS,DSTRAN,DTSTRAIN,HALF)
		CALL MAT1_MAT2(NTENS,STRESS,TSTRESS,ONE)
        CALL Aij_PLUS_Bij(TSTRAIN,DTSTRAIN,TSTRAIN)
		CALL SYMIDENDITYF(SYMIF)
C----------------------------------------------------------------
C         CALCULATE DAMAGE IN MATRIX
	 	  DAMAGEN=0.0D0	
	  CALL THERMODYNAMICS(TSTRAIN,ALPHA,BETA,DAMAGE,
     1 DAMAGEN)
	  IF (DAMAGEN.GT.DAMAGE) THEN
		DAMAGE=DAMAGEN
      END IF
	  DALAMBDA=BMO-2.0D0/3.0D0*SMO
	  SMM=SMO+2.0D0*BETA*DAMAGE
	  BMM=2.0D0/3.0D0*SMO+DALAMBDA+2.0D0*ALPHA*DAMAGE
C----------------------------------------------------------------
		CALL STIFFNESS(BMM,SMM,STIFFM)       
C----------------------------------------------------------------	  
C----------------------------------------------------------------
C    STRESS
       CALL Aijkl_Bkl(STIFFM,DTSTRAIN,DTSTRESS)
	   CALL Aij_PLUS_Bij(TSTRESS,DTSTRESS,TSTRESS)
C     OUTPUT
      CALL MAT2_MAT1(NTENS,TSTRESS,STRESS,ONE)	 
      CALL MAT4_MAT2(STIFFM,HCMATRIX,1)        
C      WRITE(6,*) 'stress=',TSTRESS
C      WRITE(6,*) STRESS		
      DO I = 1 , NTENS
        DO J = 1 , NTENS
          DDSDDE (I , J) = HCMATRIX (I , J)
        END DO
      END DO
	  STATEV(1)=DAMAGE  
	  RETURN
	 END
C
C
C
C-----------------------------------------------------------------------
C ====================================================================
C 
C               TRACE  
C
C ====================================================================
C-----------------------------------------------------------------------	 
      SUBROUTINE TRACE(TENSOR,R)
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TENSOR(3,3)
      R=0.0D0
C
      DO I=1,3
          R=R+TENSOR(I,I)
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------	 
C
      SUBROUTINE MAT2_MAT1(NTENS,TENSOR,VECTOR,FACT)
C
C ====================================================================
C
C =================== MAT2_MAT1: TENSOR TO VECTOR=====================
C
C ====================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
C
      DO I=1,NTENS
          VECTOR(I)=0.0D0
      END DO
C
      VECTOR( 1 ) = TENSOR(1 , 1)
      VECTOR( 2 ) = TENSOR(2 , 2)
      VECTOR( 3 ) = TENSOR(3 , 3)
      VECTOR( 4 ) = TENSOR(1 , 2)*FACT
C
      RETURN
      END
C
C	 
C-----------------------------------------------------------------------	 
      SUBROUTINE MAT1_MAT2(NTENS,VECTOR,TENSOR,FACT)
C ====================================================================
C 
C                MAT1_MAT2 : 2D PLAIN STRAIN VECTOR TO 3D TENOR  
C
C ====================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
C
      DO I=1,3
        DO J=1,3
          TENSOR(I,J)=0.0D0
        END DO
      END DO
C
      TENSOR(1 , 1) = VECTOR( 1 )
      TENSOR(2 , 2) = VECTOR( 2 )
	  TENSOR(3 , 3) = VECTOR( 3 )
      TENSOR(1 , 2) = VECTOR( 4 )*FACT
      TENSOR(2 , 1) = VECTOR( 4 )*FACT
C
      RETURN
      END
C-----------------------------------------------------------------------
C ====================================================================
C 
C                STIFFNESS 
C
C ====================================================================
      SUBROUTINE STIFFNESS(BU,AMU,R)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION R(3,3,3,3)
        ALAMDA=BU-2.0D0*AMU/3.0D0
        DO I=1,3
        DO J=1,3
            DO K=1,3
                DO L=1,3
                    R(I,J,K,L)=0.0D0
                END DO
            END DO
        END DO
        END DO
        R(1,1,1,1)=2.0*amu+alamda
        R(2,2,2,2)=2.0*amu+alamda
        R(3,3,3,3)=2.0*amu+alamda
        R(1,1,2,2)=alamda
        R(1,1,3,3)=alamda
        R(2,2,1,1)=alamda
        R(2,2,3,3)=alamda
        R(3,3,1,1)=alamda
        R(3,3,2,2)=alamda
        R(2,3,2,3)=amu
        R(2,3,3,2)=amu
        R(3,2,3,2)=amu
        R(3,2,2,3)=amu
        R(2,1,2,1)=amu
        R(2,1,1,2)=amu
        R(1,2,2,1)=amu
        R(1,2,1,2)=amu
        R(1,3,1,3)=amu
        R(1,3,3,1)=amu
        R(3,1,1,3)=amu
        R(3,1,3,1)=amu
		RETURN
      END
C-----------------------------------------------------------------------
C ====================================================================
C 
C                INVERSEFOURTH 
C
C ====================================================================  
      SUBROUTINE INVERSEFOURTH(F,FINV)
            INCLUDE 'ABA_PARAM.INC' 
      DIMENSION F(3,3,3,3),a(6,6),W(6,6),FINV(3,3,3,3)
	  DOUBLE PRECISION CONS
	  INTEGER IX,IY
	
		a(1,1)=F(1,1,1,1)
		a(1,2)=F(1,1,2,2)
		a(1,3)=F(1,1,3,3)
		a(1,4)=F(1,1,1,2)
		a(1,5)=F(1,1,2,3)
		a(1,6)=F(1,1,3,1)
		a(2,1)=F(2,2,1,1)
		a(2,2)=F(2,2,2,2)
		a(2,3)=F(2,2,3,3)
		a(2,4)=F(2,2,1,2)
		a(2,5)=F(2,2,2,3)
		a(2,6)=F(2,2,3,1)
		a(3,1)=F(3,3,1,1)
		a(3,2)=F(3,3,2,2)
		a(3,3)=F(3,3,3,3)
		a(3,4)=F(3,3,1,2)
		a(3,5)=F(3,3,2,3)
		a(3,6)=F(3,3,3,1)
		a(4,1)=F(1,2,1,1)
		a(4,2)=F(1,2,2,2)
		a(4,3)=F(1,2,3,3)
		a(4,4)=F(1,2,1,2)
		a(4,5)=F(1,2,2,3)
		a(4,6)=F(1,2,3,1)
		a(5,1)=F(2,3,1,1)
		a(5,2)=F(2,3,2,2)
		a(5,3)=F(2,3,3,3)
		a(5,4)=F(2,3,1,2)
		a(5,5)=F(2,3,2,3)
		a(5,6)=F(2,3,3,1)
		a(6,1)=F(3,1,1,1)
		a(6,2)=F(3,1,2,2)
		a(6,3)=F(3,1,3,3)
		a(6,4)=F(3,1,1,2)
		a(6,5)=F(3,1,2,3)
		a(6,6)=F(3,1,3,1)
      CALL INVERSE(a,6,6,W)
C	  
          DO I=1,3
              DO J=1,3
                  DO K=1,3
                      DO L=1,3
                          FINV(I,J,K,L)=0.0D0
                      END DO
                  END DO
              END DO
          END DO
C		  
          DO I=1,3
              DO J=1,3
				DO K=1,3
					DO L=1,3
						CONS=2.0D0
						IF (I.EQ.1.AND.J.EQ.1) THEN
						IX=1
						END IF
						IF (I.EQ.2.AND.J.EQ.2) THEN
						IX=2
						END IF
						IF (I.EQ.3.AND.J.EQ.3) THEN
						IX=3
						END IF
						IF (I.EQ.1.AND.J.EQ.2) THEN
						IX=4
						END IF
						IF (I.EQ.2.AND.J.EQ.1) THEN
						IX=4
						END IF
						IF (I.EQ.3.AND.J.EQ.2) THEN
						IX=5
						END IF
						IF (I.EQ.2.AND.J.EQ.3) THEN
						IX=5
						END IF
						IF (I.EQ.1.AND.J.EQ.3) THEN
						IX=6
						END IF 
						IF (I.EQ.3.AND.J.EQ.1) THEN
						IX=6
						END IF
						IF (K.EQ.1.AND.L.EQ.1) THEN
						IY=1
						END IF
						IF (K.EQ.2.AND.L.EQ.2) THEN
						IY=2
						END IF
						IF (K.EQ.3.AND.L.EQ.3) THEN
						IY=3
						END IF
						IF (K.EQ.1.AND.L.EQ.2) THEN
						IY=4
						END IF
						IF (K.EQ.2.AND.L.EQ.1) THEN
						IY=4
						END IF
						IF (K.EQ.3.AND.L.EQ.2) THEN
						IY=5
						END IF
						IF (K.EQ.2.AND.L.EQ.3) THEN
						IY=5
						END IF
						IF (K.EQ.1.AND.L.EQ.3) THEN
						IY=6
						END IF 
						IF (K.EQ.3.AND.L.EQ.1) THEN
						IY=6
						END IF
						IF (IX.LE.3.AND.IY.LE.3) THEN
						CONS=1.0D0
						END IF
						IF (IX.GT.3.AND.IY.GT.3) THEN
						CONS=4.0D0
						END IF    						
                        FINV(I,J,K,L)=W(IX,IY)/CONS
					END DO
				END DO
              END DO
          END DO
        return
      END
C-----------------------------------------------------------------------	  
C ====================================================================
C 
C                ESHELBY 
C
C ====================================================================  
      SUBROUTINE ESHELBY(ANU,R)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION R(3,3,3,3)
      DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
                      DIJ=0.0D0
                      DKL=0.0D0
                      DIK=0.0D0
                      DJL=0.0D0
                      DIL=0.0D0
                      DJK=0.0D0
              IF(I.EQ.J) THEN
              DIJ=1.0D0
              END IF
              IF(K.EQ.L) THEN
              DKL=1.0D0
              END IF              
              IF(I.EQ.K) THEN
              DIK=1.0D0
              END IF              
              IF(J.EQ.L) THEN
              DJL=1.0D0
              END IF
              IF(I.EQ.L) THEN
              DIL=1.0D0
              END IF
              IF(J.EQ.K) THEN
              DJK=1.0D0
              END IF
              R(I,J,K,L)=(5.0D0*ANU-1.0D0)/15.0D0/(1.0D0-ANU)*DIJ*DKL
     1+(-5.0D0*ANU+4.0D0)/15.0D0/(1.0D0-ANU)*(DIK*DJL+DIL*DJK)
                  END DO
              END DO 
          END DO
      END DO
	  RETURN
      END
C-----------------------------------------------------------------------	  
        SUBROUTINE TRANSFORMATION(PSI,THETA,PHI,QQ)
CC====================================================================================================
C
C                          TRANSFORMATION MATRIX Q	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION QQ(3,3)
      DOUBLE PRECISION PSI,THETA,PHI
C
		DO I=1,3
		DO J=1,3
		QQ(I,J)=0.0D0
		END DO
		END DO
      QQ(1,1)=COS(PSI)*COS(THETA)*COS(PHI)
     1-SIN(PSI)*SIN(PHI)
      QQ(1,2)=SIN(PSI)*COS(THETA)*COS(PHI)
     1+COS(PSI)*SIN(PHI)
      QQ(1,3)=-SIN(THETA)*COS(PHI)
      QQ(2,1)=-COS(PSI)*COS(THETA)*SIN(PHI)
     1-SIN(PSI)*COS(PHI)
      QQ(2,2)=-SIN(PSI)*COS(THETA)*SIN(PHI)
     1+COS(PSI)*COS(PHI)
      QQ(2,3)=SIN(THETA)*SIN(PHI)
      QQ(3,1)=COS(PSI)*SIN(THETA)
      QQ(3,2)=SIN(PSI)*SIN(THETA)
      QQ(3,3)=COS(THETA)
C
      RETURN
      END
C-----------------------------------------------------------------------
CC====================================================================================================
C
C                          SYMIDENDITYF	  
C
C=====================================================================================================
C
      SUBROUTINE  SYMIDENDITYF(R)
            INCLUDE 'ABA_PARAM.INC'
         DIMENSION R(3,3,3,3)
         R(:,:,:,:)=0.0D0
        DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
        R(1,1,1,1)=1.0D0
        R(2,2,2,2)=1.0D0
        R(3,3,3,3)=1.0D0
        R(2,3,2,3)=0.5D0
        R(2,3,3,2)=0.5D0
        R(3,2,3,2)=0.5D0
        R(3,2,2,3)=0.5D0
        R(2,1,2,1)=0.5D0
        R(2,1,1,2)=0.5D0
        R(1,2,2,1)=0.5D0
        R(1,2,1,2)=0.5D0
        R(1,3,1,3)=0.5D0
        R(1,3,3,1)=0.5D0
        R(3,1,1,3)=0.5D0
        R(3,1,3,1)=0.5D0
                  END DO
              END DO
          END DO
        END DO
      RETURN
      end 
C-----------------------------------------------------------------------
CC====================================================================================================
C
C                          TRANSPOSEF	  
C
C=====================================================================================================
C	
      SUBROUTINE TRANSPOSEF(A,B)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(3,3,3,3),B(3,3,3,3)
      DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
                      B(I,J,K,L)=A(K,L,I,J)
                  END DO
              END DO
          END DO
      END DO
      RETURN
      END 	
C-----------------------------------------------------------------------
C ====================================================================
C                        MAT4_MAT2                                                                  I
C I        THIS PROGRAM TRANSFORMS THE FOURTH ORDER COMPLIANCE       I
C I        TENSOR TO A SECOND ORDER MATRIX                           I
C I                                                                  I
C ====================================================================
C
      SUBROUTINE MAT4_MAT2(TENSOR,DMATRIX,ICOE)
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TENSOR(3,3,3,3),DMATRIX(6,6)
C
      DATA ZERO,TWO /0.0D0,2.0D0/
C
C     D2 = THE SECOND ORDER STIFFNESS MATRIX
C
      DO I=1,6
        DO J=1,6
          DMATRIX(I,J)=0.0D0
        END DO
      END DO

      IF (ICOE.EQ.1) THEN
         COE1=1.
         COE2=1.
      ELSEIF(ICOE.EQ.2) THEN
         COE1=2.
         COE2=4.
      END IF
C
      DMATRIX(1,1)=TENSOR(1,1,1,1)
      DMATRIX(1,2)=TENSOR(1,1,2,2)
      DMATRIX(1,3)=TENSOR(1,1,3,3)
      DMATRIX(1,4)=TENSOR(1,1,1,2)*COE1
      DMATRIX(1,5)=TENSOR(1,1,2,3)*COE1
      DMATRIX(1,6)=TENSOR(1,1,1,3)*COE1
C
      DMATRIX(2,1)=TENSOR(2,2,1,1)
      DMATRIX(2,2)=TENSOR(2,2,2,2)
      DMATRIX(2,3)=TENSOR(2,2,3,3)
      DMATRIX(2,4)=TENSOR(2,2,1,2)*COE1
      DMATRIX(2,5)=TENSOR(2,2,2,3)*COE1
      DMATRIX(2,6)=TENSOR(2,2,1,3)*COE1
C
      DMATRIX(3,1)=TENSOR(3,3,1,1)
      DMATRIX(3,2)=TENSOR(3,3,2,2)
      DMATRIX(3,3)=TENSOR(3,3,3,3)
      DMATRIX(3,4)=TENSOR(3,3,1,2)*COE1
      DMATRIX(3,5)=TENSOR(3,3,2,3)*COE1
      DMATRIX(3,6)=TENSOR(3,3,1,3)*COE1
C
      DMATRIX(4,1)=TENSOR(1,2,1,1)*COE1
      DMATRIX(4,2)=TENSOR(1,2,2,2)*COE1
      DMATRIX(4,3)=TENSOR(1,2,3,3)*COE1
      DMATRIX(4,4)=TENSOR(1,2,1,2)*COE2
      DMATRIX(4,5)=TENSOR(1,2,2,3)*COE2
      DMATRIX(4,6)=TENSOR(1,2,1,3)*COE2
C
      DMATRIX(5,1)=TENSOR(2,3,1,1)*COE1
      DMATRIX(5,2)=TENSOR(2,3,2,2)*COE1
      DMATRIX(5,3)=TENSOR(2,3,3,3)*COE1
      DMATRIX(5,4)=TENSOR(2,3,1,2)*COE2
      DMATRIX(5,5)=TENSOR(2,3,2,3)*COE2
      DMATRIX(5,6)=TENSOR(2,3,1,3)*COE2
C
      DMATRIX(6,1)=TENSOR(1,3,1,1)*COE1
      DMATRIX(6,2)=TENSOR(1,3,2,2)*COE1
      DMATRIX(6,3)=TENSOR(1,3,3,3)*COE1
      DMATRIX(6,4)=TENSOR(1,3,1,2)*COE2
      DMATRIX(6,5)=TENSOR(1,3,2,3)*COE2
      DMATRIX(6,6)=TENSOR(1,3,1,3)*COE2
C
C
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------	  
      SUBROUTINE A_Bij(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION B(3,3),C(3,3)
      DO I=1,3
          DO J=1,3
          C(I,J)=A*B(I,J)
          END DO
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE A_Bijkl(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION B(N,N,N,N),C(N,N,N,N)
      DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
          C(I,J,K,L)=A*B(I,J,K,L)
                  END DO
              END DO
          END DO
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Ai_Bi(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION A(N),B(N)
      C=0.0D0
      DO I=1,3
          C=C+A(I)*B(I)
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Ai_Bj(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION A(N),B(N),C(N,N)
      C(:,:)=0.0D0
      DO I=1,3
          DO J=1,3
          C(I,J)=A(I)*B(J)
          END DO
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Ai_PLUS_Bi(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION A(N),B(N),C(N)
      DO I=1,3
          C(I)=A(I)+B(I)
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aij_Bj(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION A(N,N)
      DIMENSION B(N),C(N)
      DO I=1,3
          C(I)=0.0D0
          DO J=1,3
          C(I)=C(I)+A(I,J)*B(J)
          END DO
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aij_MINUS_BijN(A,B,N,C)
C====================================================================================
C                                                                      *
C                                  Aij_MINUS_BijN                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
C      INTEGER N
      DIMENSION A(N,N),B(N,N),C(N,N)
C
      DO I=1,3
        DO J=1,3
          C(I,J)=A(I,J)-B(I,J)
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aij_PLUS_Bij(A,B,C)
C====================================================================================
C                                                                      *
C                                  Aij_PLUS_Bij                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3),B(3,3),C(3,3)
C
      DO I=1,3
        DO J=1,3
          C(I,J)=A(I,J)+B(I,J)
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aij_PLUS_BijN(A,B,N,C)
C====================================================================================
C                                                                      *
C                                  Aij_PLUS_BijN                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
C      INTEGER N
      DIMENSION A(N,N),B(N,N),C(N,N)
C
      DO I=1,3
        DO J=1,3
          C(I,J)=A(I,J)+B(I,J)
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aijkl_Bkl(A,B,C)
C========================================================================
C                                                                       =
C                              Aijkl_Bkl                                =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3),C(3,3)
      DATA ZERO /0.0D0/
C
      DO I=1,3
        DO J=1,3
          C(I,J)=ZERO
          DO K=1,3
            DO L=1,3
              C(I,J)=C(I,J)+A(I,J,K,L)*B(K,L)
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aijkl_MINUS_Bijkl(A,B,C)
C====================================================================================
C                                                                      *
C                                  Aij_MINUS_Bij                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
C
      DO I=1,3
        DO J=1,3
            DO K=1,3
                DO L=1,3
          C(I,J,K,L)=A(I,J,K,L)-B(I,J,K,L)
                END DO
            END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aijkl_PLUS_Bijkl(A,B,C)
C====================================================================================
C                                                                      *
C                                  Aij_PLUS_Bij                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
C
      DO I=1,3
        DO J=1,3
            DO K=1,3
                DO L=1,3
          C(I,J,K,L)=A(I,J,K,L)+B(I,J,K,L)
                END DO
            END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aijmn_Bmnkl(A,B,C)
CC====================================================================================================
C
C                          Aijmn_Bmnkl=Cijkl	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
C
C
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              C(I,J,K,L)=0.0D0
              DO M=1,3
                DO N=1,3
            C(I,J,K,L)=C(I,J,K,L)+A(I,J,M,N)*B(M,N,K,L)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
C
		RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aik_Bkj(A,B,C)
CC====================================================================================================
C
C                          Aik_Bkj=Cij	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3),B(3,3),C(3,3)
C
C
      DO I=1,3
        DO J=1,3
          C(I,J)=0.0D0
          DO K=1,3
            C(I,J)=C(I,J)+A(I,K)*B(K,J)
          END DO
        END DO
      END DO
C
		RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aik_BkjN(A,B,N,C)
CC====================================================================================================
C
C                          Aik_Bkj=Cij	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      INTEGER N
      DIMENSION A(N,N),B(N,N),C(N,N)
C
C
      DO I=1,3
        DO J=1,3
          C(I,J)=0.0D0
          DO K=1,3
            C(I,J)=C(I,J)+A(I,K)*B(K,J)
          END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE Aij_Bijkl(A,B,C)
C========================================================================
C                                                                       =
C                              Aij_Bijkl                                =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DOUBLE PRECISION A(3,3),B(3,3,3,3),C(3,3)
      DATA ZERO /0.0D0/
C
      DO K=1,3
        DO L=1,3
          C(K,L)=0.0D0
          DO I=1,3
            DO J=1,3
              C(K,L)=C(K,L)+A(I,J)*B(I,J,K,L)
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE INVERSE(A,N,NP,AINV)
C========================================================================
C
C    CALCULATE THE SECOND ORDER TENSOR A'S INVERSE, AINV
C    A^{-1} = AINV    
C    this subroutine inverses a (n x n) A matrix
C	 following a Gauss-Jordan elimination process
C
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C  
      DIMENSION A(NP,NP),IPIV(NP),INDXR(NP),INDXC(NP),
     1 A0(NP,NP),AINV(NP,NP)
C
      DO J=1,N
        IPIV(J)=0
      END DO
C
C
C     storage of the original A matrix
C
      DO I=1,N
        DO J=1,N
          A0(I,J)=A(I,J)
        END DO
      END DO
C
C	find a pivot among the rows of A that have not already been reduced
C
      DO I=1,N
        BIG=0.0D0
        DO J=1,N
          IF(IPIV(J).NE.1)THEN
            DO K=1,N
                IF(IPIV(K).EQ.0)THEN
                  IF(ABS(A(J,K)).GE.BIG)THEN
                   BIG=ABS(A(J,K))
                   IROW=J
                   ICOL=K
                   PIV=A(J,K)
                  END IF
                ELSEIF(IPIV(K).GT.1)THEN
                  write (7,*) 'Singular Matrix'
              END IF
            END DO
          END IF
        END DO
C
      IPIV(ICOL)=IPIV(ICOL)+1
      INDXR(I)=IROW
      INDXC(I)=ICOL
C	  
C     interchange the rows to put the pivot on the diagonal
C
      IF(irow.ne.icol)THEN
        DO L=1,N
          DUM=A(IROW,L)
          A(IROW,L)=A(ICOL,L)
          A(ICOL,L)=DUM
        END DO
      END IF
C     reduction of the row of the pivot
C       
      IF(PIV.EQ.0) write (7,*) 'Singular Matrix2'
C       
      PIVIN=1./PIV          ! numerical stabilization
C
      A(ICOL,ICOL)=1.       ! numerical stabilization
        DO L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVIN
        END DO
C
C     reduction of the column of the pivot
C
        DO LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.    ! numerical stabilization
            DO L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
            END DO
          END IF
        END DO
      END DO
C
C
C     unscramble the columns to get A-1
C		
      DO J=N,1,-1   ! reverse DO loop
        DO K=1,N
          DUM=A(K,INDXR(J))
          A(K,INDXR(J))=A(K,INDXC(J))
          A(K,INDXC(J))=DUM
        END DO
      END DO
C
C	restitution process of A and Ainv
C
      DO I=1,N
        DO J=1,N
          AINV(I,J)=A(I,J)
        END DO
      END DO
C
      DO I=1,N
        DO J=1,N
          A(I,J)=A0(I,J)
        END DO
      END DO
C     
      RETURN
      END
C----------------------------------------------------------------
C----------------------------------------------------------------
C      UNIFORM DISTRIBUTED
C----------------------------------------------------------------
C----------------------------------------------------------------
      SUBROUTINE r8_uniform_01 ( Iseed,RANUNI )
       INCLUDE 'ABA_PARAM.INC'

C      integer i4_huge 
      parameter ( i4_huge = 2147483647 )
C      integer k
      double precision RANUNI
C      integer seed

      k = Iseed / 127773

      Iseed = 16807 * ( Iseed - k * 127773 ) - k * 2836

      if ( Iseed .lt. 0 ) then
        Iseed = Iseed + i4_huge
      end if

      RANUNI = dble ( Iseed ) * 4.656612875D-10

      return
      end
C----------------------------------------------------------------
C----------------------------------------------------------------
C      LOGNORMAL DISTRIBUTED
C----------------------------------------------------------------
C----------------------------------------------------------------	
      SUBROUTINE rlognorm ( alm,ald,iseed,RANLOG )
             INCLUDE 'ABA_PARAM.INC'
      integer iseed
      a=log((alm**2)/(((alm**2)+ald)**(0.5D0)))
      b=(log(ald/(alm**2)+1.0D0))**(0.5D0)
	  CALL r8_normal_ab ( a, b, iseed,R1 )
      RANLOG=exp(R1)
	  return
      end 
C----------------------------------------------------------------
C----------------------------------------------------------------
C      NORMAL DISTRIBUTED
C----------------------------------------------------------------
C----------------------------------------------------------------	
      SUBROUTINE r8_normal_ab ( a, b, iseed ,RANNOR)
       INCLUDE 'ABA_PARAM.INC'
      parameter ( r8_pi = 3.141592653589793D00 )
      integer iseed
      CALL r8_uniform_01 ( iseed ,R1)
      CALL r8_uniform_01 ( iseed,R2 )
      x = sqrt ( -2.0D0 * log ( R1 ) ) 
     1  * cos ( 2.0D0 * r8_pi * R2 )

      RANNOR = a + b * x

      return
      end	 
C----------------------------------------------------------------
C----------------------------------------------------------------
C----------------------------------------------------------------
C
C                          P TENSOR
C
C----------------------------------------------------------------
C----------------------------------------------------------------
	  SUBROUTINE PTENSOR(A,C,BM,SM,R)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION R(3,3,3,3)
	  DOUBLE PRECISION A,C,BM,SM,W,H,PK,PL,PN,PM,PP
	  W=C/A
	  H=W*(ACOS(W)-W*SQRT(1.0D0-W**2))/(1.0D0-W**2)**(1.5D0)
	  PK=((7.0D0*H-2.0D0*W**2-4.0D0*H*W**2)*SM+3.0D0*(H-2.D0
     1 *W**2+2.0D0*H*W**2)*BM)/8.0D0/(1.0D0-W**2)/SM/(
     2 4.0D0*SM+3.0D0*BM)
	  PL=(SM+3.0D0*BM)*(2.0D0*W**2-H-2.0D0*H*W**2)/4.0D0
     1 /(1.0D0-W**2)/SM/(4.0D0*SM+3.0D0*BM)
	  PN=((6.0D0-5.0D0*H-8.0D0*W**2+8*H*W**2)*SM+3.0D0*
     1(H-2.0D0*W**2+2.0D0*H*W**2)*BM)/2.0D0/(1.0D0-W**2)/
     2 SM/(4.0D0*SM+3.0D0*BM)
	  PM=((15.0D0*H-2.0D0*W**2-12.0D0*H*W**2)*SM+3.0D0*
     1 (3.0D0*H-2.0D0*W**2)*BM)/16.0D0/(1.0D0-W**2)/SM/
     2 (4.0D0*SM+3.0D0*BM)
	  PP=(2.0D0*(4.0D0-3.0D0*H-2.0D0*W**2)*SM+3.0D0*(2.0D0-
     1 3.0D0*H+2.0D0*W**2-3.0D0*H*W**2)*BM)/8.0D0/(1.0D0-W**2)
     2 /SM/(4.0D0*SM+3.0D0*BM)
      DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
					R(I,J,K,L)=0.0D0
                  END DO
              END DO 
          END DO
      END DO
	  R(1,1,1,1)=PK+PM
	  R(2,2,2,2)=PK+PM
	  R(3,3,3,3)=PN
	  R(1,1,2,2)=PK-PM
	  R(1,1,3,3)=PL
	  R(2,2,1,1)=PK-PM
	  R(2,2,3,3)=PL	
	  R(3,3,1,1)=PL	
	  R(3,3,2,2)=PL	
	  R(2,3,2,3)=PP
	  R(2,3,3,2)=PP
	  R(3,2,2,3)=PP
	  R(3,2,3,2)=PP
	  R(3,1,1,3)=PP
	  R(3,1,3,1)=PP
	  R(1,3,1,3)=PP
	  R(1,3,3,1)=PP
	  R(1,2,1,2)=PM
	  R(1,2,2,1)=PM
	  R(2,1,1,2)=PM
	  R(2,1,2,1)=PM	  	  
	  RETURN
      END
C----------------------------------------------------------------
C----------------------------------------------------------------
C----------------------------------------------------------------
C
C                          TRANSP
C
C----------------------------------------------------------------
C----------------------------------------------------------------
      SUBROUTINE TRANSP(TQ,P,PO)
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TQ(3,3),P(3,3,3,3),PO(3,3,3,3)
C
	  PO(:,:,:,:)=0.0D0
	  DO I=1,3
	  DO J=1,3
	  DO K=1,3
	  DO L=1,3
	  DO II=1,3
	  DO JJ=1,3 
	  DO KK=1,3
	  DO LL=1,3
		PO(I,J,K,L)=PO(I,J,K,L)+TQ(I,II)*TQ(J,JJ)*TQ(K,KK)
     1 *TQ(L,LL)*P(II,JJ,KK,LL)
	  END DO
	  END DO
	  END DO
	  END DO
	  END DO
	  END DO
	  END DO 
	  END DO
	  
C
      RETURN
      END	
C----------------------------------------------------------------
C----------------------------------------------------------------
C 
C              THERMODYNAMICS
C 
C----------------------------------------------------------------
C----------------------------------------------------------------
      SUBROUTINE THERMODYNAMICS(STRAINM,A,B,DAMAGE,
     1 DAMAGEN)	  
	  INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION STRAINM(3,3)
      DOUBLE PRECISION DAMAGE,TREPS,YD,FD,DAMAGEN,C0,C1
C
      DIMENSION STEMP21(3,3),STEMP22(3,3),STEMP23(3,3),
     1 STEMP24(3,3)
C      DRIVING FORCE
	  DAMAGEN=DAMAGE
	  CALL TRACE(STRAINM,TREPS)
      CALL A_Bij((-A*TREPS),STRAINM,3,STEMP21)
      CALL A_Bij((-B*2.0D0),STRAINM,3,STEMP22)
      CALL Aik_Bkj(STEMP22,STRAINM,STEMP23)
      CALL Aij_PLUS_Bij(STEMP21,STEMP23,STEMP24)
      CALL TRACE(STEMP24,YD)
C     DAMAGE CRITERION
	  IF (STRAINM(2,2).GT.0.0D0) THEN
	  C0=0.000001
	  C1=0.0005
	  ELSE 
	  C0=0.11
	  C1=2.2
	  END IF
      FD=1.0D0/(2.0D0**0.5D0)*ABS(YD)-C0-C1*DAMAGE
      IF (FD.GT.0.0D0) THEN
		DAMAGEN=(ABS(YD)/(2.0D0**0.5D0)-C0)/C1
      END IF
	  IF (DAMAGEN.GT.0.8D0) THEN
	  DAMAGEN=0.8D0
	  END IF
	  RETURN
	  END
	  