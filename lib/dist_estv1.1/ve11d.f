* COPYRIGHT (c) 1989 AEA Technology.
*######DATE 29 Sept 1992
*######ALIAS VE11A
C######29/9/92. Modified to avoid branches into IF constructs.
      SUBROUTINE VE11AD (N,M,MEQ,A,IA,B,XL,XU,X,ACC,IACT,NACT,PAR,
     1  IPRINT,INFO,W,FGCALC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),PAR(*),W(*)
      EXTERNAL VE11TD,FGCALC
C
C     purpose.
C     This is the entry point to a package of subroutines that
C     calculate the least value of a differentiable function of
C     several variables subject to linear constraints on the
C     values of the variables.
C
C     argument list.
C     N      is the number of variables and must be set by the user.
C     M      is the number of linear constraints (excluding simple
C            bounds) and must be set by the user.
C     MEQ    is the number of constraints that are equalities and
C            must be set by the user.
C     A(.,.) is a 2-dimensional array whose columns are the gradients
C            of the M constraint functions.  Its entries must be set
C            by the user and its dimensions must be at least N by M.
C     IA     is the actual first dimension of the array A that is
C            supplied by the user, so its value may not be less than N.
C     B(.)   is a vector of constraint right hand sides that must
C            also be set by the user.  Specifically the constraints
C            on the variables X(I) I=1(1)N are
C            A(1,K)*X(1)+...+A(N,K)*X(N) .EQ. B(K)  K=1,...,MEQ
C            A(1,K)*X(1)+...+A(N,K)*X(N) .LE. B(K)  K=MEQ+1,...,M  .
C            Note that the data that define the equality constraints
C            come before the data of the inequalities.
C     XL(.)  are vectors whose components must be set to lower and
C     XU(.)  upper bounds on the variables.  Choose very large
C            negative and positive entries if a component should be
C            unconstrained, or set XL(I)=XU(I) to freeze the I-th
C            variable. Specifically these simple bounds are
C            XL(I) .LE. X(I) and X(I) .LE. XU(I)  I=1,...,N  .
C     X(.)   is the vector of variables of the optimization calculation
C            Its initial elements must be set by the user to an
C            estimate of the required solution. The subroutines can
C            usually cope with poor estimates, and there is no need
C            for X(.) to be feasible initially. These variables are
C            adjusted automatically and the values that give the least
C            feasible calculated value of the objective functions are
C            available in X(.) on the return from VE11AD.
C     ACC    is a tolerance on the first order conditions at the
C            calculated solution of the optimization problem. These
C            first order conditions state that, if X(.) is a solution,
C            then there is a set of active constraints with indices
C            IACT(K) K=1(1)NACT, say, such that X(.) is on the
C            boundaries of these constraints, and the gradient of the
C            objective function can be expressed in the form
C            GRAD(F)=PAR(1)*GRAD(C(IACT(1)))+...
C                        ...+PAR(NACT)*GRAD(C(IACT(NACT)))  .
C            Here PAR(K) K=1(1)NACT are Lagrange multipliers that are
C            nonpositive for inequality constraints, and
C            GRAD(C(IACT(K))) is the gradient of the IACT(K)-th
C            constraint function, so it is A(.,IACT(K)) if IACT(K)
C            .LE. M, and it is minus or plus the J-th coordinate vector
C            if the constraint is the lower or upper bound on X(J)
C            respectively.  The normal return from th calculation
C            occurs when X(.) is feasible and the sum of squares of
C            components of the vector RESKT(.) is at most ACC**2,
C            where RESKT(.) is the N-component vector of residuals of
C            the first order condition that is displayed above.
C            Sometimes the package cannot satisfy this condition,
C            because noise in the function values can prevent a change
C            to the variables, no line search being allowed to increas
C            the objective function.
C     IACT(.) is a working space array of integer variables that must
C            be provided by the user.  Its length must be at least
C            (M+2*N).  Its leading entries on the return from the
C            subroutine are the indices IACT(k) K=1(1)NACT that are
C            mentioned in the previous paragraph: in other words they
C            are the indices of the final active constraints.
C            Here the indices M+1,...,M+N and M+N+1,...,M+2*N denote
C            the lower and upper bounds respectively.
C     NACT   is set automatically to the integer variable of this ilk
C            that has been introduced already.
C     PAR    is a one-dimensional array that will hold the Lagrange
C            multipliers PAR(K) K=1(1)NACT on the return from VE11AD,
C            these parameters being defined in the above paragraph
C            on ACC. The length of PAR should be at least N.
C     IPRINT must be set by the user to specify the frequency of
C            printing during the execution of the optimization package.
C            There is no printed output if IPRINT=0. Otherwise, after
C            ensuring feasibility, information is given every
C            IABS(IPRINT) iterations and whenever a parameter called
C            TOL is reduced.  The printing provides the values of X(.),
C            F(.) and G(.)=GRAD(F) if IPRINT is positive, while if
C            IPRINT is negative this information is augmented by the
C            current values of IACT(K) K=1(1)NACT, PAR(K) K=1(1)NACT
C            and RESKT(I) I=1(1)N.  The reason for returning to the
C            calling program is also displayed when IPRINT is nonzero.
C     INFO   is an integer variable that should be set to zero
C            initially, unless the user wishes to impose an upper
C            bound on the number of evaluations of the objective
C            function and its gradient, in which case INFO should be
C            set to the value of this bound.  On the exit from VE11AD
C            it will have one of the following integer values to
C            indicate the reason for leaving the optimization package:
C            INFO=1   X(.) is feasible and the condition that depends
C                     on ACC is satisfied.
C            INFO=2   X(.) is feasible and rounding errors are
C                     preventing further progress.
C            INFO=3   X(.) is feasible but the objective function fails
C                     to decrease although a decrease is predicted by
C                     the current gradient vector. If this return
C                     occurs and KTRES(.) has large components then the
C                     user's calculation of the gradient of the
C                     objective function may be incorrect.  One should
C                     also question the coding of the gradient when
C                     the final rate of convergence is slow.
C            INFO=4   In this case the calculation cannot begin because
C                     IA is less than N or because the lower bound on a
C                     variable is greater than the upper bound.
C            INFO=5   This value indicates that the equality
C                     constraints are inconsistent.   These constraints
C                     include any components of X(.) that are frozen by
C                     setting XL(I)=XU(I).
C            INFO=6   In this case there is an error return because the
C                     equality constraints and the bounds on the
C                     variables are found to be inconsistent.
C            INFO=7   This value indicates that there is no vector of
C                     variables that satisfies all of the constraints.
C                     Specifically, when this return or an INFO=6
C                     return occurs, the current active constraints
C                     (whose indices are IACT(K) K=1(1)NACT) prevent
C                     any change in X(.) that reduces the sum of
C                     constraint violations, where only bounds are
C                     included in this sum if INFO=6.
C            INFO=8   In this case the limit on the number of calls of
C                     subroutine FGCALC (see below) has been reached,
C                     and there would have been further calculation
C                     otherwise.
C     W(.)   is a working space array of real variables that must be
C            provided by the user. Its length must be at least
C            (M+11*N+N**2).  On exit from the package one can find the
C            final components of GRAD(F) and RESKT(.) in W(1),...,W(N)
C            and W(N+1),...,W(2*N) respectively.
C     FGCALC is a user supplied subroutine to calculate the objective
C            function and its gradient.Its first two lines being
C            SUBROUTINE FGCALC (N,X,F,G)
C            DIMENSION X(*),G(*)   .
C            It is called automatically with N set as above and with
C            X(.) set to a feasible vector of variables.  It should
C            calculate the value of the objective function and its
C            gradient for this X(.) and should set them in F and G(I)
C            I=1(1)N respectively, without disturbing N or any
C            of the components of X(.).
C
C Note 1.    The variables N, M, MEQ, IA, ACC and IPRINT and the
C            elements of the arrays A(,.,), B(.), XL(.) and XU(.) are
C            not altered by the optimization procedure.  Their values,
C            the value of INFO and the initial components of X(.)
C            must be set on entry to VE11AD.
C Note 2.    A paper on the method of calculation and a report on the
C            main features of the computer code are available from
C            the author M.J.D.Powell (D.A.M.T.P., University of
C            Cambridge, Silver Street, Cambridge CB3 9EW, England).
C
C     Partition the workspace array.
C
      IG=1
      IRESKT=IG+N
      IZ=IRESKT+N
      IU=IZ+N*N
      IXBIG=IU+N
      IBRES=IXBIG+N
      ID=IBRES+M+N+N
      IZTG=ID+N
      IGM=IZTG+N
      IXS=IGM+N
      IGS=IXS+N
C
C     Call the optimization package.
C
      CALL VE11BD (N,M,MEQ,A,IA,B,XL,XU,X,ACC,IACT,NACT,PAR,IPRINT,
     1  INFO,W(IG),W(IZ),W(IU),W(IXBIG),W(IRESKT),W(IBRES),W(ID),
     2  W(IZTG),W(IGM),W(IXS),W(IGS),FGCALC)
      RETURN
      END
      SUBROUTINE VE11CD (N,M,XL,XU,X,IACT,MEQL,INFO,Z,U,XBIG,RELACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XL(*),XU(*),X(*),IACT(*),Z(*),U(*),XBIG(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
C
C     Set RELACC.
C
      ZTPAR=100.0
      RELACC=ONE
   10 RELACC=0.5*RELACC
      TEMPA=ZTPAR+0.5*RELACC
      TEMPB=ZTPAR+RELACC
      IF (ZTPAR .LT. TEMPA .AND. TEMPA .LT. TEMPB) GOTO 10
C
C     Seek bound inconsistencies and bound equality constraints.
C
      MEQL=0
      DO 20 I=1,N
      IF (XL(I) .GT. XU(I)) GOTO 50
      IF (XL(I) .EQ. XU(I)) MEQL=MEQL+1
   20 CONTINUE
C
C     Initialize U, Z and XBIG.
C
      JACT=0
      NN=N*N
      DO 30 I=1,NN
      Z(I)=ZERO
   30 CONTINUE
      IZ=0
      DO 40 I=1,N
      IF (XL(I) .EQ. XU(I)) THEN
          X(I)=XU(I)
          JACT=JACT+1
          U(JACT)=ONE
          IACT(JACT)=I+M+N
          J=JACT
      ELSE
          J=I+MEQL-JACT
      END IF
      Z(IZ+J)=ONE
      IZ=IZ+N
      XBIG(I)=DABS(X(I))
   40 CONTINUE
      INFO=1
   50 RETURN
      END
      SUBROUTINE VE11HD (N,M,A,IA,IACT,NACT,Z,U,RELACC,INDXBD,ZTC,
     1  CGRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),Z(*),U(*),ZTC(*),CGRAD(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
      NP=NACT+1
      ICON=IACT(INDXBD)
      IACT(INDXBD)=IACT(NP)
      IACT(NP)=ICON
C
C     Form ZTC when the new constraint is a bound.
C
      IF (ICON .GT. M) THEN
          INEWBD=ICON-M
          IF (INEWBD .LE. N) THEN
              TEMP=-ONE
          ELSE
              INEWBD=INEWBD-N
              TEMP=ONE
          END IF
          IZNBD=INEWBD*N-N
          DO 10 J=1,N
          ZTC(J)=TEMP*Z(IZNBD+J)
   10     CONTINUE
C
C     Else form ZTC for an ordinary constraint.
C
      ELSE
          DO 20 I=1,N
          CGRAD(I)=A(I,ICON)
   20     CONTINUE
          DO 35 J=1,N
          ZTC(J)=ZERO
          IZ=J
          DO 30 I=1,N
          ZTC(J)=ZTC(J)+Z(IZ)*CGRAD(I)
          IZ=IZ+N
   30     CONTINUE
   35     CONTINUE
      END IF
C
C     Find any Givens rotations to apply to the last columns of Z.
C
      J=N
   40 JP=J
      J=J-1
      IF (J .GT. NACT) THEN
          IF (ZTC(JP) .EQ. ZERO) GOTO 40
          IF (DABS(ZTC(JP)) .LE. RELACC*DABS(ZTC(J))) THEN
              TEMP=DABS(ZTC(J))
          ELSE IF (DABS(ZTC(J)) .LE. RELACC*DABS(ZTC(JP))) THEN
              TEMP=DABS(ZTC(JP))
          ELSE
              TEMP=DABS(ZTC(JP))*DSQRT(ONE+(ZTC(J)/ZTC(JP))**2)
          END IF
          WCOS=ZTC(J)/TEMP
          WSIN=ZTC(JP)/TEMP
          ZTC(J)=TEMP
C
C     Apply the rotation when the new constraint is a bound.
C
          IZ=J
          IF (ICON .GT. M) THEN
              DO 50 I=1,N
              TEMP=WCOS*Z(IZ+1)-WSIN*Z(IZ)
              Z(IZ)=WCOS*Z(IZ)+WSIN*Z(IZ+1)
              Z(IZ+1)=TEMP
              IZ=IZ+N
   50         CONTINUE
              Z(IZNBD+JP)=ZERO
C
C     Else apply the rotation for an ordinary constraint.
C
          ELSE
              WPIV=ZERO
              DO 60 I=1,N
              TEMPA=WCOS*Z(IZ+1)
              TEMPB=WSIN*Z(IZ)
              TEMP=DABS(CGRAD(I))*(DABS(TEMPA)+DABS(TEMPB))
              IF (TEMP .GT. WPIV) THEN
                  WPIV=TEMP
                  IPIV=I
              END IF
              Z(IZ)=WCOS*Z(IZ)+WSIN*Z(IZ+1)
              Z(IZ+1)=TEMPA-TEMPB
              IZ=IZ+N
   60         CONTINUE
C
C     Ensure orthogonality of Z(.,JP) to CGRAD.
C
              SUM=ZERO
              IZ=JP
              DO 70 I=1,N
              SUM=SUM+Z(IZ)*CGRAD(I)
              IZ=IZ+N
   70         CONTINUE
              IF (SUM .NE. ZERO) THEN
                  IZ=IPIV*N-N+JP
                  Z(IZ)=Z(IZ)-SUM/CGRAD(IPIV)
              END IF
          END IF
          GO TO 40
      END IF
C
C     Test for linear independence in the proposed new active set.
C
      IF (ZTC(NP) .EQ. ZERO) GOTO 90
      IF (ICON .LE. M) THEN
          SUM=ZERO
          SUMABS=ZERO
          IZ=NP
          DO 80 I=1,N
          TEMP=Z(IZ)*CGRAD(I)
          SUM=SUM+TEMP
          SUMABS=SUMABS+DABS(TEMP)
          IZ=IZ+N
   80     CONTINUE
          IF (DABS(SUM) .LE. RELACC*SUMABS) GOTO 90
      END IF
C
C     Set the new diagonal element of U and return.
C
      U(NP)=ONE/ZTC(NP)
      NACT=NP
   90 RETURN
      END
      SUBROUTINE VE11GD (N,M,A,IA,B,XL,XU,X,IACT,NACT,XBIG,RELACC,TOL,
     1  MEQL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),XBIG(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
C
C     Set VIOL to the greatest relative constraint residual of the first
C       NACT constraints.
C
      VIOL=ZERO
      IF (NACT .GT. MEQL) THEN
          KL=MEQL+1
          DO 20 K=KL,NACT
          J=IACT(K)
          IF (J .LE. M) THEN
              RES=B(J)
              RESABS=DABS(B(J))
              DO 10 I=1,N
              RES=RES-A(I,J)*X(I)
              RESABS=RESABS+DABS(A(I,J)*XBIG(I))
   10         CONTINUE
          ELSE
              JM=J-M
              IF (JM .LE. N) THEN
                  RES=X(JM)-XL(JM)
                  RESABS=XBIG(JM)+DABS(XL(JM))
              ELSE
                  JM=JM-N
                  RES=XU(JM)-X(JM)
                  RESABS=XBIG(JM)+DABS(XU(JM))
              END IF
          END IF
          IF (RES .GT. ZERO) VIOL=DMAX1(VIOL,RES/RESABS)
   20     CONTINUE
      END IF
C
C     Adjust TOL.
C
      TOL=0.1*DMIN1(TOL,VIOL)
      IF (TOL .LE. RELACC+RELACC) THEN
          TOL=RELACC
          DO 30 I=1,N
          XBIG(I)=DABS(X(I))
   30     CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE VE11JD (N,M,A,IA,B,XL,XU,X,IACT,NACT,PAR,G,Z,U,XBIG,
     1  BRES,D,ZTG,RELACC,TOL,STEPCB,SUMRES,MEQL,MSAT,MTOT,INDXBD,
     2  GM,GMNEW,PARNEW,CGRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),PAR(*),G(*),
     1  Z(*),U(*),XBIG(*),BRES(*),D(*),ZTG(*),GM(*),GMNEW(*),PARNEW(*),
     2  CGRAD(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
      IDIFF=MTOT-MSAT
C
C     Calculate and partition the residuals of the inactive constraints,
C       and set the gradient vector when seeking feasibility.
C
      IF (IDIFF .GT. 0) THEN
          DO 10 I=1,N
          G(I)=ZERO
   10     CONTINUE
          SUMRES=ZERO
      END IF
      MSATK=MSAT
      MDEG=NACT
      MSAT=NACT
      KL=MEQL+1
      DO 50 K=KL,MTOT
      J=IACT(K)
C
C     Calculate the residual of the current constraint.
C
      IF (J .LE. M) THEN
          RES=B(J)
          RESABS=DABS(B(J))
          DO 20 I=1,N
          RES=RES-X(I)*A(I,J)
          RESABS=RESABS+DABS(XBIG(I)*A(I,J))
   20     CONTINUE
      ELSE
          JM=J-M
          IF (JM .LE. N) THEN
              RES=X(JM)-XL(JM)
              RESABS=DABS(XBIG(JM))+DABS(XL(JM))
          ELSE
              JM=JM-N
              RES=XU(JM)-X(JM)
              RESABS=DABS(XBIG(JM))+DABS(XU(JM))
          END IF
      END IF
      BRES(J)=RES
C
C     Set TEMP to the relative residual.
C
      TEMP=ZERO
      IF (RESABS .NE. ZERO) TEMP=RES/RESABS
      IF (K .GT. MSATK .AND. TEMP .LT. ZERO) THEN
          IF (TEMP+RELACC .GE. ZERO) THEN
              IF (J .LE. M) THEN
                  SUM=DABS(B(J))
                  DO 30 I=1,N
                  SUM=SUM+DABS(X(I)*A(I,J))
   30             CONTINUE
              ELSE
                  JM=J-M
                  IF (JM .LE. N) THEN
                      SUM=DABS(X(JM))+DABS(XL(JM))
                  ELSE
                      SUM=DABS(X(JM-N))+DABS(XU(JM-N))
                  END IF
              END IF
              IF (DABS(RES) .LE. SUM*RELACC) TEMP=ZERO
          END IF
      END IF
C
C     Place the residual in the appropriate position.
C
      IF (K .LE. NACT) GOTO 50
      IF (K .LE. MSATK .OR. TEMP .GE. ZERO) THEN
          MSAT=MSAT+1
          IF (MSAT .LT. K) THEN
              IACT(K)=IACT(MSAT)
          END IF
          IF (TEMP .GT. TOL) THEN
              IACT(MSAT)=J
          ELSE
              MDEG=MDEG+1
              IACT(MSAT)=IACT(MDEG)
              IACT(MDEG)=J
          END IF
C
C     Update the gradient and SUMRES if the constraint is violated when
C       seeking feasibility.
C
      ELSE
          IF (J .LE. M) THEN
              DO 40 I=1,N
              G(I)=G(I)+A(I,J)
   40         CONTINUE
          ELSE
              J=J-M
              IF (J .LE. N) THEN
                  G(J)=G(J)-ONE
              ELSE
                  G(J-N)=G(J-N)+ONE
              END IF
          END IF
          SUMRES=SUMRES+DABS(RES)
      END IF
   50 CONTINUE
C
C    Seek the next search direction unless VE11JD was called from VE11ED
C       and feasibility has been achieved.
C
      STEPCB=ZERO
      IF (IDIFF .GT. 0 .AND. MSAT .EQ. MTOT) GOTO 60
      CALL VE11OD (N,M,A,IA,IACT,NACT,PAR,G,Z,U,D,ZTG,RELACC,DDOTG,MEQL,
     1  MDEG,GM,GMNEW,PARNEW,CGRAD)
C
C     Calculate the (bound on the) step-length due to the constraints.
C
      IF (DDOTG .LT. ZERO) THEN
          CALL VE11PD (N,M,A,IA,IACT,BRES,D,STEPCB,DDOTG,MDEG,MSAT,
     1      MTOT,INDXBD)
      END IF
      IF (IDIFF .EQ. 0) SUMRES=DDOTG
   60 RETURN
      END
      SUBROUTINE VE11ND (N,M,A,IA,IACT,NACT,Z,U,RELACC,IDROP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),Z(*),U(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
      NM=NACT-1
      IF (IDROP .EQ. NACT) GOTO 60
      ISAVE=IACT(IDROP)
C
C     Cycle through the constraint exchanges that are needed.
C
      DO 50 J=IDROP,NM
      JP=J+1
      ICON=IACT(JP)
      IACT(J)=ICON
C
C     Calculate the (J,JP) element of R.
C
      IF (ICON .LE. M) THEN
          RJJP=ZERO
          IZ=J
          DO 10 I=1,N
          RJJP=RJJP+Z(IZ)*A(I,ICON)
          IZ=IZ+N
   10     CONTINUE
      ELSE
          IBD=ICON-M
          IF (IBD .LE. N) THEN
              IZBD=IBD*N-N
              RJJP=-Z(IZBD+J)
          ELSE
              IBD=IBD-N
              IZBD=IBD*N-N
              RJJP=Z(IZBD+J)
          END IF
      END IF
C
C     Calculate the parameters of the next rotation.
C
      UJP=U(JP)
      TEMP=RJJP*UJP
      DENOM=DABS(TEMP)
      IF (DENOM*RELACC .LT. ONE) DENOM=DSQRT(ONE+DENOM*DENOM)
      WCOS=TEMP/DENOM
      WSIN=ONE/DENOM
C
C     Rotate Z when a bound constraint is promoted.
C
      IZ=J
      IF (ICON .GT. M) THEN
          DO 20 I=1,N
          TEMP=WCOS*Z(IZ+1)-WSIN*Z(IZ)
          Z(IZ)=WCOS*Z(IZ)+WSIN*Z(IZ+1)
          Z(IZ+1)=TEMP
          IZ=IZ+N
   20     CONTINUE
          Z(IZBD+JP)=ZERO
C
C     Rotate Z when an ordinary constraint is promoted.
C
      ELSE
          WPIV=ZERO
          DO 30 I=1,N
          TEMPA=WCOS*Z(IZ+1)
          TEMPB=WSIN*Z(IZ)
          TEMP=DABS(A(I,ICON))*(DABS(TEMPA)+DABS(TEMPB))
          IF (TEMP .GT. WPIV) THEN
              WPIV=TEMP
              IPIV=I
          END IF
          Z(IZ)=WCOS*Z(IZ)+WSIN*Z(IZ+1)
          Z(IZ+1)=TEMPA-TEMPB
          IZ=IZ+N
   30     CONTINUE
C
C     Ensure orthogonality to promoted constraint.
C
          SUM=ZERO
          IZ=JP
          DO 40 I=1,N
          SUM=SUM+Z(IZ)*A(I,ICON)
          IZ=IZ+N
   40     CONTINUE
          IF (SUM .NE. ZERO) THEN
              IZ=IPIV*N-N+JP
              Z(IZ)=Z(IZ)-SUM/A(IPIV,ICON)
          END IF
      END IF
C
C     Set the new diagonal elements of U.
C
      U(JP)=-DENOM*U(J)
      U(J)=UJP/DENOM
   50 CONTINUE
C
C     Return.
C
      IACT(NACT)=ISAVE
   60 NACT=NM
      RETURN
      END
      SUBROUTINE VE11DD (N,M,MEQ,A,IA,B,XU,IACT,MEQL,INFO,Z,U,RELACC,
     1  AM,CGRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XU(*),IACT(*),Z(*),U(*),AM(*),CGRAD(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
C
C     Try to add the next equality constraint to the active set.
C
      DO 50 KEQ=1,MEQ
      IF (MEQL .LT. N) THEN
          NP=MEQL+1
          IACT(NP)=KEQ
          CALL VE11HD (N,M,A,IA,IACT,MEQL,Z,U,RELACC,NP,AM,CGRAD)
          IF (MEQL .EQ. NP) GOTO 50
      END IF
C
C     If linear dependence occurs then find the multipliers of the
C       dependence relation and apply them to the right hand sides.
C
      SUM=B(KEQ)
      SUMABS=DABS(B(KEQ))
      IF (MEQL .GT. 0) THEN
          DO 10 I=1,N
          AM(I)=A(I,KEQ)
   10     CONTINUE
          K=MEQL
   20     VMULT=ZERO
          IZ=K
          DO 30 I=1,N
          VMULT=VMULT+Z(IZ)*AM(I)
          IZ=IZ+N
   30     CONTINUE
          VMULT=VMULT*U(K)
          J=IACT(K)
          IF (J .LE. M) THEN
              DO 40 I=1,N
              AM(I)=AM(I)-VMULT*A(I,J)
   40         CONTINUE
              RHS=B(J)
          ELSE
              JM=J-M-N
              AM(JM)=AM(JM)-VMULT
              RHS=XU(JM)
          END IF
          SUM=SUM-RHS*VMULT
          SUMABS=SUMABS+DABS(RHS*VMULT)
          K=K-1
          IF (K .GE. 1) GOTO 20
      END IF
C
C     Error return if the constraints are inconsistent.
C
      IF (DABS(SUM) .GT. RELACC*SUMABS) THEN
          INFO=5
          GOTO 60
      END IF
   50 CONTINUE
   60 RETURN
      END
      SUBROUTINE VE11OD (N,M,A,IA,IACT,NACT,PAR,G,Z,U,D,ZTG,RELACC,
     1  DDOTG,MEQL,MDEG,GM,GMNEW,PARNEW,CGRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),PAR(*),G(*),Z(*),U(*),D(*),ZTG(*),
     1  GM(*),GMNEW(*),PARNEW(*),CGRAD(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
C
C     Initialize GM and cycle backwards through the active set.
C
   10 DO 20 I=1,N
      GM(I)=G(I)
   20 CONTINUE
      K=NACT
   30 IF (K .GT. 0) THEN
C
C     Set TEMP to the next multiplier, but reduce the active set if
C       TEMP has an unacceptable sign.
C
          TEMP=ZERO
          IZ=K
          DO 40 I=1,N
          TEMP=TEMP+Z(IZ)*GM(I)
          IZ=IZ+N
   40     CONTINUE
          TEMP=TEMP*U(K)
          IF (K .GT. MEQL .AND. TEMP .GT. ZERO) THEN
              CALL VE11ND (N,M,A,IA,IACT,NACT,Z,U,RELACC,K)
              GOTO 10
          END IF
C
C     Update GM using the multiplier that has just been calculated.
C
          J=IACT(K)
          IF (J .LE. M) THEN
              DO 50 I=1,N
              GM(I)=GM(I)-TEMP*A(I,J)
   50         CONTINUE
          ELSE
              JM=J-M
              IF (JM .LE. N) THEN
                  GM(JM)=GM(JM)+TEMP
              ELSE
                  GM(JM-N)=GM(JM-N)-TEMP
              END IF
          END IF
          PAR(K)=TEMP
          K=K-1
          GOTO 30
      END IF
C
C     Calculate the search direction and DDOTG.
C
      DDOTG=ZERO
      IF (NACT .LT. N) THEN
          CALL VE11QD (N,M,A,IA,IACT,NACT,PAR,Z,U,D,ZTG,GM,RELACC,
     1      DDOTGM,MEQL,MDEG,GMNEW,PARNEW,CGRAD)
          IF (DDOTGM .LT. ZERO) THEN
              DO 60 I=1,N
              DDOTG=DDOTG+D(I)*G(I)
   60         CONTINUE
          END IF
      END IF
      RETURN
      END
      SUBROUTINE VE11ED (N,M,A,IA,B,XL,XU,X,IACT,NACT,PAR,INFO,G,Z,
     1  U,XBIG,RELACC,TOL,MEQL,MSAT,MTOT,BRES,D,ZTG,GM,GMNEW,PARNEW,
     2  CGRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),PAR(*),G(*),Z(*),
     1  U(*),XBIG(*),BRES(*),D(*),ZTG(*),GM(*),GMNEW(*),PARNEW(*),
     2  CGRAD(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
C
C     Make the correction to X for the active constraints.
C
      INFO=0
   10 CALL VE11ID (N,M,A,IA,B,XL,XU,X,IACT,NACT,INFO,Z,U,XBIG,RELACC,
     1  TOL,MEQL)
      IF (INFO .GT. 0) MSAT=NACT
      IF (MSAT .EQ. MTOT) GOTO 60
C
C     Try to correct the infeasibility.
C
   20 MSATK=MSAT
      SUMRSK=ZERO
   30 CALL VE11JD (N,M,A,IA,B,XL,XU,X,IACT,NACT,PAR,G,Z,U,XBIG,BRES,
     1  D,ZTG,RELACC,TOL,STEPCB,SUMRES,MEQL,MSAT,MTOT,INDXBD,GM,GMNEW,
     2  PARNEW,CGRAD)
C
C     Include the new constraint in the active set.
C
      IF (STEPCB .GT. ZERO) THEN
          DO 40 I=1,N
          X(I)=X(I)+STEPCB*D(I)
          XBIG(I)=DMAX1(XBIG(I),DABS(X(I)))
   40     CONTINUE
          CALL VE11HD (N,M,A,IA,IACT,NACT,Z,U,RELACC,INDXBD,GMNEW,CGRAD)
      END IF
C
C     Test whether to continue the search for feasibility.
C
      IF (MSAT .LT. MTOT) THEN
          IF (STEPCB .EQ. ZERO) GOTO 50
          IF (MSATK .LT. MSAT) GOTO 20
          IF (SUMRSK .EQ. ZERO .OR. SUMRES .LT. SUMRSK) THEN
              SUMRSK=SUMRES
              ITEST=0
          END IF
          ITEST=ITEST+1
          IF (ITEST .LE. 2) GOTO 30
C
C     Reduce TOL if it may be too large to allow feasibility.
C
   50     IF (TOL .GT. RELACC) THEN
              CALL VE11GD (N,M,A,IA,B,XL,XU,X,IACT,NACT,XBIG,RELACC,
     1          TOL,MEQL)
              GOTO 10
          END IF
      END IF
   60 RETURN
      END
      SUBROUTINE VE11KD (N,M,A,IA,IACT,NACT,PAR,G,RESKT,Z,U,BRES,RELAXF,
     1  MEQL,SSQKT,PARW,RESKTW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),PAR(*),G(*),RESKT(*),Z(*),U(*),
     1  BRES(*),PARW(*),RESKTW(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
C
C     Calculate the Lagrange parameters and the residual vector.
C
      DO 10 I=1,N
      RESKT(I)=G(I)
   10 CONTINUE
      IF (NACT .GT. 0) THEN
          ICASE=0
   20     DO 50 KK=1,NACT
          K=NACT+1-KK
          J=IACT(K)
          TEMP=ZERO
          IZ=K
          DO 30 I=1,N
          TEMP=TEMP+Z(IZ)*RESKT(I)
          IZ=IZ+N
   30     CONTINUE
          TEMP=TEMP*U(K)
          IF (ICASE .EQ. 0) PAR(K)=ZERO
          IF (K .LE. MEQL .OR. PAR(K)+TEMP .LT. ZERO) THEN
              PAR(K)=PAR(K)+TEMP
          ELSE
              TEMP=-PAR(K)
              PAR(K)=ZERO
          END IF
          IF (TEMP .NE. ZERO) THEN
              IF (J .LE. M) THEN
                  DO 40 I=1,N
                  RESKT(I)=RESKT(I)-TEMP*A(I,J)
   40             CONTINUE
              ELSE
                  JM=J-M
                  IF (JM .LE. N) THEN
                      RESKT(JM)=RESKT(JM)+TEMP
                  ELSE
                      RESKT(JM-N)=RESKT(JM-N)-TEMP
                  END IF
              END IF
          END IF
   50     CONTINUE
C
C     Calculate the sum of squares of the KT residual vector.
C
          SSQKT=ZERO
          IF (NACT .EQ. N) GOTO 130
          DO 60 I=1,N
          SSQKT=SSQKT+RESKT(I)**2
   60     CONTINUE
C
C     Apply iterative refinement to the residual vector.
C
          IF (ICASE .EQ. 0) THEN
              ICASE=1
              DO 70 K=1,NACT
              PARW(K)=PAR(K)
   70         CONTINUE
              DO 80 I=1,N
              RESKTW(I)=RESKT(I)
   80         CONTINUE
              SSQKTW=SSQKT
              GOTO 20
          END IF
C
C     Undo the iterative refinement if it does not reduce SSQKT.
C
          IF (SSQKTW .LT. SSQKT) THEN
              DO 90 K=1,NACT
              PAR(K)=PARW(K)
   90         CONTINUE
              DO 100 I=1,N
              RESKT(I)=RESKTW(I)
  100         CONTINUE
              SSQKT=SSQKTW
          END IF
C
C     Calculate SSQKT when there are no active constraints.
C
      ELSE
          SSQKT=ZERO
          DO 110 I=1,N
          SSQKT=SSQKT+G(I)**2
  110     CONTINUE
      END IF
C
C     Predict the reduction in F if one corrects any positive residuals
C       of active inequality constraints.
C
      RELAXF=ZERO
      IF (MEQL. LT. NACT) THEN
          KL=MEQL+1
          DO 120 K=KL,NACT
          J=IACT(K)
          IF (BRES(J) .GT. ZERO) THEN
              RELAXF=RELAXF-PAR(K)*BRES(J)
          END IF
  120     CONTINUE
      END IF
  130 RETURN
      END
      SUBROUTINE VE11LD (N,X,G,D,XS,GS,RELACC,STEPCB,DDOTG,F,STEP,
     1  NFVALS,NFMAX,GOPT,FGCALC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),G(*),D(*),XS(*),GS(*),GOPT(*)
      EXTERNAL FGCALC
      DATA ZERO/0.0D0/,ONE/1.0D0/
C
C     Initialization.
C
      RELINT=0.9
      ICOUNT=0
      RATIO=-ONE
      DO 10 I=1,N
      XS(I)=X(I)
      GS(I)=G(I)
      GOPT(I)=G(I)
      IF (D(I) .NE. ZERO) THEN
          TEMP=DABS(X(I)/D(I))
          IF (RATIO .LT. ZERO .OR. TEMP .LT. RATIO) RATIO=TEMP
      END IF
   10 CONTINUE
      STEP=DMIN1(ONE,STEPCB)
      STPMIN=DMAX1(RELACC*RATIO,1.0E-12*STEP)
      STEP=DMAX1(STPMIN,STEP)
      SBASE=ZERO
      FBASE=F
      DDOTGB=DDOTG
      STPLOW=ZERO
      FLOW=F
      DGLOW=DDOTG
      STPHGH=ZERO
      STPOPT=ZERO
      FOPT=F
      DGOPT=DABS(DDOTG)
C
C     Calculate another function and gradient value.
C
   20 DO 30 I=1,N
      X(I)=XS(I)+STEP*D(I)
   30 CONTINUE
      CALL FGCALC (N,X,F,G)
      ICOUNT=ICOUNT+1
      DGMID=ZERO
      DO 40 I=1,N
      DGMID=DGMID+D(I)*G(I)
   40 CONTINUE
      IF (F .LE. FOPT) THEN
          IF (F .LT. FOPT .OR. DABS(DGMID) .LT. DGOPT) THEN
              STPOPT=STEP
              FOPT=F
              DO 50 I=1,N
              GOPT(I)=G(I)
   50         CONTINUE
              DGOPT=DABS(DGMID)
          END IF
      END IF
      IF (NFVALS+ICOUNT .EQ. NFMAX) GOTO 70
C
C      Modify the bounds on the steplength or convergence.
C
      IF (F .GE. FBASE+0.1*(STEP-SBASE)*DDOTGB) THEN
          IF (STPHGH .GT. ZERO .OR. F .GT. FBASE .OR. DGMID .GT.
     1      0.5*DDOTG) THEN
              STPHGH=STEP
              FHGH=F
              DGHGH=DGMID
              GOTO 60
          END IF
          SBASE=STEP
          FBASE=F
          DDOTGB=DGMID
      END IF
      IF (DGMID .GE. 0.7*DDOTGB) GOTO 70
      STPLOW=STEP
      FLOW=F
      DGLOW=DGMID
   60 IF (STPHGH .GT. ZERO .AND. STPLOW .GE. RELINT*STPHGH) GOTO 70
C
C     Calculate the next step length or end the iterations.
C
      IF (STPHGH .EQ. ZERO) THEN
          IF (STEP .EQ. STEPCB) GOTO 70
          TEMP=10.0
          IF (DGMID .GT. 0.9*DDOTG) TEMP=DDOTG/(DDOTG-DGMID)
          STEP=DMIN1(TEMP*STEP,STEPCB)
          GOTO 20
      ELSE IF (ICOUNT .EQ. 1 .OR. STPLOW .GT. ZERO) THEN
          DGKNOT=2.0*(FHGH-FLOW)/(STPHGH-STPLOW)-0.5*(DGLOW+DGHGH)
          IF (DGKNOT .GE. ZERO) THEN
              RATIO=DMAX1(0.1D0,0.5*DGLOW/(DGLOW-DGKNOT))
          ELSE
              RATIO=(0.5*DGHGH-DGKNOT)/(DGHGH-DGKNOT)
          END IF
          STEP=STPLOW+RATIO*(STPHGH-STPLOW)
          GOTO 20
      ELSE
          STEP=0.1*STEP
          IF (STEP .GE. STPMIN) GOTO 20
      END IF
C
C     Return from subroutine.
C
   70 IF (STEP .NE. STPOPT) THEN
          STEP=STPOPT
          F=FOPT
          DO 80 I=1,N
          X(I)=XS(I)+STEP*D(I)
          G(I)=GOPT(I)
   80     CONTINUE
      END IF
      NFVALS=NFVALS+ICOUNT
      RETURN
      END
      SUBROUTINE VE11BD (N,M,MEQ,A,IA,B,XL,XU,X,ACC,IACT,NACT,PAR,
     1  IPRINT,INFO,G,Z,U,XBIG,RESKT,BRES,D,ZTG,GM,XS,GS,FGCALC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),PAR(*),G(*),
     1  Z(*),U(*),XBIG(*),RESKT(*),BRES(*),D(*),ZTG(*),GM(*),XS(*),
     2  GS(*)
      COMMON /VE11UD/LP
      EXTERNAL FGCALC
      DATA ZERO/0.0D0/,ONE/1.0D0/
C
C     Initialize ZZNORM, ITERC, NFVALS and NFMAX.
C
      ZZNORM=-ONE
      ITERC=0
      NFVALS=0
      NFMAX=0
      IF (INFO .GT. 0) NFMAX=INFO
C
C     Check the bounds on N, M and MEQ.
C
      INFO=4
      IF (MAX0(1-N,1-M,MEQ*(MEQ-M)) .GT. 0) THEN
          IF (IPRINT .NE. 0 .AND. LP.GT.0)WRITE(LP,1010)
 1010     FORMAT (/5X,'ERROR RETURN FROM VE11AD BECAUSE A CONDITION',
     1      ' ON N, M OR MEQ IS VIOLATED')
          GOTO 40
      END IF
C
C     Initialize RELACC, Z, U and TOL.
C
      CALL VE11CD (N,M,XL,XU,X,IACT,MEQL,INFO,Z,U,XBIG,RELACC)
      TOL=DMAX1(0.01D0,10.0*RELACC)
      IF (INFO .EQ. 4) THEN
          IF (IPRINT .NE. 0 .AND. LP.GT.0)WRITE(LP,1020)
 1020     FORMAT (/5X,'ERROR RETURN FROM VE11AD BECAUSE A LOWER',
     1      ' BOUND EXCEEDS AN UPPER BOUND')
          GOTO 40
      END IF
C
C     Add any equality constraints to the active set.
C
      IF (MEQ .GT. 0) THEN
          CALL VE11DD (N,M,MEQ,A,IA,B,XU,IACT,MEQL,INFO,Z,U,RELACC,XS,
     1      GS)
          IF (INFO .EQ. 5) THEN
              IF (IPRINT .NE. 0 .AND. LP.GT.0)WRITE(LP,1030)
 1030         FORMAT (/5X,'ERROR RETURN FROM VE11AD BECAUSE THE',
     1          ' EQUALITY CONSTRAINTS ARE INCONSISTENT')
              GOTO 40
          END IF
      END IF
      NACT=MEQL
      MSAT=MEQL
C
C     Add the bounds to the list of constraints.
C
      MTOT=NACT
      DO 10 I=1,N
      IF (XL(I) .LT. XU(I)) THEN
          MTOT=MTOT+2
          IACT(MTOT-1)=M+I
          IACT(MTOT)=M+N+I
      END IF
   10 CONTINUE
C
C     Try to satisfy the bound constraints.
C
      CALL VE11ED (N,M,A,IA,B,XL,XU,X,IACT,NACT,PAR,INFO,G,Z,U,XBIG,
     1  RELACC,TOL,MEQL,MSAT,MTOT,BRES,D,ZTG,GM,RESKT,XS,GS)
      IF (MSAT .LT. MTOT) THEN
          IF (IPRINT .NE. 0 .AND. LP.GT.0)WRITE(LP,1040)
 1040     FORMAT (/5X,'ERROR RETURN FROM VE11AD BECAUSE THE',
     1      ' EQUALITIES AND BOUNDS ARE INCONSISTENT')
          INFO=6
          GOTO 40
      END IF
C
C     Add the ordinary inequalities to the list of constraints.
C
      IF (M .GT. MEQ) THEN
          MP=MEQ+1
          DO 20 K=MP,M
          MTOT=MTOT+1
          IACT(MTOT)=K
   20     CONTINUE
      END IF
C
C     Correct any constraint violations.
C
   30 CALL VE11ED (N,M,A,IA,B,XL,XU,X,IACT,NACT,PAR,INFO,G,Z,U,XBIG,
     1  RELACC,TOL,MEQL,MSAT,MTOT,BRES,D,ZTG,GM,RESKT,XS,GS)
      IF (MSAT .LT. MTOT) THEN
          IF (IPRINT .NE. 0 .AND. LP.GT.0)WRITE(LP,1050)
 1050     FORMAT (/5X,'ERROR RETURN FROM VE11AD BECAUSE THE',
     1      ' CONSTRAINTS ARE INCONSISTENT')
          INFO=7
          GOTO 40
      ELSE IF (MEQL .EQ. N) THEN
          IF (IPRINT .NE. 0 .AND. LP.GT.0)WRITE(LP,1060)
 1060     FORMAT (/5X,'VE11AD FINDS THAT THE VARIABLES ARE',
     1      ' DETERMINED BY THE EQUALITY CONSTRAINTS')
          GOTO 40
      END IF
C
C     Minimize the objective function in the case when constraints are
C       treated as degenerate if their residuals are less than TOL.
C
      CALL VE11FD (N,M,A,IA,B,XL,XU,X,ACC,IACT,NACT,PAR,IPRINT,INFO,G,Z,
     1  U,XBIG,RELACC,ZZNORM,TOL,MEQL,MTOT,ITERC,NFVALS,NFMAX,RESKT,
     2  BRES,D,ZTG,GM,XS,GS,FGCALC)
C
C     Reduce TOL if necessary.
C
      IF (TOL .GT. RELACC .AND. NACT .GT. 0) THEN
          IF (NFVALS .NE. NFMAX) THEN
              CALL VE11GD (N,M,A,IA,B,XL,XU,X,IACT,NACT,XBIG,RELACC,TOL,
     1          MEQL)
              GOTO 30
          ELSE
              INFO=8
          END IF
      END IF
      IF (IPRINT .NE. 0) THEN
          IF (INFO .EQ. 1 .AND. LP.GT.0)WRITE(LP,1070)
 1070     FORMAT (/5X,'VE11AD HAS ACHIEVED THE REQUIRED ACCURACY')
          IF (INFO .EQ. 2 .AND. LP.GT.0)WRITE(LP,1080)
 1080     FORMAT (/5X,'VE11AD CAN MAKE NO FURTHER PROGRESS BECAUSE',
     1      ' OF ROUNDING ERRORS')
          IF (INFO .EQ. 3 .AND. LP.GT.0)WRITE(LP,1090)
 1090     FORMAT (/5X,'VE11AD CAN MAKE NO FURTHER PROGRESS BECAUSE',
     1      ' F WILL NOT DECREASE ANY MORE')
          IF (INFO .EQ. 8 .AND. LP.GT.0)WRITE(LP,1100)
 1100     FORMAT (/5X,'VE11AD HAS REACHED THE GIVEN LIMIT ON THE',
     1      ' NUMBER OF CALLS OF FGCALC')
      END IF
   40 RETURN
      END
      SUBROUTINE VE11FD (N,M,A,IA,B,XL,XU,X,ACC,IACT,NACT,PAR,IPRINT,
     1  INFO,G,Z,U,XBIG,RELACC,ZZNORM,TOL,MEQL,MTOT,ITERC,NFVALS,
     2  NFMAX,RESKT,BRES,D,ZTG,GM,XS,GS,FGCALC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),PAR(*),G(*),Z(*),
     1  U(*),XBIG(*),RESKT(*),BRES(*),D(*),ZTG(*),GM(*),XS(*),GS(*)
      COMMON /VE11UD/LP
      EXTERNAL FGCALC
      SAVE F
      DATA ZERO/0.0D0/,ONE/1.0D0/
C
C     Initialize the minimization calculation.
C
      MSAT=MTOT
      ITERK=ITERC
      NFVALK=NFVALS
      IF (NFVALS .EQ. 0 .OR. INFO .EQ. 1) THEN
          CALL FGCALC (N,X,F,G)
          NFVALS=NFVALS+1
      END IF
      FPREV=DABS(F+F+ONE)
      ITERP=-1
      IF (IPRINT .NE. 0) THEN
          IF(LP.GT.0)WRITE(LP,1000) TOL
 1000     FORMAT (/5X,'NEW VALUE OF TOL =',1PE13.5)
          ITERP=ITERC
      END IF
C
C     Calculate the next search direction.
C
   10 CALL VE11JD (N,M,A,IA,B,XL,XU,X,IACT,NACT,PAR,G,Z,U,XBIG,BRES,D,
     1  ZTG,RELACC,TOL,STEPCB,DDOTG,MEQL,MSAT,MTOT,INDXBD,GM,RESKT,XS,
     2  GS)
C
C     Calculate the Kuhn Tucker residual vector.
C
      CALL VE11KD (N,M,A,IA,IACT,NACT,PAR,G,RESKT,Z,U,BRES,RELAXF,MEQL,
     1  SSQKT,XS,GS)
C
C     Test for convergence.
C
      IF (SSQKT .LE. ACC*ACC) THEN
          INFO=1
          GOTO 70
      END IF
      IF (DDOTG .GE. ZERO) THEN
          INFO=2
          GOTO 70
      END IF
C
C     Test for termination due to no decrease in F.
C
      IF (F .GE. FPREV) THEN
          IF (TOL .EQ. RELACC .OR. NACT .EQ. 0) THEN
              IF (DIFF .GT. ZERO) GOTO 20
          END IF
          INFO=3
          GOTO 70
      END IF
   20 DIFF=FPREV-F
      FPREV=F
C
C     Test that more calls of FGCALC are allowed.
C
      IF (NFVALS .EQ. NFMAX) THEN
          INFO=8
          GOTO 70
      END IF
C
C     Test whether to reduce TOL and to provide printing.
C
      IF (TOL .GT. RELACC .AND. ITERC .GT. ITERK .AND.
     1  0.1*RELAXF .GE. DMAX1(DIFF,-0.5*DDOTG)) GOTO 70
      IF (ITERP .EQ. ITERC) THEN
          ITERP=ITERC+IABS(IPRINT)
          GOTO 80
      END IF
C
C     Calculate the step along the search direction.
C
   40 ITERC=ITERC+1
      CALL VE11LD (N,X,G,D,XS,GS,RELACC,STEPCB,DDOTG,F,STEP,NFVALS,
     1  NFMAX,BRES,FGCALC)
      IF (STEP .EQ. ZERO) THEN
          INFO=3
          SUM=ZERO
          DO 50 I=1,N
          SUM=SUM+DABS(D(I)*GS(I))
   50     CONTINUE
          IF (DDOTG+RELACC*SUM .GE. ZERO) INFO=2
          GOTO 70
      END IF
C
C     Revise XBIG.
C
      DO 60 I=1,N
      XBIG(I)=DMAX1(XBIG(I),DABS(X(I)))
   60 CONTINUE
C
C     Revise the second derivative approximation.
C
      CALL VE11MD (N,X,NACT,G,Z,ZTG,XS,GS,ZZNORM)
C
C     Add a constraint to the active set if it restricts the step.
C
      IF (STEP .EQ. STEPCB) THEN
          K=IACT(INDXBD)
          IF (K .GT. M) THEN
              K=K-M
              IF (K .LE. N) THEN
                  X(K)=XL(K)
              ELSE
                  X(K-N)=XU(K-N)
              END IF
          END IF
          CALL VE11HD (N,M,A,IA,IACT,NACT,Z,U,RELACC,INDXBD,XS,GS)
      END IF
      GOTO 10
C
C     Printing from the subroutine.
C
   70 IF (IPRINT .EQ. 0) GOTO 90
      ITERP=-1
   80 IF(LP.GT.0)WRITE(LP,1010) ITERC,NFVALS,F
 1010 FORMAT (/5X,'ITERS =',I4,5X,'F.VALS =',I4,5X,'F =',1PE15.7)
      IF (NFVALS .GT. NFVALK) THEN
          IF(LP.GT.0)WRITE(LP,1020) (X(I),I=1,N)
 1020     FORMAT ('  X =',(1P,5E14.5))
          IF(LP.GT.0)WRITE(LP,1030) (G(I),I=1,N)
 1030     FORMAT ('  G =',(1P,5E14.5))
      ELSE
          IF(LP.GT.0)WRITE(LP,1040)
 1040     FORMAT (5X,'NO CHANGE TO X AND G SINCE PREVIOUS OUTPUT')
      END IF
      IF (IPRINT .LT. 0) THEN
          IF (NACT .EQ. 0) THEN
              IF(LP.GT.0)WRITE(LP,1050)
 1050         FORMAT (5X,'NO ACTIVE CONSTRAINTS')
          ELSE
              IF(LP.GT.0)WRITE(LP,1060) (IACT(I),I=1,NACT)
 1060         FORMAT (' IA =',(14I5))
              IF(LP.GT.0)WRITE(LP,1070) (PAR(I),I=1,NACT)
 1070         FORMAT (' LP =',(1P,5E14.5))
          END IF
          IF (NACT .EQ. N) THEN
              IF(LP.GT.0)WRITE(LP,1080)
 1080         FORMAT (5X,'KT RESIDUAL VECTOR IS ZERO')
          ELSE
              IF(LP.GT.0)WRITE(LP,1090)(RESKT(I),I=1,N)
 1090         FORMAT (' KT =',(1P,5E14.5))
          END IF
      END IF
      IF (ITERP .GE. 0) GOTO 40
   90 RETURN
      END
       SUBROUTINE VE11SD (N,M,A,IA,IACT,NACT,Z,U,D,RELACC,MDEG,ZZDIAG,
     1  GMNEW,CGRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),Z(*),U(*),D(*),ZZDIAG(*),GMNEW(*),
     1  CGRAD(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
C
C     Initialization.
C
      NP=NACT+1
      KHIGH=MDEG
      IZ=0
      DO 20 I=1,N
      ZZDIAG(I)=ZERO
      DO 10 J=NP,N
      ZZDIAG(I)=ZZDIAG(I)+Z(IZ+J)**2
   10 CONTINUE
      IZ=IZ+N
   20 CONTINUE
C
C     Calculate the scalar products of D with its constraints.
C
   30 CVMAX=ZERO
      DO 50 K=NP,KHIGH
      J=IACT(K)
      IF (J .LE. M) THEN
          SUM=ZERO
          SUMABS=ZERO
          SUMD=ZERO
          DO 40 I=1,N
          TEMP=D(I)*A(I,J)
          SUM=SUM+TEMP
          SUMABS=SUMABS+DABS(TEMP)
          SUMD=SUMD+ZZDIAG(I)*A(I,J)**2
   40     CONTINUE
      ELSE
          JM=J-M
          IF (JM .LE. N) THEN
              SUM=-D(JM)
          ELSE
              JM=JM-N
              SUM=D(JM)
          END IF
          SUMABS=DABS(SUM)
          SUMD=ZZDIAG(JM)
      END IF
C
C     Pick out the most violated constraint, or return if the
C       violation is negligible.
C
      IF (SUM .GT. RELACC*SUMABS) THEN
          CVIOL=SUM*SUM/SUMD
          IF (CVIOL .GT. CVMAX) THEN
              CVMAX=CVIOL
              IADD=K
              SAVSUM=SUM
              SAVABS=SUMABS
          END IF
      END IF
   50 CONTINUE
      IF (CVMAX .LE. ZERO) GOTO 140
      IF (NACT .EQ. 0) GOTO 120
C
C     Set GMNEW to the gradient of the most violated constraint.
C
      J=IACT(IADD)
      IF (J .LE. M) THEN
          JMV=0
          DO 60 I=1,N
          GMNEW(I)=A(I,J)
   60     CONTINUE
      ELSE
          JMV=J-M
          DO 70 I=1,N
          GMNEW(I)=ZERO
   70     CONTINUE
          IF (JMV .LE. N) THEN
              GMNEW(JMV)=-ONE
          ELSE
              JMV=JMV-N
              GMNEW(JMV)=ONE
          END IF
      END IF
C
C     Modify GMNEW for the next active constraint.
C
      K=NACT
   80 TEMP=ZERO
      IZ=K
      DO 90 I=1,N
      TEMP=TEMP+Z(IZ)*GMNEW(I)
      IZ=IZ+N
   90 CONTINUE
      TEMP=TEMP*U(K)
      J=IACT(K)
      IF (J .LE. M) THEN
          DO 100 I=1,N
          GMNEW(I)=GMNEW(I)-TEMP*A(I,J)
  100     CONTINUE
      ELSE
          JM=J-M
          IF (JM .LE. N) THEN
              GMNEW(JM)=GMNEW(JM)+TEMP
          ELSE
              GMNEW(JM-N)=GMNEW(JM-N)-TEMP
          END IF
      END IF
C
C     Revise the values of SAVSUM and SAVABS.
C
      SUM=ZERO
      SUMABS=ZERO
      DO 110 I=1,N
      TEMP=D(I)*GMNEW(I)
      SUM=SUM+TEMP
      SUMABS=SUMABS+DABS(TEMP)
  110 CONTINUE
      SAVSUM=DMIN1(SAVSUM,SUM)
      SAVABS=DMAX1(SAVABS,SUMABS)
      K=K-1
      IF (K .GE. 1) GOTO 80
C
C     Add the new constraint to the active set if the constraint
C       violation is still significant.
C
      IF (JMV .GT. 0) D(JMV)=ZERO
      IF (SAVSUM .LE. RELACC*SAVABS) GOTO 130
  120 K=NACT
      CALL VE11HD (N,M,A,IA,IACT,NACT,Z,U,RELACC,IADD,GMNEW,CGRAD)
      IF (NACT .GT. K) GOTO 140
C
C     Seek another constraint violation.
C
      IADD=NP
  130 IF (NP .LT. KHIGH) THEN
          K=IACT(KHIGH)
          IACT(KHIGH)=IACT(IADD)
          IACT(IADD)=K
          KHIGH=KHIGH-1
          GOTO 30
      END IF
  140 RETURN
      END
      SUBROUTINE VE11ID (N,M,A,IA,B,XL,XU,X,IACT,NACT,INFO,Z,U,XBIG,
     1  RELACC,TOL,MEQL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),B(*),XL(*),XU(*),X(*),IACT(*),Z(*),U(*),
     1  XBIG(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
      IF (NACT .EQ. 0) GOTO 50
      DO 30 K=1,NACT
C
C     Calculate the next constraint residual.
C
      J=IACT(K)
      IF (J .LE. M) THEN
          RES=B(J)
          RESABS=DABS(B(J))
          RESBIG=RESABS
          DO 10 I=1,N
          TEMPA=A(I,J)
          TEMP=TEMPA*X(I)
          RES=RES-TEMP
          RESABS=RESABS+DABS(TEMP)
          RESBIG=RESBIG+DABS(TEMPA)*XBIG(I)
   10     CONTINUE
      ELSE
          JX=J-M
          IF (JX .LE. N) THEN
              RES=X(JX)-XL(JX)
              RESABS=DABS(X(JX))+DABS(XL(JX))
              RESBIG=XBIG(JX)+DABS(XL(JX))
              SAVEX=XL(JX)
          ELSE
              JX=JX-N
              RES=XU(JX)-X(JX)
              RESABS=DABS(X(JX))+DABS(XU(JX))
              RESBIG=XBIG(JX)+DABS(XU(JX))
              SAVEX=XU(JX)
          END IF
      END IF
C
C     Shift X if necessary.
C
      IF (RES .NE. ZERO) THEN
          TEMP=RES/RESABS
          IF (K .LE. MEQL) TEMP=-DABS(TEMP)
          IF (TOL .EQ. RELACC .OR. TEMP+RELACC .LT. ZERO) THEN
              INFO=1
              SCALE=RES*U(K)
              IZ=K
              DO 20 I=1,N
              X(I)=X(I)+SCALE*Z(IZ)
              IZ=IZ+N
              XBIG(I)=DMAX1(XBIG(I),DABS(X(I)))
   20         CONTINUE
              IF (J .GT. M) X(JX)=SAVEX
C
C     Else flag a constraint deletion if necessary.
C
          ELSE IF (RES/RESBIG .GT. TOL) THEN
              IACT(K)=-IACT(K)
          END IF
      END IF
   30 CONTINUE
C
C     Delete any flagged constraints and then return.
C
      IDROP=NACT
  40  IF (IACT(IDROP) .LT. 0) THEN
          IACT(IDROP)=-IACT(IDROP)
          CALL VE11ND (N,M,A,IA,IACT,NACT,Z,U,RELACC,IDROP)
      END IF
      IDROP=IDROP-1
      IF (IDROP .GT. MEQL) GOTO 40
   50 RETURN
      END
      SUBROUTINE VE11QD (N,M,A,IA,IACT,NACT,PAR,Z,U,D,ZTG,GM,RELACC,
     1  DDOTGM,MEQL,MDEG,GMNEW,PARNEW,CGRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),PAR(*),Z(*),U(*),D(*),ZTG(*),GM(*),
     1  GMNEW(*),PARNEW(*),CGRAD(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
      MP=MEQL+1
      DTEST=ZERO
C
C     Calculate the search direction and branch if it is not downhill.
C
   10 CALL VE11RD (N,NACT,Z,D,ZTG,GM,RELACC,DDOTGM)
      IF (DDOTGM .EQ. ZERO) GOTO 120
C
C     Branch if there is no need to consider any degenerate constraints.
C     The test gives termination if two consecutive additions to the
C       active set fail to increase the predicted new value of F.
C
      IF (NACT .EQ. MDEG) GOTO 120
      NP=NACT+1
      SUM=ZERO
      DO 20 J=NP,N
      SUM=SUM+ZTG(J)**2
   20 CONTINUE
      IF (DTEST .GT. ZERO .AND. SUM .GE. DTEST) THEN
          IF (ITEST .EQ. 1) GOTO 120
          ITEST=1
      ELSE
          DTEST=SUM
          ITEST=0
      END IF
C
C     Add a constraint to the active set if there are any significant
C       violations of degenerate constraints.
C
      K=NACT
      CALL VE11SD (N,M,A,IA,IACT,NACT,Z,U,D,RELACC,MDEG,GMNEW,PARNEW,
     1  CGRAD)
      IF (NACT .EQ. K) GOTO 120
      PAR(NACT)=ZERO
C
C     Calculate the new reduced gradient and Lagrange parameters.
C
   30 DO 40 I=1,N
      GMNEW(I)=GM(I)
   40 CONTINUE
      K=NACT
   50 TEMP=ZERO
      IZ=K
      DO 60 I=1,N
      TEMP=TEMP+Z(IZ)*GMNEW(I)
      IZ=IZ+N
   60 CONTINUE
      TEMP=TEMP*U(K)
      PARNEW(K)=PAR(K)+TEMP
      IF (K .EQ. NACT) PARNEW(K)=DMIN1(PARNEW(K),ZERO)
      J=IACT(K)
      IF (J .LE. M) THEN
          DO 70 I=1,N
          GMNEW(I)=GMNEW(I)-TEMP*A(I,J)
   70     CONTINUE
      ELSE
          JM=J-M
          IF (JM .LE. N) THEN
              GMNEW(JM)=GMNEW(JM)+TEMP
          ELSE
              GMNEW(JM-N)=GMNEW(JM-N)-TEMP
          END IF
      END IF
      K=K-1
      IF (K .GT. MEQL) GOTO 50
C
C     Set RATIO for linear interpolation between PAR and PARNEW.
C
      RATIO=ZERO
      IF (MP .LT. NACT) THEN
          KU=NACT-1
          DO 80 K=MP,KU
          IF (PARNEW(K) .GT. ZERO) THEN
              RATIO=PARNEW(K)/(PARNEW(K)-PAR(K))
              IDROP=K
          END IF
   80     CONTINUE
      END IF
C
C     Apply the linear interpolation.
C
      THETA=ONE-RATIO
      DO 90 K=MP,NACT
      PAR(K)=DMIN1(THETA*PARNEW(K)+RATIO*PAR(K),ZERO)
   90 CONTINUE
      DO 100 I=1,N
      GM(I)=THETA*GMNEW(I)+RATIO*GM(I)
  100 CONTINUE
C
C     Drop a constraint if RATIO is positive.
C
      IF (RATIO .GT. ZERO) THEN
          CALL VE11ND (N,M,A,IA,IACT,NACT,Z,U,RELACC,IDROP)
          DO 110 K=IDROP,NACT
          PAR(K)=PAR(K+1)
  110     CONTINUE
          GOTO 30
      END IF
C
C     Return if there is no freedom for a new search direction.
C
      IF (NACT .LT. N) GOTO 10
      DDOTGM=ZERO
  120 RETURN
      END
      SUBROUTINE VE11RD (N,NACT,Z,D,ZTG,GM,RELACC,DDOTGM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(*),D(*),ZTG(*),GM(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
      DDOTGM=ZERO
      IF (NACT .GE. N) GOTO 60
C
C     Premultiply GM by the transpose of Z.
C
      NP=NACT+1
      DO 20 J=NP,N
      SUM=ZERO
      SUMABS=ZERO
      IZ=J
      DO 10 I=1,N
      TEMP=Z(IZ)*GM(I)
      SUM=SUM+TEMP
      SUMABS=SUMABS+DABS(TEMP)
      IZ=IZ+N
   10 CONTINUE
      IF (DABS(SUM) .LE. RELACC*SUMABS) SUM=ZERO
      ZTG(J)=SUM
   20 CONTINUE
C
C     Form D by premultiplying ZTG by -Z.
C
      IZ=0
      DO 40 I=1,N
      SUM=ZERO
      SUMABS=ZERO
      DO 30 J=NP,N
      TEMP=Z(IZ+J)*ZTG(J)
      SUM=SUM-TEMP
      SUMABS=SUMABS+DABS(TEMP)
   30 CONTINUE
      IF (DABS(SUM) .LE. RELACC*SUMABS) SUM=ZERO
      D(I)=SUM
      IZ=IZ+N
   40 CONTINUE
C
C     Test that the search direction is downhill.
C
      SUMABS=ZERO
      DO 50 I=1,N
      TEMP=D(I)*GM(I)
      DDOTGM=DDOTGM+TEMP
      SUMABS=SUMABS+DABS(TEMP)
   50 CONTINUE
      IF (DDOTGM+RELACC*SUMABS .GE. ZERO) DDOTGM=ZERO
   60 RETURN
      END
      SUBROUTINE VE11PD (N,M,A,IA,IACT,BRES,D,STEPCB,DDOTG,MDEG,MSAT,
     1  MTOT,INDXBD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IA,*),IACT(*),BRES(*),D(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/

C    Modified 29/9/92 by John Reid to avoid branches to within IF
C    constructs.
C    The IF constructs at labels 10 and 40 made into IF statements.

C
C    Set steps to constraint boundaries and find the least positive one.
C
      IFLAG=0
      STEPCB=ZERO
      INDXBD=0
      K=MDEG
   10 K=K+1
      IF (K .GT. MTOT) GO TO 40
C
C     Form the scalar product of D with the current constraint normal.
C
   20     J=IACT(K)
          IF (J .LE. M) THEN
              SP=ZERO
              DO 30 I=1,N
              SP=SP+D(I)*A(I,J)
   30         CONTINUE
          ELSE
              JM=J-M
              IF (JM .LE. N) THEN
                  SP=-D(JM)
              ELSE
                  SP=D(JM-N)
              END IF
          END IF
C
C     The next branch is taken if label 20 was reached via label 40.
C
          IF (IFLAG .EQ. 1) GOTO 50
C
C     Set BRES(J) to indicate the status of the j-th constraint.
C
          IF (SP*BRES(J) .LE. ZERO) THEN
              BRES(J)=ZERO
          ELSE
              BRES(J)=BRES(J)/SP
              IF (STEPCB .EQ. ZERO .OR. BRES(J) .LT. STEPCB) THEN
                  STEPCB=BRES(J)
                  INDXBD=K
              END IF
          END IF
          GO TO 10
C
C     Try to pass through the boundary of a violated constraint.
C
   40 IF (INDXBD .LE. MSAT) RETURN
          IFLAG=1
          K=INDXBD
          GOTO 20
   50     MSAT=MSAT+1
          IACT(INDXBD)=IACT(MSAT)
          IACT(MSAT)=J
          BRES(J)=ZERO
          INDXBD=MSAT
          DDOTG=DDOTG-SP
          IF (DDOTG .LT. ZERO .AND. MSAT .LT. MTOT) THEN
C
C     Seek the next constraint boundary along the search direction.
C
              TEMP=ZERO
              KL=MDEG+1
              DO 60 K=KL,MTOT
              J=IACT(K)
              IF (BRES(J) .GT. ZERO) THEN
                  IF (TEMP .EQ. ZERO .OR. BRES(J) .LT. TEMP) THEN
                      TEMP=BRES(J)
                      INDXBD=K
                  END IF
              END IF
   60         CONTINUE
              IF (TEMP .GT. ZERO) THEN
                  STEPCB=TEMP
                  GOTO 40
              END IF
          END IF
      RETURN
      END
      SUBROUTINE VE11MD (N,X,NACT,G,Z,ZTG,XS,GS,ZZNORM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),G(*),Z(*),ZTG(*),XS(*),GS(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/
C
C     Test if there is sufficient convexity for the update.
C
      DD=ZERO
      DG=ZERO
      TEMP=ZERO
      DO 10 I=1,N
      XS(I)=X(I)-XS(I)
      DD=DD+XS(I)**2
      TEMP=TEMP+GS(I)*XS(I)
      GS(I)=G(I)-GS(I)
      DG=DG+GS(I)*XS(I)
   10 CONTINUE
      IF (DG .LT. 0.1*DABS(TEMP)) GOTO 90
C
C     Transform the Z matrix.
C
      K=N
   20 KP=K
      K=K-1
      IF (K .GT. NACT) THEN
          IF (ZTG(KP) .EQ. ZERO) GOTO 20
          TEMP=DABS(ZTG(KP))*DSQRT(ONE+(ZTG(K)/ZTG(KP))**2)
          WCOS=ZTG(K)/TEMP
          WSIN=ZTG(KP)/TEMP
          ZTG(K)=TEMP
          IZ=K
          DO 30 I=1,N
          TEMP=WCOS*Z(IZ+1)-WSIN*Z(IZ)
          Z(IZ)=WCOS*Z(IZ)+WSIN*Z(IZ+1)
          Z(IZ+1)=TEMP
          IZ=IZ+N
   30     CONTINUE
          GOTO 20
      END IF
C
C     Update the value of ZZNORM.
C
      IF (ZZNORM .LT. ZERO) THEN
          ZZNORM=DD/DG
      ELSE
          TEMP=DSQRT(ZZNORM*DD/DG)
          ZZNORM=DMIN1(ZZNORM,TEMP)
          ZZNORM=DMAX1(ZZNORM,0.1*TEMP)
      END IF
C
C     Complete the updating of Z.
C
      NP=NACT+1
      TEMP=DSQRT(DG)
      IZ=NP
      DO 40 I=1,N
      Z(IZ)=XS(I)/TEMP
      IZ=IZ+N
   40 CONTINUE
      IF (NP .LT. N) THEN
          KM=NP+1
          DO 80 K=KM,N
          TEMP=ZERO
          IZ=K
          DO 50 I=1,N
          TEMP=TEMP+GS(I)*Z(IZ)
          IZ=IZ+N
   50     CONTINUE
          TEMP=TEMP/DG
          SUM=ZERO
          IZ=K
          DO 60 I=1,N
          Z(IZ)=Z(IZ)-TEMP*XS(I)
          SUM=SUM+Z(IZ)**2
          IZ=IZ+N
   60     CONTINUE
          IF (SUM .LT. ZZNORM) THEN
              TEMP=DSQRT(ZZNORM/SUM)
              IZ=K
              DO 70 I=1,N
              Z(IZ)=TEMP*Z(IZ)
              IZ=IZ+N
   70         CONTINUE
          END IF
   80     CONTINUE
      END IF
   90 RETURN
      END
      BLOCK DATA VE11TD
      COMMON /VE11UD/LP
      DATA LP/6/
      END
