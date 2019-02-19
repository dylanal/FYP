

      Program CGL
      
      use variables_module   
      
      Real*8 :: Wkkr,Wkki,K0r,K0i,Gammar,Gammai
      Complex*16 :: zi
      integer NXI,NO
      
      zi=cmplx(0,1.d0)
      
!     Read the input file.
      open(1,file='input.dat')
      
      read(1,*) 
      read(1,*)
      read(1,*)
         
      read(1,*) NXI
      IF (NXI .NE. NX) THEN
      	PRINT*,'ERROR: Dimension does not match!'
      	PRINT *,'NX=', NX, 'NXI=', NXI
      	stop
      END IF
      read(1,*)
      read(1,*) Xs,Xf
      
      read(1,*) 
      read(1,*) 
      
      read(1,*) W0IMAX,W0XXI,W0r
      
      read(1,*) 
      read(1,*) K0R,K0I,WKKR,WKKI,GammaR,GammaI
      
      Print*,'------------     CGL solver is started  -----------------'
      print*,'Wkkr=',Wkkr,'Wkki=',Wkki
      print*,'K0r=',K0r,'K0i=',K0i
      print*,'GammaR=',Gammar,'GammaI=',Gammai
      print*,'---------------------------------------------------------'
      
      Wkk=wkkr+wkki*zi
      k0=k0r+k0i*zi
      Gamma=gammar+zi*gammai
      
      read(1,*)
      read(1,*)
      read(1,*) dt,CFL
      PRINT *,'dt=', dt, 'CFL=', CFL
      
      read(1,*)
      read(1,*) NT,NT_fld,NT_trace
     
      
      read(1,*)
      read(1,*)
      
      read(1,*) NO
      read(1,*)
      read(1,*) infilename
      read(1,*) intimename

      close(1)
      
!     Set grid information

      call grid     
      PRINT *,'Recommended dt=',CFL*DX/ABS(K0*WKK)    
      PRINT *,'Total time=',float(NT)*dt
         
      call infield(NO)
      
      open(2,file='trace.dat')


      do IT=1,NT       
       	time=time+dt
      	call RK3CN2
       if (mod(it,NT_trace) .eq. 0) call trace
       if (mod(it,NT_fld) .eq. 0) 	call fld_plot    	
      end do
      
      close(2)	
    
      END PROGRAM CGL
      
!-----------------------------------------------------------------------

      Subroutine grid
      
      use variables_module  
      
      COMPLEX*16 ZI 
      INTEGER I
       
      ZI=CMPLX(0.,1.D0)
       
       DO I=1,NX
       	X(I)=XS+FLOAT(I-1)/FLOAT(NX-1)*(XF-XS)
       	W0(I)=ZI*(W0IMAX-0.5D0*W0XXI*X(I)**2)+W0r
       END DO
       
       DX=ABS(X(2)-X(1))
       
     
       
      RETURN
      END SUBROUTINE GRID
      
!-----------------------------------------------------------------------

      Subroutine infield(NO)
      
      use variables_module  
       
       INTEGER I
       INTEGER NO,Ntmp
       REAL*8 :: PSIR(NX),PSII(NX)
       
       IF (NO .EQ. 0) THEN
       	TIME=0.D0
        DO I=1,NX
       	 PSI(I)=CMPLX(0.1D0,0.D0)
        END DO
       ELSE
        OPEN(2,FILE=infieldname)
        DO I=1,NX
        	READ(2,101) X(I),PSIR(I),PSII(I)
        	PSI(I)=CMPLX(PSIR(I),PSII(I))
        END DO
        CLOSE(2)
        OPEN(3,FILE=intimename)
        read(3,*) Ntmp,time
        close(3)
       end if        
       
101   format(100(1x,e16.10))        
      RETURN
      END SUBROUTINE infield      

!-----------------------------------------------------------------------

      Subroutine RK3CN2
      
      use variables_module  
      
      COMPLEX*16 RHS(NX),RXL(NX),RX(NX),RXP(NX)
      COMPLEX*16 LHSA(NX),LHSB(NX),LHSC(NX)
      complex*16 delta(nx)
      
      REAL*8 ALP(3),GAM(3),RHO(3)
      COMPLEX*16 ZI
      INTEGER IK
      
      integer i
        
      PI=acos(-1.)
      ZI=CMPLX(0.D0,1.D0)
      
      ALP(1)=4./15.
      ALP(2)=1./15.
      ALP(3)=1./6.
      
      GAM(1)=8./15.
      GAM(2)=5./12.
      GAM(3)=3./4.
      
      RHO(1)=0.d0
      RHO(2)=-17./60.
      RHO(3)=-5./12.
      
      DO i=1,NX
      	 RXL(I)=0.D0
         RX(I)=0.D0
         RXP(I)=0.D0
      END DO   
      
      DO IK=1,3
      	DO i=1,NX
         LHSA(i)=0.d0
         LHSB(i)=0.d0
         LHSC(i)=0.d0 
         RHS(I)=0.D0  
         delta(i)=0.d0 
         !if (i .eq. nx/4+10) delta(i)=0.01d0     
        END DO    
         
        DO I=2,NX-1
         RXL(I)=
     &  -ZI*(W0(I)+0.5D0*WKK*K0*K0)*PSI(I)
     &  +WKK*K0*(PSI(I+1)-PSI(I-1))/(2.D0*DX)
     &  +ZI*WKK/2.D0*(PSI(I+1)-2*PSI(I)+PSI(I-1))/(DX*DX)  
      
         RX(I)= 	 
     &  -ZI*GAMMA*(ABS(PSI(I))**2)*PSI(I)   	 
         	 
         RHS(I)=PSI(I)+DT*(ALP(IK)*RXL(I)+GAM(IK)*RX(I)+RHO(IK)*RXP(I))

         RXP(I)=RX(I)
        
         LHSA(I)=-1.*DT*ALP(IK)*(-1.D0*WKK*K0/(2.D0*DX)
     &     +ZI*0.5D0*WKK/(DX*DX))        
      	 LHSB(I)=1.-DT*ALP(IK)*(-ZI*(W0(I)+0.5D0*WKK*K0*K0)
     &     +(-2.D0*ZI*0.5D0*WKK/(DX*DX)))      	
      	 LHSC(I)=-1.*DT*ALP(IK)*(1.D0*WKK*K0/(2.D0*DX)
     &     +ZI*0.5D0*WKK/(DX*DX))
        END DO
            
      ! Boundary condition for U
      LHSB(1)=1.D0
      LHSB(NX)=1.D0
      	
      RHS(1)=0.D0
      RHS(NX)=0.D0
      
      ! SOLVE LHS PSI = RHS
       CALL TDMA_COMPLEX(LHSA,LHSB,LHSC,RHS,NX)
      
      ! UPDATE PSI
      DO I=1,NX
      	PSI(I)=RHS(I)
      END DO 	
      
      END DO  
              
      end subroutine RK3CN2		
      
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|    
      SUBROUTINE TDMA_COMPLEX(A,B,C,G,NY)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Uses the Thomas algorithm to solve Ax=b for tridiagonal A
! The RHS vector and solution is complex
! Input lower, main, and upper diagonals, ld, md, ud, and rhs x
! Returns solution in x
! The indexing should be done by ROW, ie.
! [ b1  c1   0   0   0 ...
! [ a2  b2  c2   0   0 ...
! [  0  a3  b3   c3  0 ...

      INTEGER I, J, NY
      COMPLEX*16 A(NY), B(NY), C(NY)
      COMPLEX*16 G(NY)

      DO J=1,NY-1
          A(J+1)=-A(J+1)/B(J)
          B(J+1)=B(J+1)+A(J+1)*C(J)
          G(J+1)=G(J+1)+A(J+1)*G(J)
      END DO
       
       G(NY)=G(NY)/B(NY)
      
        DO J=NY-1,1,-1
          G(J)=(G(J)-C(J)*G(J+1))/B(J)
        END DO

      RETURN
      END SUBROUTINE TdMA_COMPLEX
      

!-----------------------------------------------------------------------        

      SUBROUTINE FLD_PLOT
      
      use variables_module  
      
      INTEGER I
      real*8 pi
      
      character*16 tname,tname2
      character*3 tfn1,tfn2
      integer idg1,idg2,idg3,idg4,nn
      
      pi=acos(-1.)
      
      tfn1='fld.dat'
      tfn2='time.tec'
      nn=it/100
      idg1=nn/1000
      idg2=(nn-idg1*1000)/100
      idg3=(nn-idg1*1000-idg2*100)/10
      idg4=nn-idg1*1000-idg2*100-idg3*10
      tname=tfn1//char(idg1+48)//char(idg2+48)//
     &      char(idg3+48)//char(idg4+48)//'.tec'
      tname2=tfn2//char(idg1+48)//char(idg2+48)//
     &      char(idg3+48)//char(idg4+48)//'.tec'
      
      
      open(317,file=tname, status='unknown',form='formatted')
        do i=1,nx
        	write(317,200) x(i),dble(psi(i)),dimag(psi(i)),abs(psi(i)) 	
        ENDDO
200   format(100(1x,e16.10))    
      close(317)	    
      
      open(318,file=tname2, status='unknown',form='formatted')
        	write(318,*) it,time  
      close(318)	 
          
      return
      end subroutine FLD_PLOT  

!-----------------------------------------------------------------------

      Subroutine trace
      
      use variables_module  
      real*8 energy
      
      energy=0.d0
      
      do i=2,Nx
      	energy=energy+0.5*(abs(psi(i))+abs(psi(i-1)))*dx
      end do
      
      write(2,100) time,energy,real(psi(nx/2))

100   format(100(1x,e16.10))   

      return
      end subroutine trace 