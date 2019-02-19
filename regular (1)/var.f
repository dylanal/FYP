      
      
      module variables_module

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Spatial resolution parameters
! (We hardwire these into the code so that the compiler may perform
!  optimizations based on the grid size at compile time).
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INTEGER   Nx
      INCLUDE   'grid_def'
      
      COMPLEX*16 :: PSI(NX)
      COMPLEX*16 :: W0(NX)
      REAL*8 :: X(NX)
      REal*8 :: Xs,Xf,dx
      Real*8 :: W0IMAX,W0XXI,W0r
      COMPLEX*16 :: WKK,K0,Gamma
      
      REAL*8 :: DT,CFL
      
      INTEGER :: IT,NT,NT_fld,NT_trace     
      REAL*8 :: TIME
      
      character*16 infilename,intimename
      
      
      end module variables_module
      
      