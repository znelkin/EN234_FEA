!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also contains the following subrouines:
!          abq_UEL_2D_integrationpoints           - defines integration points for 2D continuum elements
!          abq_UEL_2D_shapefunctions              - defines shape functions for 2D continuum elements
!          abq_UEL_1D_integrationpoints(n_points, n_nodes, xi, w)  = defines integration points for 1D line integral
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UELch(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
    !  INCLUDE 'ABA_PARAM.INC'
       implicit double precision (a-h,o-z)
    !
    !
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

    !
    !       Variables that must be computed in this routine
    !       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
    !       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
    !       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
    !       ENERGY(1:8)
    !                                  Energy(1) Kinetic Energy
    !                                  Energy(2) Elastic Strain Energy
    !                                  Energy(3) Creep Dissipation
    !                                  Energy(4) Plastic Dissipation
    !                                  Energy(5) Viscous Dissipation
    !                                  Energy(6) Artificial strain energy
    !                                  Energy(7) Electrostatic energy
    !                                  Energy(8) Incremental work done by loads applied to the element
    !       PNEWDT                     Allows user to control ABAQUS time increments.
    !                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
    !                                  a time increment equal to PNEWDT*DTIME
    !                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
    !
    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       NRHS                       Dimension variable
    !       NSVARS                     Total # element state variables
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U                          Vector of DOF at the end of the increment
    !       DU                         Vector of DOF increments
    !       V                          Vector of velocities (defined only for implicit dynamics)
    !       A                          Vector of accelerations (defined only for implicit dynamics)
    !       JTYPE                      Integer identifying element type (the number n in the Un specification in the input file)
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
    !       NDLOAD                     Number of user-defined distributed loads defined for this element
    !       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(1:NDLOAD)           Distributed load magnitudes
    !       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
    !       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
    !       PREDEF(1,...)              Value of predefined field
    !       PREDEF(2,...)              Increment in predefined field
    !       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
    !       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable
    !       LFLAGS(1)                  Defines procedure type
    !       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
    !       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
    !                                   2 => Subroutine must return stiffness AMATRX = -dF/du
    !                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
    !                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
    !                                   5 => Define the RHS only
    !                                   6 => Define the mass matrix for the initial acceleration calculation
    !                                   100 => Define perturbation quantities for output
    !       LFLAGS(4)                   0 => General step   1 => linear perturbation step
    !       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
    !       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
    !       PERIOD                      Time period of the current step
    !
    !
    ! Local Variables
      integer      :: i,j,n_points,kint, nfacenodes, ipoin, ksize
      integer      :: face_node_list(3)                       ! List of nodes on an element face
    !
      double precision  ::  xi(2,9)                          ! Area integration points
      double precision  ::  w(9)                             ! Area integration weights
      double precision  ::  N(9)                             ! 2D shape functions
      double precision  ::  Nbar(4)                          ! 2D shape functions for mu and c 
      double precision  ::  dNdxi(9,2)                       ! 2D shape function derivatives
      double precision  ::  dNdxibar(4,2)                    ! 2D shape function derivatives for mu and c 
      double precision  ::  dNdx(9,2)                        ! Spatial derivatives
      double precision  ::  dNdxbar(4,2)                     ! Spatial derivatives for mu and c
      double precision  ::  dxdxi(2,2)                       ! Derivative of spatial coords wrt normalized coords
      double precision  ::  dxdxibar(2,2)                    ! Derivative of spatial coords wrt normalized coords

    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(2,3)                  ! Coords of nodes on an element face
      double precision  ::  xi1(6)                            ! 1D integration points
      double precision  ::  w1(6)                              ! Integration weights
      double precision  ::  N1(3)                             ! 1D shape functions
      double precision  ::  dN1dxi(3)                         ! 1D shape function derivatives
      double precision  ::  norm(2)                           ! Normal to an element face
      double precision  ::  dxdxi1(2)                         ! Derivative of 1D spatial coord wrt normalized areal coord
    !
      double precision  ::  strain(4)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
      double precision  ::  stress(4)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  Del(4,4)                          ! stress = D*(strain)  (NOTE FACTOR OF 2 in shear strain)
      double precision  ::  B(9,24)                           ! strain = B*(dof_total)
      double precision  ::  ktemp(22,22)                      ! Temporary stiffness (for incompatible mode elements)
      double precision  ::  rhs_temp(22)                      ! Temporary RHS vector (for incompatible mode elements)
      double precision  ::  dxidx(2,2), dxidxbar(2,2)         ! Jacobian inverse
      double precision  ::  determinant, determinantbar       ! determinant
      double precision  ::  E, xnu, D44, D11, D12             ! Material properties
      double precision  ::  Omega, WGibbs, Kappa, Diffc, Theta! Diffusion Properties
      double precision  ::  F, dFdc                           ! Gibbs Free Energy & its derivative
      double precision  ::  D(9,9)                            ! New D matrix
      double precision  ::  SumD
      double precision  ::  qvec(9)                              ! qvec

    !
    !     Example ABAQUS UEL implementing 2D linear elastic elements
    !     Includes option for incompatible mode elements
    !     El props are:

    !     PROPS(1)         Young's modulus
    !     PROPS(2)         Poisson's ratio


      if (NNODE == 3) n_points = 1              ! Linear triangle
      if (NNODE == 4) n_points = 4               ! Linear rectangle
      if (NNODE == 6) n_points = 4              ! Quadratic triangle
      if (NNODE == 8) n_points = 9               ! Serendipity rectangle
      if (NNODE == 9) n_points = 9             ! Quadratic rect

    ! Write your code for a 2D element below

      call abq_UEL_2D_integrationpoints(n_points, NNODE, xi, w)
      
      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      

      
      
      ktemp(1:2*NNODE+4,1:2*NNODE+4) = 0.d0
      rhs_temp(1:2*NNODE+4) = 0.d0
      Del = 0.d0
      E = PROPS(1)
      xnu = PROPS(2)
      !New Diffusion-Related Parameters
      Omega = PROPS(3)
      WGibbs= PROPS(4)
      Kappa = PROPS(5)
      Diffc = PROPS(6)
      Theta = PROPS(7)
      !Define elasic D matrix outside the integration loop so that the Complete matrix can be assembled inside 
      D12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
      D11 = (1-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
      D44 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
      Del(1:3,1:2) = D12
      Del(1,1) = D11
      Del(2,2) = D11
      Del(3,3) = D11
      Del(4,4) = D44
      ! New D matrix
      D = 0.d0
      
      Energy(1:8) = 0.d0
      
      call abq_UEL_2D_shapefunctions([0.d0,0.d0],NNODE,N,dNdxi)
          dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
          det0 = dxdxi(1,1)*dxdxi(2,2) - dxdxi(2,1)*dxdxi(1,2)

  
          
          
      do kint = 1, 4! n_points
          call abq_UEL_2D_shapefunctions(xi(1:2,kint),8,N,dNdxi)
          dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
          determinant = dxdxi(1,1)*dxdxi(2,2) - dxdxi(2,1)*dxdxi(1,2)
          dxidx(1,1) = dxdxi(2,2)
          dxidx(2,2) = dxdxi(1,1)
          dxidx(1,2) = -dxdxi(1,2)
          dxidx(2,1) = -dxdxi(2,1)
          dxidx = dxidx/determinant

          dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)
          
          ! Get the Interpolation functions for the four noded element
          call abq_UEL_2D_shapefunctions(xi(1:2,kint),4,Nbar,dNdxibar)
          dxdxibar(1:2,1:2) = matmul(coords(1:2,1:4),dNdxibar(1:4,1:2))
          determinantbar = dxdxibar(1,1)*dxdxibar(2,2) - 
     1    dxdxibar(2,1)*dxdxibar(1,2)
          dxidxbar(1,1) = dxdxibar(2,2)
          dxidxbar(2,2) = dxdxibar(1,1)
          dxidxbar(1,2) = -dxdxibar(1,2)
          dxidxbar(2,1) = -dxdxibar(2,1)
          dxidxbar = dxidxbar/determinantbar

          dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)
          
          ! Get dNdxbar
          
          ! Define the modified B matrix
          B = 0.d0
          ! Elasticity related parts of the B Matrix
          B(1,1)  = dNdx(1,1)
!          B(1,2)  =
!          B(1,3)  =
!          B(1,4)  =
          B(1,5)  = dNdx(2,1)
!          B(1,6)  =
!          B(1,7)  =
!          B(1,8)  =
          B(1,9)  = dNdx(3,1)
!          B(1,10) = 
!          B(1,11) =
!          B(1,12) =
          B(1,13) = dNdx(4,1)
!          B(1,14) =
!          B(1,15) =
!          B(1,16) =
          B(1,17) = dNdx(5,1)
!          B(1,18) =
          B(1,19) = dNdx(6,1)
!          B(1,20) =
          B(1,21) = dNdx(7,1)
!          B(1,22) =
          B(1,23) = dNdx(8,1)
!          B(1,24) =
          
!          B(2,1)  =
          B(2,2)  = dNdx(1,2)
!          B(2,3)  =
!          B(2,4)  =
!          B(2,5)  =
          B(2,6)  = dNdx(2,2)
!          B(2,7)  =
!          B(2,8)  =
!          B(2,9)  =
          B(2,10) = dNdx(3,2) 
!          B(2,11) =
!          B(2,12) =
!          B(2,13) =
          B(2,14) = dNdx(4,2)
!          B(2,15) =
!          B(2,16) =
!          B(2,17) =
          B(2,18) = dNdx(5,2)
!          B(2,19) =
          B(2,20) = dNdx(6,2)
!          B(2,21) =
          B(2,22) = dNdx(7,2)
!          B(2,23) =
          B(2,24) = dNdx(8,2)
          
          B(3,1)  = dNdx(1,2)
          B(3,2)  = dNdx(1,1)
!          B(3,3)  =
!          B(3,4)  =
          B(3,5)  = dNdx(2,2)
          B(3,6)  = dNdx(2,1)
!          B(3,7)  =
!          B(3,8)  =
          B(3,9)  = dNdx(3,2)
          B(3,10) = dNdx(3,1)
!          B(3,11) =
!          B(3,12) =
          B(3,13) = dNdx(4,2)
          B(3,14) = dNdx(4,1)
!          B(3,15) =
!          B(3,16) =
          B(3,17) = dNdx(5,2)
          B(3,18) = dNdx(5,1)
          B(3,19) = dNdx(6,2)
          B(3,20) = dNdx(6,1)
          B(3,21) = dNdx(7,2)
          B(3,22) = dNdx(7,1)
          B(3,23) = dNdx(8,2)
          B(3,24) = dNdx(8,1)
          ! Mu and C related parts of the B Matrix
!          B(4,1)  = 
!          B(4,2)  =
          B(4,3)  = Nbar(1)
!          B(4,4)  =
!          B(4,5)  = 
!          B(4,6)  =
          B(4,7)  = Nbar(2)
!          B(4,8)  =
!          B(4,9)  = 
!          B(4,10) = 
          B(4,11) = Nbar(3)
!          B(4,12) =
!          B(4,13) = 
!          B(4,14) =
          B(4,15) = Nbar(4)
!          B(4,16) =
!          B(4,17) = 
!          B(4,18) =
!          B(4,19) = 
!          B(4,20) =
!          B(4,21) = 
!          B(4,22) =
!          B(4,23) = 
!          B(4,24) =
          
!          B(5,1)  = 
!          B(5,2)  =
!          B(5,3)  =
          B(5,4)  = Nbar(1)
!          B(5,5)  = 
!          B(5,6)  =
!          B(5,7)  =
          B(5,8)  = Nbar(2)
!          B(5,9)  = 
!          B(5,10) = 
!          B(5,11) =
          B(5,12) = Nbar(3)
!          B(5,13) = 
!          B(5,14) =
!          B(5,15) =
          B(5,16) = Nbar(4)
!          B(5,17) = 
!          B(5,18) =
!          B(5,19) = 
!          B(5,20) =
!          B(5,21) = 
!          B(5,22) =
!          B(5,23) = 
!          B(5,24) =

          
          ! Parts of the B Matrix related to the Derivatives of Mu and C
!          B(6,1)  = 
!          B(6,2)  =
          B(6,3)  = dNdxbar(1,1)
!          B(6,4)  =
!          B(6,5)  = 
!          B(6,6)  =
          B(6,7)  = dNdxbar(2,1)
!          B(6,8)  =
!          B(6,9)  = 
!          B(6,10) = 
          B(6,11) = dNdxbar(3,1)
!          B(6,12) =
!          B(6,13) = 
!          B(6,14) =
          B(6,15) = dNdxbar(4,1)
!          B(6,16) =
!          B(6,17) = 
!          B(6,18) =
!          B(6,19) = 
!          B(6,20) =
!          B(6,21) = 
!          B(6,22) =
!          B(6,23) = 
!          B(6,24) =
          
!          B(7,1)  = 
!          B(7,2)  =
          B(7,3)  = dNdxbar(1,2)
!          B(7,4)  =
!          B(7,5)  = 
!          B(7,6)  =
          B(7,7)  = dNdxbar(2,2)
!          B(7,8)  =
!          B(7,9)  = 
!          B(7,10) = 
          B(7,11) = dNdxbar(3,2)
!          B(7,12) =
!          B(7,13) = 
!          B(7,14) =
          B(7,15) = dNdxbar(4,2)
!          B(7,16) =
!          B(7,17) = 
!          B(7,18) =
!          B(7,19) = 
!          B(7,20) =
!          B(7,21) = 
!          B(7,22) =
!          B(7,23) = 
!          B(7,24) =
          
!          B(8,1)  = 
!          B(8,2)  =
!          B(8,3)  =
          B(8,4)  = dNdxbar(1,1)
!          B(8,5)  = 
!          B(8,6)  =
!          B(8,7)  =
          B(8,8)  = dNdxbar(2,1)
!          B(8,9)  = 
!          B(8,10) = 
!          B(8,11) =
          B(8,12) = dNdxbar(3,1)
!          B(8,13) = 
!          B(8,14) =
!          B(8,15) =
          B(8,16) = dNdxbar(4,1)
!          B(8,17) = 
!          B(8,18) =
!          B(8,19) = 
!          B(8,20) =
!          B(8,21) = 
!          B(8,22) =
!          B(8,23) = 
!          B(8,24) =
          
!          B(9,1)  = 
!          B(9,2)  =
!          B(9,3)  =
          B(9,4)  = dNdxbar(1,2)
!          B(9,5)  = 
!          B(9,6)  =
!          B(9,7)  =
          B(9,8)  = dNdxbar(2,2)
!          B(9,9)  = 
!          B(9,10) = 
!          B(9,11) =
          B(9,12) = dNdxbar(3,2)
!          B(9,13) = 
!          B(9,14) =
!          B(9,15) =
          B(9,16) = dNdxbar(4,2)
!          B(9,17) = 
!          B(9,18) =
!          B(9,19) = 
!          B(9,20) =
!          B(9,21) = 
!          B(9,22) =
!          B(9,23) = 
!          B(9,24) =
          
          ! Define the Gibbs Free energy as a function of concentration
          F = 2.d0*WGibbs*U(5)
     1     *(U(5) - 1.d0)
     2      *(2.d0*U(5) - 1.d0)
          
          ! Define the derivative of the Gibbs Free Energy
          dFdc = WGibbs*(12.d0*(U(5) + DU(5,1))**2.d0
     1     - 12.d0*(U(5) + DU(5,1)) + 2.d0)
          
          ! Define the summation for D(4,5)
          SumD = (1.d0/3.d0)*(D11 + D22 + D33 + 6.d0*D12)
          !Define the new D matrix
          D(1:3,1:3) = Del(1:3,1:3)
          D(1,5) = (-1.d0/3.d0)*Omega*(D11 + 2.d0*D12)
          D(2,5) = (-1.d0/3.d0)*Omega*(D22 + 2.d0*D12)
          D(4,1) = (-1.d0/3.d0)*Omega*(D11 + 2.d0*D12)
          D(4,2) = (-1.d0/3.d0)*Omega*(D22 + 2.d0*D12)
          D(4,4) = 1.d0
          D(4,5) = -dFdc + (Omega**2.d0)*SumD
          D(6,8) = -Kappa
          D(7,9) = -Kappa
          D(8,6) = Theta*Diffc
          D(8,7) = Theta*Diffc
      end do
        
      do kint = 1, n_points
          call abq_UEL_2D_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)
          dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
          determinant = dxdxi(1,1)*dxdxi(2,2) - dxdxi(2,1)*dxdxi(1,2)
          dxidx(1,1) = dxdxi(2,2)
          dxidx(2,2) = dxdxi(1,1)
          dxidx(1,2) = -dxdxi(1,2)
          dxidx(2,1) = -dxdxi(2,1)
          dxidx = dxidx/determinant

          dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)
          B = 0.d0
          B(1,1:2*NNODE-1:2) = dNdx(1:NNODE,1)
          B(1,2*NNODE+1) = (det0/determinant)*xi(1,kint)*dxidx(1,1)
          B(1,2*NNODE+3) = (det0/determinant)*xi(2,kint)*dxidx(2,1)
          B(2,2:2*NNODE:2) = dNdx(1:NNODE, 2)
          B(2,2*NNODE+2) = (det0/determinant)*xi(1,kint)*dxidx(1,2)
          B(2,2*NNODE+4) = (det0/determinant)*xi(2,kint)*dxidx(2,2)
          B(4,1:2*NNODE-1:2) = dNdx(1:NNODE, 2)
          B(4,2:2*NNODE:2) = dNdx(1:NNODE, 1)
          B(4,2*NNODE+1) = (det0/determinant)*xi(1,kint)*dxidx(1,2)
          B(4,2*NNODE+2) = (det0/determinant)*xi(1,kint)*dxidx(1,1)
          B(4,2*NNODE+3) = (det0/determinant)*xi(2,kint)*dxidx(2,2)
          B(4,2*NNODE+4) = (det0/determinant)*xi(2,kint)*dxidx(2,1)
          
          ! Calculate the stress and strain 
          strain(1:3) = matmul(B(1:3,1:24),U(1:24))
          stress(1:3) = matmul(Del(1:3,1:3),strain(1,3)
          
          ! Calculae vector q
          qvec = 0.d0
          qvec(1) = stress(1)
          qvec(2) = stress(2)
          qvec(3) = stress(3)
          
          qvec(4) = U(4) + DU(4,1)
     1     -(2.d0*WGibbs*(U(5) + DU(5,1))
     2      *(U(5,1) + DU(5,1) -1.d0)*(2.d0*(U(5,1) + DU(5,1))-1.d0))
     3       -Omega*(stress(1)+stress(2))
          qvec(5) = DU(5,1)/DTIME
          
          qvec(6) = -Kappa*(U(8) + DU(8,1))
          qvec(7) = -Kappa*(U(9) + DU(9,1))
          
          qvec(8) = Diffc*(U(6) + Theta*DU(6,1))
          qvec(9) = Diffc*(U(7) + Theta*DU(7,1))
!          rhs_temp(1:2*NNODE+4) = rhs_temp(1:2*NNODE+4) -
!     1     matmul(transpose(B(1:4,1:2*NNODE+4)),stress)*w(kint)*
!     2      determinant
!          SVARS(4*kint-3: 4*kint) = stress(1:4)
!          Energy(2) = Energy(2) + 0.5*dot_product(stress,strain)*w(kint)
!     1       *determinant
      end do
      
!          AMATRX(1:2*NNODE,1:2*NNODE) = kuu(1:2*NNODE,1:2*NNODE)
!     1     -matmul(kua(1:2*NNODE,1:4),
!     2     matmul(kaainv(1:4,1:4),kau(1:4,1:2*NNODE)))
!          RHS(1:2*NNODE,1) = rhs_temp(1:2*NNODE) -
!     1     matmul(kua(1:2*NNODE,1:4),matmul(kaainv(1:4,1:4),rhs_temp
!     2     (2*NNODE+1:2*NNODE+4)))
          
          
      PNEWT = 1.d0
      
      return
      

      END SUBROUTINE UELch