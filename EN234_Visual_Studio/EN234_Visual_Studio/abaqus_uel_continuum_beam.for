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

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
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
      double precision  ::  dNdxi(9,2)                       ! 2D shape function derivatives
      double precision  ::  dNdx(9,2)                        ! Spatial derivatives
      double precision  ::  dxdxi(2,2)                       ! Derivative of spatial coords wrt normalized coords

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
      double precision  ::  D(4,4)                            ! stress = D*(strain)  (NOTE FACTOR OF 2 in shear strain)
      double precision  ::  B(4,22)                           ! strain = B*(dof_total)
      double precision  ::  ktemp(22,22)                      ! Temporary stiffness (for incompatible mode elements)
      double precision  ::  rhs_temp(22)                      ! Temporary RHS vector (for incompatible mode elements)
      double precision  ::  kuu(18,18)                        ! Upper diagonal stiffness
      double precision  ::  kaa(4,4),kaainv(4,4)              ! Lower diagonal stiffness
      double precision  ::  kau(4,18)                         ! Lower quadrant of stiffness
      double precision  ::  kua(18,4)                         ! Upper quadrant of stiffness
      double precision  ::  alpha(4)                          ! Internal DOF for incompatible mode element
      double precision  ::  dxidx(2,2), determinant, det0     ! Jacobian inverse and determinant
      double precision  ::  E, xnu, D44, D11, D12             ! Material properties
      double precision  ::  h, width                              ! Cross Section of beam
      double precision  ::  L, delta_x, delta_y
      double precision  ::  Trig_Theta(2,2)                   ! row is the ith theta, column is sine cosine
      double precision  ::  s_coords(4, 2)                    ! s node positions
      double precision  ::  T(8,6)
      double precision  ::  wt(20)                            ! Integration weights for trapazoidal integration
      double precision  ::  e_hat_One(2), e_hat_Two(2)        ! Laminar Basis Vectors
      double precision  ::  M_One, M_Two                      ! Magnitude for the dxdxi vectors
      double precision  ::  R(2,4)                            ! Rotation Matrix
      double precision  ::  epsilon_three(3)                  ! Temporary variable for strain vector
      double precision  ::  Strain_hat(2), Stress_hat(2)      ! Stress and Strain in new coordinate system

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
      
      RHS(1:6,1) = 0.d0
      AMATRX(1:6,1:6) = 0.d0
      SVARS(1:3) = 0.d0
      
      D = 0.d0
      E = PROPS(3)
      xnu = PROPS(4)
      h = PROPS(1)
      w = PROPS(2)
      
      ! Calculate L
      delta_x = coords(1,2) - coords(1,1)
      delta_y = coords(2,2) - coords(2,1)
      L = sqrt(delta_x**2.d0 + delta_y**2.d0)
      ! Calculate trigonometric quantities for s_nodes
      Trig_Theta = 0.d0
      ! Sine
      Trig_Theta(1,1) = (coords(2,2) - coords(2,1))/L
      Trig_Theta(2,1) = (coords(2,2) - coords(2,1))/L
      ! Cosine
      Trig_Theta(1,2) = (coords(1,2) - coords(1,1))/L 
      Trig_Theta(2,2) = (coords(1,2) - coords(1,1))/L
      
      ! Calculate the coordinates of snodes
      s_coords = 0.d0 
      
      s_coords(1,1) = coords(1,1) + (h/2.d0)*Trig_Theta(1,1)
      s_coords(2,1) = coords(1,2) + (h/2.d0)*Trig_Theta(2,1)
      s_coords(3,1) = coords(1,2) - (h/2.d0)*Trig_Theta(2,1)
      s_coords(4,1) = coords(1,1) - (h/2.d0)*Trig_Theta(1,1)
      
      s_coords(1,2) = coords(2,1) - (h/2.d0)*Trig_Theta(1,2)
      s_coords(2,2) = coords(2,2) - (h/2.d0)*Trig_Theta(2,2)
      s_coords(3,2) = coords(2,2) + (h/2.d0)*Trig_Theta(2,2)
      s_coords(4,2) = coords(2,1) + (h/2.d0)*Trig_Theta(1,2)
      
      ! Calculate T to map displacements
      T = 0.d0
      
      T(1,1) = 1.d0
      T(1,3) = coords(2,1) - s_coords(1,2)
      
      T(2,2) = 1.d0
      T(2,3) = -(coords(1,1) - s_coords(1,1))
      
      T(3,4) = 1.d0
      T(3,6) = coords(2,2) - s_coords(2,2)
      
      T(4,5) = 1.d0
      T(4,6) = coords(1,2) - s_coords(2,1)
      
      T(5,4) = 1.d0
      T(5,6) = coords(2,2) - s_coords(3,2)
      
      T(6,5) = 1.d0
      T(6,6) = coords(1,2) - s_coords(3,1)
      
      T(7,1) = 1.d0
      T(7,3) = coords(2,1) - s_coords(4,2)
      
      T(8,2) = 1.d0
      T(8,3) = -(coords(1,1) - s_coords(4,1))
      
      ! Create D matrix
      D = 0.d0
      D(1,1) = E
      D(2,2) = E/(2.d0*(1.d0+xnu))
      !SVARS(1:3) = 0.d0

      call Trapazoidal_Integration_Points(wt)      
      Energy(1:8) = 0.d0
      
      do kint = 1, 20

          ! Get Laminar Basis Vectors
          call abq_UEL_2D_shapefunctions([0,xitrap(kint)],4,N,dNdxi)
          dxdxi = matmul(transpose(s_coords(1:4,1:2)),dNdxi(1:4,1:2))
          determinant = dxdxi(1,1)*dxdxi(2,2) - dxdxi(2,1)*dxdxi(1,2)
          dxidx(1,1) = dxdxi(2,2)
          dxidx(2,2) = dxdxi(1,1)
          dxidx(1,2) = -dxdxi(1,2)
          dxidx(2,1) = -dxdxi(2,1)
          dxidx = dxidx/determinant
          
          e_hat_One(1) = dxdxi(1,1)
          e_hat_One(2) = dxdxi(1,2)
          M_One = sqrt(dxdxi(1,1)**2.d0 + dxdxi(1,2)**2.d0)
          e_hat_One = e_hat_One/M_One
          
          e_hat_Two(1) = dxdxi(2,1)
          e_hat_Two(2) = dxdxi(2,2)
          M_Two = sqrt(dxdxi(2,1)**2.d0 + dxdxi(2,2)**2.d0)
          e_hat_Two = e_hat_Two/M_Two
          
          ! Get R and B matrices
          B = 0.d0
          B(1,1:2*4-1:2) = dNdx(1:4,1)
          B(2,2:2*4:2) = dNdx(1:4, 2)
          B(3,1:2*4-1:2) = dNdx(1:4, 2)
          B(3,2:2*4:2) = dNdx(1:4, 1)
          
          R = 0.d0
          R(1,1) = e_hat_One(1)**2.d0
          R(1,2) = e_hat_One(2)**2.d0
          R(1,3) = e_hat_One(1)*e_hat_One(2)
          R(1,4) = R(1,3)
          R(2,1) = 2.d0*e_hat_One(1)*e_hat_Two(1)
          R(2,2) = 2.d0*e_hat_One(2)*e_hat_Two(2)
          R(2,3) = e_hat_One(1)*e_hat_Two(2) + e_hat_One(2)*e_hat_Two(1)
          R(2,4) = R(2,3)
          
          ! Calculate Strain
          strain = 0.d0
          strain(1:3) = matmul(B(1:3,1:8),matmul(T(1:8,1:6),U(1:6)))
          epsilon_three(1) = strain(1)
          epsilon_three(2) = strain(2)
          epsilon_three(3) = strain(3)
          
          ! Calculate Strain_hat
          Strain_hat = 0.d0
          Strain_hat(1:2) = matmul(R(1:2,1:3),epsilon_three(1:3))
          
          ! Calculate Stress_hat
          Stress_hat = 0.d0
          Stress_hat(1:2) = matmul(D(1:2,1:2),Strain_hat(1:2))
          
      
          RHS(1:6, 1) = RHS(1:6,1)
     1       -width*L*h*5.d-1*matmul(transpose(matmul(R(1:2,1:3)
     2        ,matmul(B(1:3,1:8),T(1:8,1:6)))),Stress_hat(1:2))*wt(kint)

          
          AMATRX(1:6,1:6) = AMATRX(1:6,1:6) + width*L*h*5.d-1*matmul(
     1     transpose(matmul(R(1:2,1:3),matmul(B(1:3,1:8),T(1:8,1:6)))),
     2      matmul(D(1:2,1:2),matmul(R(1:2,1:3),matmul(B(1:3,1:8),
     3       T(1:8,1:6)))))*wt(kint)
          
          call psi(kint,value)
          SVARS(1) = SVARS(1) + width*h*5.d-1*Stress_hat(1)
          SVARS(2) = SVARS(2) + width*h*5.d-1*Stress_hat(2)
          SVARS(3) = SVARS(3) + width*h**2.d0*25.d-2*Stress_hat(1)*value
     1     *wt(kint)
      
          
      end do
      PNEWT = 1.d0
      
      return
      

      END SUBROUTINE UEL
    
      subroutine psi(step,value)
      
      implicit none 
      
      integer, intent(in) :: step
      
      double precision, intent(out) :: value

      value = (1.d0/10.d0)*step - 1.d0
      end subroutine psi

      subroutine Trapazoidal_Integration_Points(wt)
      
      implicit none
      
      double precision, intent(out) :: wt(*)
      
      wt(1)  = 1.d0/20.d0
      wt(2)  = 1.d0/10.d0
      wt(3)  = 1.d0/10.d0 
      wt(4)  = 1.d0/10.d0
      wt(5)  = 1.d0/10.d0
      wt(6)  = 1.d0/10.d0
      wt(7)  = 1.d0/10.d0
      wt(8)  = 1.d0/10.d0
      wt(9)  = 1.d0/10.d0
      wt(10) = 1.d0/10.d0
      wt(11) = 1.d0/10.d0
      wt(12) = 1.d0/10.d0
      wt(13) = 1.d0/10.d0
      wt(14) = 1.d0/10.d0
      wt(15) = 1.d0/10.d0
      wt(16) = 1.d0/10.d0
      wt(17) = 1.d0/10.d0
      wt(18) = 1.d0/10.d0
      wt(19) = 1.d0/10.d0
      wt(20) = 1.d0/20.d0
      
      end subroutine Trapazoidal_Integration_Points
      