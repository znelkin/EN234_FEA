!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also contains the following subrouines:
!          abq_UEL_3D_integrationpoints           - defines integration ponits for 3D continuum elements
!          abq_UEL_3D_shapefunctions              - defines shape functions for 3D continuum elements
!          abq_UEL_invert3D                       - computes the inverse and determinant of a 3x3 matrix
!          abq_facenodes_3D                       - returns list of nodes on the face of a 3D element
!
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL_3D(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
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
      integer      :: i,j,n_points,kint, nfacenodes, ipoin
      integer      :: face_node_list(8)                       ! List of nodes on an element face
    !
      double precision  ::  xi(3,64)                          ! Volumetric Integration points
      double precision  ::  w(64)                             ! Integration weights
      double precision  ::  N(20)                             ! 3D Shape functions
      double precision  ::  dNdxi(20,3)                       ! 3D Shape function derivatives
      double precision  ::  dxdxi(3,3)                        ! Derivative of position wrt normalized coords
      double precision  ::  dNdx(20,3)                        ! Derivative of shape functions wrt spatial coords
    !
    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(3,8)                  ! Coords of nodes on an element face
      double precision  ::  xi2(2,9)                          ! Area integration points
      double precision  ::  N2(9)                             ! 2D shape functions
      double precision  ::  dNdxi2(9,2)                       ! 2D shape function derivatives
      double precision  ::  norm(3)                           ! Normal to an element face
      double precision  ::  dxdxi2(3,2)                       ! Derivative of spatial coord wrt normalized areal coord
    !
      double precision  ::  strain(6)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
      double precision  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  D(6,6)                            ! stress = D*(strain)  (NOTE FACTOR OF 2 in shear strain)
!      double precision  ::  B(6,60)                           ! strain = B*(dof_total)
      double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
      double precision  ::  E, xnu, D44, D11, D12             ! Material properties
      double precision  ::  F(3,3)                            ! Deformation Gradient
      double precision  ::  Finv(3,3)                         ! Inverse of Deformation Gradient
      double precision  ::  JJ                                ! Determinant of F
      double precision  ::  B(3,3)                            ! Cauchy-Green Deformation Tensor    
    !
    !     Example ABAQUS UEL implementing 3D linear elastic elements
    !     El props are:

    !     PROPS(1)         Young's modulus
    !     PROPS(2)         Poisson's ratio

      if (NNODE == 4) n_points = 1               ! Linear tet
      if (NNODE == 10) n_points = 4              ! Quadratic tet
      if (NNODE == 8) n_points = 8               ! Linear Hex
      if (NNODE == 20) n_points = 27             ! Quadratic hex

      call abq_UEL_3D_integrationpoints(n_points, NNODE, xi, w)

      if (MLVARX<3*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 3*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0


      ENERGY(1:8) = 0.d0

    !     --  Loop over integration points
      do kint = 1, n_points
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)
        dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)
        
    !    Calculate F
         do i = 1,3
            ie = 3*(NNODE-1)+i
            F(i,1:3) = matmul(U(i:ie:3),dNdx(1:NNODE,1:3))
            F(i,i) = F(i,i) + 1.d0
         end do
    !    Calculate Cauchy-Green Tensor
         B = matmul(F,transpose(F))
    !    Calculate Finv and the Jacobian
         call abq_UEL_invert3d(F,Finv,JJ)
         
        
!        B = 0.d0
!        B(1,1:3*NNODE-2:3) = dNdx(1:NNODE,1)
!        B(2,2:3*NNODE-1:3) = dNdx(1:NNODE,2)
!        B(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)
!        B(4,1:3*NNODE-2:3) = dNdx(1:NNODE,2)
!        B(4,2:3*NNODE-1:3) = dNdx(1:NNODE,1)
!        B(5,1:3*NNODE-2:3) = dNdx(1:NNODE,3)
!        B(5,3:3*NNODE:3)   = dNdx(1:NNODE,1)
!        B(6,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
!        B(6,3:3*NNODE:3)   = dNdx(1:NNODE,2)

!        strain = matmul(B(1:6,1:3*NNODE),U(1:3*NNODE))

!        stress = matmul(D,strain)
!        RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1)
!     1   - matmul(transpose(B(1:6,1:3*NNODE)),stress(1:6))*
!     2                                          w(kint)*determinant

!        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
!     1  + matmul(transpose(B(1:6,1:3*NNODE)),matmul(D,B(1:6,1:3*NNODE)))
!     2                                             *w(kint)*determinant

!        ENERGY(2) = ENERGY(2)
!     1   + 0.5D0*dot_product(stress,strain)*w(kint)*determinant           ! Store the elastic strain energy

!        if (NSVARS>=n_points*6) then   ! Store stress at each integration point (if space was allocated to do so)
!            SVARS(6*kint-5:6*kint) = stress(1:6)
!        endif
      end do


      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
    !
    !   Apply distributed loads
    !
    !   Distributed loads are specified in the input file using the Un option in the input file.
    !   n specifies the face number, following the ABAQUS convention
    !
    !
      do j = 1,NDLOAD

        call abq_facenodes_3D(NNODE,iabs(JDLTYP(j,1)),
     1                                     face_node_list,nfacenodes)

        do i = 1,nfacenodes
            face_coords(1:3,i) = coords(1:3,face_node_list(i))
        end do

        if (nfacenodes == 3) n_points = 3
        if (nfacenodes == 6) n_points = 4
        if (nfacenodes == 4) n_points = 4
        if (nfacenodes == 8) n_points = 9

        call abq_UEL_2D_integrationpoints(n_points, nfacenodes, xi2, w)

        do kint = 1,n_points
            call abq_UEL_2D_shapefunctions(xi2(1:2,kint),
     1                        nfacenodes,N2,dNdxi2)
            dxdxi2 = matmul(face_coords(1:3,1:nfacenodes),
     1                           dNdxi2(1:nfacenodes,1:2))
            norm(1)=(dxdxi2(2,1)*dxdxi2(3,2))-(dxdxi2(2,2)*dxdxi2(3,1))
            norm(2)=(dxdxi2(1,1)*dxdxi2(3,2))-(dxdxi2(1,2)*dxdxi2(3,1))
            norm(3)=(dxdxi2(1,1)*dxdxi2(2,2))-(dxdxi2(1,2)*dxdxi2(2,1))

            do i = 1,nfacenodes
                ipoin = 3*face_node_list(i)-2
                RHS(ipoin:ipoin+2,1) = RHS(ipoin:ipoin+2,1)
     1                 - N2(1:nfacenodes)*adlmag(j,1)*norm(1:3)*w(kint)      ! Note determinant is already in normal
            end do
        end do
      end do

      return

      END SUBROUTINE UEL_3D


      subroutine hyper(element_properties,n_properties,F,JJ,B,Stress,D)
       
       implicit none
       
       integer, intent(in)           :: n_properties
       double precision, intent(in)  :: element_properties(n_properties)
       double precision, intent(in)  :: F(3,3)
       double precision, intent(in)  :: JJ
       double precision, intent(in)  :: B(3,3)
       double precision, intent(out) :: Stress(6)
       double precision, intent(out) :: D(6,6)
       
       double precision :: mu
       double precision :: K
       double precision :: G(6,6)
       double precision :: C(6)
       double precision :: Binv(3,3)
       double precision :: JB       ! Jacobian of B
       double precision :: Cinv(6)
       double precision :: Cstar(6)
       double precision :: Cbar(6)
       double precision :: Cstarbar(6)
       double precision :: Aye(6) ! Identity Vector
       double precision :: Pvec(6)
       double precision :: Q
       double precision :: Omega(6,6)
       
       ! Subroutine to calculate the material stress Sigma = Finverse*J*cauchystress*Finversetranspose as well as
       ! the tangent matrix D
       ! Define the element Properties
       mu = element_properties(1)
       K  = element_properties(2)
       G = 0.d0
       G(1,1) = element_properties(3)
       G(2,2) = element_properties(4)
       G(3,3) = element_properties(5)
       G(4,4) = element_properties(6)
       G(5,5) = G(4,4)
       G(6,6) = G(4,4)
       
       ! Define C vector From Cauchy-Green Tensor
       C(1) = B(1,1)
       C(2) = B(2,2)
       C(3) = B(3,3)
       C(4) = B(1,2)
       C(5) = B(1,3)
       C(6) = B(2,3)
       ! Define Cinv vector from the inverse of the Cauchy-Green Tensor
       call abq_UEL_invert3d(B,Binv,JB)
       Cinv(1) = Binv(1,1)
       Cinv(2) = Binv(2,2)
       Cinv(3) = Binv(3,3)
       Cinv(4) = Binv(1,2)
       Cinv(5) = Binv(1,3)
       Cinv(6) = Binv(2,3)
       ! Define Cstar with shear terms doubled
       Cstar(1) = B(1,1)
       Cstar(2) = B(2,2)
       Cstar(3) = B(3,3)
       Cstar(4) = 2.d0*B(1,2)
       Cstar(5) = 2.d0*B(1,3)
       Cstar(6) = 2.d0*B(2,3)
       ! Define Cbar by normalizing C with the Jacobian
       Cbar(1) = C(1)/JJ**(2.d0/3.d0)
       Cbar(2) = C(2)/JJ**(2.d0/3.d0)
       Cbar(3) = C(3)/JJ**(2.d0/3.d0)
       Cbar(4) = C(4)/JJ**(2.d0/3.d0)
       Cbar(5) = C(5)/JJ**(2.d0/3.d0)
       Cbar(6) = C(6)/JJ**(2.d0/3.d0)
       ! Define Cstarbar by normalizing Cstar with the Jacobian
       Cstarbar(1) = Cstar(1)/JJ**(2.d0/3.d0)
       Cstarbar(2) = Cstar(2)/JJ**(2.d0/3.d0)
       Cstarbar(3) = Cstar(3)/JJ**(2.d0/3.d0)
       Cstarbar(4) = Cstar(4)/JJ**(2.d0/3.d0)
       Cstarbar(5) = Cstar(5)/JJ**(2.d0/3.d0)
       Cstarbar(6) = Cstar(6)/JJ**(2.d0/3.d0)
       
       ! Define the Pvec using G, Cstarbar, Cstar, Cinv, and Identity matrix
       Aye = 0.d0
       Aye(1) = 1.d0
       Aye(2) = 1.d0
       Aye(3) = 1.d0
       Aye(4) = 0.d0
       Aye(5) = 0.d0
       Aye(6) = 0.d0
       
       Pvec = (1.d0/(2.d0*JJ**(2.d0/3.d0)))*(matmul(G,(Cstarbar-Aye))
     1  -(1.d0/3.d0)*(dot_product(Cstar,matmul(G,(Cstarbar-Aye))))*Cinv)
       
       ! Define Q
       Q = (1.d0/4.d0)*dot_product((Cbar-Aye),matmul(G,(Cbar-Aye)))
       
       ! Define Omega
       Omega = 0.d0
       Omega(1,1) = Cinv(1)*Cinv(1)
       Omega(1,2) = Cinv(4)*Cinv(4)
       Omega(1,3) = Cinv(5)*Cinv(5)
       Omega(1,4) = Cinv(1)*Cinv(4)
       Omega(1,5) = Cinv(1)*Cinv(5)
       Omega(1,6) = Cinv(4)*Cinv(5)
       
       Omega(2,1) = Omega(1,2)
       Omega(2,2) = Cinv(2)*Cinv(2)
       Omega(2,3) = Cinv(6)*Cinv(6)
       Omega(2,4) = Cinv(4)*Cinv(2)
       Omega(2,5) = Cinv(4)*Cinv(6)
       Omega(2,6) = Cinv(2)*Cinv(6)
       
       Omega(3,1) = Omega(1,3)
       Omega(3,2) = Omega(2,3)
       Omega(3,3) = Cinv(3)*Cinv(3)
       Omega(3,4) = Cinv(5)*Cinv(6)
       Omega(3,5) = Cinv(5)*Cinv(3)
       Omega(3,6) = Cinv(6)*Cinv(3)
       
       Omega(4,1) = Omega(1,4)
       Omega(4,2) = Omega(2,4)
       Omega(4,3) = Omega(3,4)
       Omega(4,4) = (Cinv(1)*Cinv(2)+ Cinv(4)*Cinv(4))/2.d0
       Omega(4,5) = (Cinv(1)*Cinv(6)+ Cinv(5)*Cinv(4))/2.d0
       Omega(4,6) = (Cinv(4)*Cinv(6)+ Cinv(5)*Cinv(2))/2.d0
       
       Omega(5,1) = Omega(1,5)
       Omega(5,2) = Omega(2,5)
       Omega(5,3) = Omega(3,5)
       Omega(5,4) = Omega(4,5)
       Omega(5,5) = (Cinv(1)*Cinv(3)+ Cinv(5)*Cinv(5))/2.d0
       Omega(5,6) = (Cinv(4)*Cinv(3)+ Cinv(5)*Cinv(6))/2.d0
       
       Omega(6,1) = Omega(1,6)
       Omega(6,2) = Omega(2,6)
       Omega(6,3) = Omega(3,6)
       Omega(6,4) = Omega(4,6)
       Omega(6,5) = Omega(5,6)
       Omega(6,6) = (Cinv(2)*Cinv(3)+ Cinv(6)*Cinv(6))/2.d0
       
       ! Define the Material Stress
       
       Stress = 0.d0
       Stress = Stress + mu*exp(Q)*Pvec
       Stress = Stress + K*JJ*(JJ-1.d0)*Cinv
       
       ! Define the D matrix
      end subroutine hyper

