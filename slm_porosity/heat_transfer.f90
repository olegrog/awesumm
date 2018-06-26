MODULE user_case
  ! CASE 2/3 D channel periodic in x-direction and if dim=3 in z-derection
  USE precision
  USE wlt_vars
  USE wavelet_filters_mod
  USE elliptic_mod
  USE elliptic_vars
  USE wlt_trns_mod
  USE wlt_trns_vars
  USE wlt_trns_util_mod
  USE io_3d_vars
  USE util_mod
  USE util_vars
  USE share_consts
  USE pde
  USE variable_mapping
  USE sizes
  USE share_kry
  USE vector_util_mod
  USE field
  USE input_file_reader
  USE debug_vars
  !
  ! case specific variables
  !
  INTEGER n_var_pressure
  INTEGER n_var_temp
  INTEGER n_var_enthalpy
  INTEGER n_var_lfrac
  INTEGER n_var_porosity
  INTEGER smoothing_method

  REAL (pr) :: fusion_delta         ! = liquidus - solidus
  REAL (pr) :: fusion_heat          ! the latent heat of fusion
  REAL (pr) :: convective_transfer  ! convective heat transfer coefficient
  REAL (pr) :: radiative_transfer   ! radiative heat transfer coefficient
  REAL (pr) :: power                ! laser power
  REAL (pr) :: absorptivity
  REAL (pr) :: scanning_speed
  REAL (pr) :: initial_porosity
  REAL (pr) :: initial_enthalpy
  REAL (pr) :: smoothing_width      ! width of spline in terms of fusion_delta
  REAL (pr) :: eps_zero
  REAL (pr) :: Dconductivity_solid
  REAL (pr) :: Dconductivity_liquid
  REAL (pr) :: conductivity_fusion
  REAL (pr) :: Dcapacity_solid
  REAL (pr) :: Dcapacity_liquid
  REAL (pr) :: capacity_fusion
  REAL (pr) :: emissivity

  REAL (pr) :: enthalpy_S
  REAL (pr) :: enthalpy_L
  REAL (pr), DIMENSION(3) :: x0     ! Initial coordinates of the center of the laser beam
CONTAINS

  !
  ! In user_setup_pde() we setup how many variables are integrated, which are interpolated
  ! to the next times step and if any exeact solution exists to compare to during the run.
  ! We also set the variable names that are used when the result files are written.
  !
  ! The following variables must be setup in this routine:
  !
  !
  ! n_integrated     ! first n_integrated eqns will be acted on for time integration
  ! n_var_additional !
  ! n_var
  !
  !
  !
  !
  SUBROUTINE  user_setup_pde ( VERB )
    USE variable_mapping
    IMPLICIT NONE
    LOGICAL, OPTIONAL :: VERB         ! debug info
    INTEGER :: i

    !register_var( var_name, integrated, adapt, interpolate, exact, saved, req_restart, FRONT_LOAD )
    CALL register_var( 'enthalpy ', &
        integrated=.TRUE., &
        adapt=(/.TRUE.,.TRUE./), &
        interpolate=(/.TRUE.,.TRUE./), &
        exact=(/.FALSE.,.FALSE./), &
        saved=.TRUE., &
        req_restart=.TRUE.)
    CALL register_var( 'temperature', &
        integrated=.FALSE., &
        adapt=(/.FALSE.,.FALSE./), &
        interpolate=(/.FALSE.,.FALSE./), &
        exact=(/.FALSE.,.FALSE./), &
        saved=.TRUE., &
        req_restart=.FALSE.)
    CALL register_var( 'liquid_fraction', &
        integrated=.FALSE., &
        adapt=(/.FALSE.,.FALSE./), &
        interpolate=(/.FALSE.,.FALSE./), &
        exact=(/.FALSE.,.FALSE./), &
        saved=.TRUE., &
        req_restart=.FALSE.)
    CALL register_var( 'porosity', &
        integrated=.FALSE., &
        adapt=(/.TRUE.,.TRUE./), &
        interpolate=(/.FALSE.,.FALSE./), &
        exact=(/.FALSE.,.FALSE./), &
        saved=.TRUE., &
        req_restart=.FALSE.)
    CALL register_var( 'pressure ', &
        integrated=.FALSE., &
        adapt=(/.FALSE.,.FALSE./), &
        interpolate=(/.FALSE.,.FALSE./), &
        exact=(/.FALSE.,.FALSE./), &
        saved=.TRUE., &
        req_restart=.FALSE.)
    CALL setup_mapping()
    CALL print_variable_registery( FULL=.TRUE.)

    n_var_pressure  = get_index('pressure ')
    n_var_temp      = get_index('temperature')
    n_var_enthalpy  = get_index('enthalpy')
    n_var_lfrac     = get_index('liquid_fraction')
    n_var_porosity  = get_index('porosity')

    ALLOCATE ( Umn(1:n_var) )
    Umn = 0.0_pr !set up here if mean quantities are not zero and used in scales or equation

    IF (verb_level.GT.0) THEN
       PRINT *, 'n_integrated = ',n_integrated
       PRINT *, 'n_var = ',n_var
       PRINT *, 'n_var_exact = ',n_var_exact
       PRINT *, '*******************Variable Names*******************'
       DO i = 1,n_var
          WRITE (*, u_variable_names_fmt) u_variable_names(i)
       END DO
       PRINT *, '****************************************************'
    END IF

  END SUBROUTINE  user_setup_pde
  !
  ! Set the exact solution for comparison to the simulated solution
  !
  ! u          - array to fill in the exact solution
  ! nlocal       - number of active wavelets
  ! ne_local        - total number of equations
  ! t          - time of current time step
  ! l_n_var_exact_soln_index - index into the elements of u for which we need to
  !                            find the exact solution
  SUBROUTINE  user_exact_soln (u, nlocal,  t_local, l_n_var_exact_soln)
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: nlocal
    REAL (pr), INTENT (IN) ::  t_local
    REAL (pr), DIMENSION (nlocal,n_var_exact), INTENT (INOUT) :: u
    LOGICAL , INTENT (IN) :: l_n_var_exact_soln(n_var)
    REAL (pr) :: t_zero
    INTEGER :: i
    DOUBLE PRECISION :: DERF
    EXTERNAL :: DERF
    u = 0
  END SUBROUTINE  user_exact_soln

  SUBROUTINE user_initial_conditions (u, nlocal, ne_local, t_local, scl, scl_fltwt, iter)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nlocal, ne_local
    INTEGER, INTENT (INOUT) :: iter ! iteration of call while adapting initial grid
    REAL (pr), DIMENSION (nlocal,ne_local), INTENT (INOUT) :: u
    REAL (pr)  :: scl(1:n_var),scl_fltwt
    REAL (pr), INTENT (IN) :: t_local
    REAL (pr) :: t_zero
    INTEGER :: i

    IF (dim.EQ.2) x0(dim) = xyzlimits(1,dim)
    IF ( IC_restart_mode.EQ.0 ) THEN
       DO i = 1, nlocal
          u(i,n_var_enthalpy) = initial_enthalpy*EXP(-SUM((x(i,:)-x0)**2))*EXP(-(x(i,dim)-x0(dim))**2*power*absorptivity)
       END DO
       !WHERE (x(:,dim).NE.x0(dim))
       !   u(:,n_var_enthalpy) = 0
       !END WHERE
    END IF

  END SUBROUTINE user_initial_conditions

!--********************************
  ! Arguments
  ! u         - field on adaptive grid
  ! nlocal      - number of active points
  ! ne_local       - number of equations
  ! t         - current time
  !
  ! Global variables used
  ! ifn       - number of equations
  ! x()
  ! xO
  ! y0
  !--********************************
  SUBROUTINE user_algebraic_BC (Lu, u, nlocal, ne_local, jlev, meth)
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: jlev, meth, ne_local, nlocal
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: Lu
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (IN)    :: u

    INTEGER :: i, ie, shift, face_type, nloc
    REAL (pr), DIMENSION (ne_local,nlocal,dim) :: du, d2u
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc

    CALL c_diff_fast (u, du, d2u, jlev, nlocal, meth, 10, ne_local, 1, ne_local)

    DO ie = 1, ne_local
       shift = nlocal*(ie-1)
       i_p_face(0) = 1
       DO i=1,dim
          i_p_face(i) = i_p_face(i-1)*3
       END DO
       DO face_type = 0, 3**dim - 1
          face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
          IF( ANY( face(1:dim) /= 0) .AND. ie == n_var_enthalpy ) THEN
             CALL get_all_indices_by_face (face_type, jlev, nloc, iloc)
             iloc(1:nloc) = shift + iloc(1:nloc)
             IF( nloc > 0 ) THEN
                IF( face(dim) < 0 ) THEN    ! z=0 face
                  ! dependence on temperature should be linear; therefore, we use only u_prev_timestep here
                  Lu(iloc(1:nloc)) = Neumann_bc(u_prev_timestep(iloc(1:nloc))) * du(ie, iloc(1:nloc), dim) + &
                     Dirichlet_bc(u_prev_timestep(iloc(1:nloc))) * u(iloc(1:nloc))
                ELSE                        ! other faces
                   Lu(iloc(1:nloc)) = u(iloc(1:nloc))
                END IF
             END IF
          END IF
       END DO
    END DO
    !PRINT *, 'user_algebraic_BC', Lu(shift+iloc(1:nloc))

  END SUBROUTINE user_algebraic_BC

  SUBROUTINE user_algebraic_BC_diag (Lu_diag, nlocal, ne_local, jlev, meth)
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: jlev, meth, ne_local, nlocal
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: Lu_diag

    INTEGER :: i, ie, shift, face_type, nloc
    REAL (pr), DIMENSION (nlocal,dim) :: du, d2u
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc

    CALL c_diff_diag (du, d2u, jlev, nlocal, meth, meth, 10)

    DO ie = 1, ne_local
       shift = nlocal*(ie-1)
       i_p_face(0) = 1
       DO i=1,dim
          i_p_face(i) = i_p_face(i-1)*3
       END DO
       DO face_type = 0, 3**dim - 1
          face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
          IF( ANY( face(1:dim) /= 0) .AND. ie == n_var_enthalpy ) THEN
             CALL get_all_indices_by_face (face_type, jlev, nloc, iloc)
             iloc(1:nloc) = shift + iloc(1:nloc)
             IF( nloc > 0 ) THEN
                IF( face(dim) < 0 ) THEN    ! z=0 face
                   Lu_diag(iloc(1:nloc)) = Neumann_bc(u_prev_timestep(iloc(1:nloc))) * du(iloc(1:nloc), dim) + &
                      Dirichlet_bc(u_prev_timestep(iloc(1:nloc)))
                ELSE                        ! other faces
                   Lu_diag(iloc(1:nloc)) = 1.0_pr
                END IF
             END IF
          END IF
       END DO
    END DO
    !PRINT *, 'user_algebraic_BC_diag', Lu_diag(shift+iloc(1:nloc))

  END SUBROUTINE user_algebraic_BC_diag

  SUBROUTINE user_algebraic_BC_rhs (rhs, ne_local, nlocal, jlev)
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: ne_local, nlocal, jlev
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: rhs
    REAL (pr), DIMENSION (nlocal*ne_local) :: T, DT

    INTEGER :: i, ie, shift, face_type, nloc, meth=1
    REAL (pr), DIMENSION (ne_local,nlocal,dim) :: du, d2u
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc

    DO ie = 1, ne_local
       shift = nlocal*(ie-1)
       i_p_face(0) = 1
       DO i=1,dim
          i_p_face(i) = i_p_face(i-1)*3
       END DO
       DO face_type = 0, 3**dim - 1
          face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
          IF( ANY( face(1:dim) /= 0) .AND. ie == n_var_enthalpy ) THEN
             CALL get_all_indices_by_face (face_type, jlev, nloc, iloc)
             iloc(1:nloc) = shift + iloc(1:nloc)
             IF( nloc > 0 ) THEN
                IF( face(dim) < 0 ) THEN    ! z=0 face
                   IF( dim == 3 ) THEN
                      rhs(iloc(1:nloc)) = -exp(-(x(iloc(1:nloc), 1) - scanning_speed*t - x0(1))**2 - (x(iloc(1:nloc), 2) - x0(2))**2 ) / pi
                   ELSE
                      rhs(iloc(1:nloc)) = -exp(-(x(iloc(1:nloc), 1) - scanning_speed*t - x0(1))**2 ) / pi**.5
                   END IF
                   T = temperature(u_prev_timestep)
                   DT = Dtemperature(u_prev_timestep, T)
                   rhs(iloc(1:nloc)) = rhs(iloc(1:nloc)) * absorptivity * power - F_heat_flux(T(iloc(1:nloc))) + &
                      F_heat_flux(T(iloc(1:nloc)), .TRUE.)*DT(iloc(1:nloc))*u_prev_timestep(iloc(1:nloc))
                ELSE                        ! other faces
                   rhs(iloc(1:nloc)) = 0
                END IF
             END IF
          END IF
       END DO
    END DO

  END SUBROUTINE user_algebraic_BC_rhs

  SUBROUTINE user_project (u, p, nlocal, meth)
    !--Makes u divergence free
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: meth, nlocal
    REAL (pr), DIMENSION (nlocal),     INTENT(INOUT) :: p
    REAL (pr), DIMENSION (nlocal,n_integrated), INTENT(INOUT) :: u
  END SUBROUTINE user_project

  FUNCTION user_rhs (u_integrated, p)
    IMPLICIT NONE
    REAL (pr), DIMENSION (ng,ne), INTENT(IN) :: u_integrated
    REAL (pr), DIMENSION (ng), INTENT(IN) :: p
    REAL (pr), DIMENSION (n) :: user_rhs
    REAL (pr), DIMENSION (ng) :: T, phi
    INTEGER :: ie, shift, i
    INTEGER, PARAMETER :: meth = 1
    REAL (pr), DIMENSION (ne,ng,dim) :: du, du_dummy
    REAL (pr), DIMENSION (dim,ng,dim):: d2u, d2u_dummy

    ie = n_var_enthalpy
    shift = ng*(ie-1)
    IF (IMEXswitch.LE.0) THEN
       user_rhs(shift+1:shift+ng) = 0.0_pr
    END IF
    IF (IMEXswitch .GE. 0) THEN
       CALL c_diff_fast(enthalpy_fusion(u_integrated(:,ie)), du, du_dummy, j_lev, ng, meth, 10, ne, 1, ne)
       T = temperature(u_integrated(:,ie))
       phi = liquid_fraction(u_integrated(:,ie))
       du(ie,:,:) = du(ie,:,:) * SPREAD(conductivity(T, phi) / capacity(T, phi) * porosity_term(), 2, dim)
       CALL c_diff_fast(du(ie,:,:), d2u, d2u_dummy, j_lev, ng, meth, 10, dim, 1, dim)
       DO i = 1, dim
          user_rhs(shift+1:shift+ng) = user_rhs(shift+1:shift+ng) + d2u(i,:,i)
       END DO
    END IF
  END FUNCTION user_rhs

  FUNCTION user_Drhs (pert_u, u_prev, meth)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: meth
    REAL (pr), DIMENSION (ng,ne) :: pert_u, u_prev
    REAL (pr), DIMENSION (n) :: user_Drhs
    REAL (pr), DIMENSION (ng) :: T, phi, Dphi, DT, k_0, c_p, Dk_0, Dc_p
    INTEGER :: ie, shift, i
    REAL (pr), DIMENSION (ng,2*ne) :: for_du
    REAL (pr), DIMENSION (2*ne,ng,dim) :: du, du_dummy ! der(pert_u) in du(1:ne,:,:) and der(u) in du(ne+1:2*ne,:,:)
    REAL (pr), DIMENSION (dim,ng,dim) :: d2u, d2u_dummy
    REAL (pr), DIMENSION (ng,dim) :: part1, part2

    ie = n_var_enthalpy
    shift = ng*(ie-1)
    IF (IMEXswitch.LE.0) THEN
       user_Drhs(shift+1:shift+ng) = 0.0_pr
    END IF
    IF (IMEXswitch .GE. 0) THEN
       T = temperature(u_prev(:,ie))
       DT = Dtemperature(u_prev(:,ie), T)
       phi = liquid_fraction(u_prev(:,ie))
       Dphi = liquid_fraction(u_prev(:,ie), .TRUE.)
       k_0 = conductivity(T, phi)
       c_p = capacity(T, phi)
       Dk_0 = conductivity(T, phi, 1)*Dphi + conductivity(T, phi, 2)*DT
       Dc_p = capacity(T, phi, 1)*Dphi + capacity(T, phi, 2)*DT
       for_du = RESHAPE((/ enthalpy_fusion(u_prev(:,ie), .TRUE.)*pert_u(:,ie), enthalpy_fusion(u_prev(:,ie)) /), SHAPE(for_du))
       CALL c_diff_fast(for_du, du, du_dummy, j_lev, ng, meth, 10, 2*ne, 1, 2*ne)
       part1 = du(2*ie,:,:) * SPREAD((Dk_0*c_p - Dc_p*k_0) / c_p**2 * enthalpy_fusion(u_prev(:,ie), .TRUE.)*pert_u(:,ie), 2, dim)
       part2 = du(ie,:,:) * SPREAD(k_0 / c_p, 2, dim)
       CALL c_diff_fast((part1 + part2) * SPREAD(porosity_term(), 2, dim), d2u, d2u_dummy, j_lev, ng, meth, 10, dim, 1, dim)
       DO i = 1, dim
          user_Drhs(shift+1:shift+ng) = user_Drhs(shift+1:shift+ng) + d2u(i,:,i)
       END DO
    END IF
    IF (user_Drhs(1).NE.user_Drhs(1)) THEN
       PRINT *, '--- NaN in user_Drhs ---'
       CALL ABORT
    END IF

  END FUNCTION user_Drhs

  FUNCTION user_Drhs_diag (meth)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: meth
    REAL (pr), DIMENSION (n) :: user_Drhs_diag
    REAL (pr), DIMENSION (ng) :: T, phi, Dphi, DT, k_0, c_p, Dk_0, dc_p
    INTEGER :: ie, shift, i
    REAL (pr), DIMENSION (ng,dim) :: du, d2u, part1, part2
    REAL (pr), DIMENSION (ne,ng,dim) :: du_prev, du_dummy

    ie = n_var_enthalpy
    shift = ng*(ie-1)
    IF (IMEXswitch.LE.0) THEN
       user_Drhs_diag(shift+1:shift+ng) = 1.0_pr
    END IF
    IF (IMEXswitch .GE. 0) THEN
       T = temperature(u_prev_timestep(shift+1:shift+ng))
       phi = liquid_fraction(u_prev_timestep(shift+1:shift+ng))
       Dphi = liquid_fraction(u_prev_timestep(shift+1:shift+ng), .TRUE.)
       DT = Dtemperature(u_prev_timestep(shift+1:shift+ng), T)
       k_0 = conductivity(T, phi)
       c_p = capacity(T, phi)
       Dk_0 = conductivity(T, phi, 1)*Dphi + conductivity(T, phi, 2)*DT
       Dc_p = capacity(T, phi, 1)*Dphi + capacity(T, phi, 2)*DT
       CALL c_diff_fast(enthalpy_fusion(u_prev_timestep(shift+1:shift+ng)), du_prev, du_dummy, j_lev, ng, meth, 10, ne, 1, ne)
       CALL c_diff_diag(du, d2u, j_lev, ng, meth, meth, -11)
       part1 = du * du_prev(ie,:,:) * SPREAD((Dk_0*c_p - Dc_p*k_0) / c_p**2, 2, dim)
       part2 = d2u * SPREAD(k_0 / c_p, 2, dim)
       user_Drhs_diag(shift+1:shift+ng) = SUM(part1 + part2, 2) * &
          enthalpy_fusion(u_prev_timestep(shift+1:shift+ng), .TRUE.) * porosity_term()
    END IF

  END FUNCTION user_Drhs_diag

  FUNCTION user_chi (nlocal, t_local )
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nlocal
    REAL (pr), INTENT (IN) :: t_local
    REAL (pr), DIMENSION (nlocal) :: user_chi

    user_chi = 0.0_pr
  END FUNCTION user_chi


  FUNCTION user_mapping ( xlocal, nlocal, t_local )
    USE curvilinear
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nlocal
    REAL (pr), INTENT (IN) :: t_local
    REAL (pr), DIMENSION (nlocal,dim), INTENT(IN) :: xlocal
    REAL (pr), DIMENSION (nlocal,dim) :: user_mapping

    user_mapping(:,1:dim) = xlocal(:,1:dim)

  END FUNCTION user_mapping
  !
  ! Calculate any statitics
  !
  ! startup_flag - 0 when adapting to IC,then 1 inmain integration loop
  !
  SUBROUTINE user_stats ( u ,j_mn, startup_flag)
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: startup_flag
    INTEGER , INTENT (IN) :: j_mn
    REAL (pr), DIMENSION (nwlt,1:n_var), INTENT (IN) :: u
  END SUBROUTINE user_stats


  SUBROUTINE user_cal_force (u, n, t_local, force, drag, lift)
    !--Calculates drag and lift on obstacle using penalization formula
    IMPLICIT NONE

    INTEGER, INTENT (IN) :: n
    REAL (pr), INTENT (IN) :: t_local
    REAL (pr), DIMENSION (dim), INTENT (INOUT) :: force
    REAL (pr), INTENT (OUT) :: drag, lift
    REAL (pr), DIMENSION (n,dim) :: u
    drag = 0.0_PR
    lift = 0.0_PR

    !
    ! There is no obstacle in flow
    !

  END SUBROUTINE user_cal_force

  !
  ! Read input from "case_name"_pde.inp file
  ! case_name is in string file_name read from command line
  ! in read_command_line_input()
  !
  SUBROUTINE user_read_input()
    IMPLICIT NONE
    REAL(pr), DIMENSION(1) :: tmp

    call input_real ('fusion_delta', fusion_delta, 'stop')
    call input_real ('fusion_heat', fusion_heat, 'stop')
    call input_real ('convective_transfer', convective_transfer, 'stop')
    call input_real ('power', power, 'stop')
    call input_real ('absorptivity', absorptivity, 'stop')
    call input_real ('scanning_speed', scanning_speed, 'stop')
    call input_real ('initial_porosity', initial_porosity, 'stop')
    call input_real ('initial_enthalpy', initial_enthalpy, 'stop')
    call input_real ('smoothing_width', smoothing_width, 'stop')

    call input_real ('Dconductivity_liquid', Dconductivity_liquid, 'stop')
    call input_real ('Dconductivity_solid', Dconductivity_solid, 'stop')
    call input_real ('conductivity_fusion', conductivity_fusion, 'stop')
    call input_real ('Dcapacity_liquid', Dcapacity_liquid, 'stop')
    call input_real ('Dcapacity_solid', Dcapacity_solid, 'stop')
    call input_real ('capacity_fusion', capacity_fusion, 'stop')
    call input_real ('eps_zero', eps_zero, 'stop')
    call input_real ('emissivity', emissivity, 'stop')
    call input_real ('radiative_transfer', radiative_transfer, 'stop')

    call input_real_vector ('x0', x0, 3, 'stop')

    call input_integer ('smoothing_method', smoothing_method, 'stop')

    enthalpy_S = 0.0_pr
    enthalpy_L = 1.0_pr
    tmp = enthalpy((/ 1.0_pr - fusion_delta/2 /), (/ liquid_fraction((/ enthalpy_S /)) /))
    enthalpy_S = tmp(1)
    tmp = enthalpy((/ 1.0_pr + fusion_delta/2 /), (/ liquid_fraction((/ enthalpy_L /)) /))
    enthalpy_L = tmp(1)
    PRINT *, 'enthalpy_S', enthalpy_S
    PRINT *, 'enthalpy_L', enthalpy_L
  END SUBROUTINE user_read_input



  !
  ! calculate any additional variables
  !
  ! arg
  ! flag  - 0 calledwhile adapting to IC, 1 - called in main time integration loop
  !
  ! These additional variables are calculated and left in real space.
  !
  SUBROUTINE user_additional_vars( t_local, flag )
    IMPLICIT NONE
    REAL (pr), INTENT (IN) ::  t_local
    INTEGER , INTENT(IN) :: flag ! 0- called during adaption to IC, 1 called during main integration loop
    IF (.NOT.flag) THEN
       u(:,n_var_porosity) = initial_porosity
    END IF
    u(:,n_var_temp) = temperature(u(:,n_var_enthalpy))
    u(:,n_var_lfrac) = liquid_fraction(u(:,n_var_enthalpy))
    u(:,n_var_porosity) = porosity(u(:,n_var_enthalpy))
    u(:,n_var_pressure) = 0.0_pr
  END SUBROUTINE user_additional_vars


  !
  ! calculate any additional scalar variables
  !
  ! arg
  ! flag  - 0 calledwhile adapting to IC, 1 - called in main time integration loop
  !
  SUBROUTINE user_scalar_vars( flag )
    IMPLICIT NONE
    INTEGER , INTENT(IN) :: flag ! 0- called during adaption to IC, 1 called during main integration loop



  END SUBROUTINE user_scalar_vars

  !
  !************ Calculating Scales ***************************
  !
  ! Note the order of the components in the scl array
  ! correspond to u_tn, v_tn, w_tn, u_tn-1, v_tn-1, w_tn-1, u_tn-2, v_tn-2, w_tn-2
  !
  SUBROUTINE user_scales(flag, use_default, u, nlocal, ne_local, l_n_var_adapt , l_n_var_adapt_index, &
       scl, scl_fltwt) !add
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: flag ! 0 during initial adaptation, then 1 during time advancement
    LOGICAL , INTENT(INOUT) :: use_default
    INTEGER, INTENT (IN) :: nlocal, ne_local
    REAL (pr), DIMENSION (1:nlocal,1:ne_local), INTENT (IN) :: u
    LOGICAL , INTENT (IN) :: l_n_var_adapt(ne_local)
    INTEGER , INTENT (IN) :: l_n_var_adapt_index(1:ne_local)
    REAL (pr), DIMENSION (1:ne_local), INTENT (INOUT) :: scl
    REAL (pr) ,  INTENT(IN) ::scl_fltwt !weight for temporal filter on scl

    !
    ! Ignore the output of this routine and use default scales routine
    !
    use_default = .TRUE.

    !
    ! NOTE: For a parallel run, synchronize scl(:) across the processors in user_scales.
    !       Use the subroutine scales of default_util.f90 as an example,
    !       subroutine parallel_global_sum of parallel.f90 may be helpful.
    !       Already synchronized: sumdA_global
  END SUBROUTINE user_scales

  SUBROUTINE user_cal_cfl (use_default, u, cfl_out)
    USE precision
    USE sizes
    USE pde
    IMPLICIT NONE
    LOGICAL , INTENT(INOUT) :: use_default
    REAL (pr),                                INTENT (INOUT) :: cfl_out
    REAL (pr), DIMENSION (nwlt,n_integrated), INTENT (IN)    :: u
    INTEGER                    :: i
    REAL (pr), DIMENSION (dim) :: cfl
    REAL (pr), DIMENSION(dim,nwlt) :: h_arr

    use_default = .FALSE.
    CALL get_all_local_h (h_arr)
    cfl_out = MAXVAL(dt/h_arr(1,:)*scanning_speed)
    PRINT *, "cfl_convective", cfl_out
    DO i = 1, nwlt
       cfl(1:dim) = ABS (u(i,1:dim)+Umn(1:dim)) * dt/h_arr(1:dim,i)
       cfl_out = MAX (cfl_out, MAXVAL(cfl))
    END DO

  !  cfl_out = 1.0_pr ! no CFL condition
    IF (u(1,1).NE.u(1,1)) THEN
      PRINT *, '--- INFINITE ---'
      CALL ABORT
    END IF

  END SUBROUTINE user_cal_cfl

  SUBROUTINE user_init_sgs_model( )
    IMPLICIT NONE
  END SUBROUTINE user_init_sgs_model

  SUBROUTINE  user_sgs_force (u_loc, nlocal)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nlocal
    REAL (pr), DIMENSION (nlocal,n_integrated), INTENT (INOUT) :: u_loc
  END SUBROUTINE  user_sgs_force

  FUNCTION user_sound_speed (u, neq, nwlt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nwlt, neq
    REAL (pr), DIMENSION (nwlt,neq), INTENT(IN) :: u
    REAL (pr), DIMENSION (nwlt) :: user_sound_speed

    user_sound_speed(:) = 0.0_pr
    !user_sound_speed(:) = SQRT( gamma*(gamma-1.0_pr)* &
    !       ( u(:,4)-0.5_pr*SUM(u(:,n_var_mom(1:dim))**2,DIM=2)/u(:,1)) ) ! pressure

  END FUNCTION user_sound_speed

  SUBROUTINE  user_pre_process
    IMPLICIT NONE
  END SUBROUTINE user_pre_process

  SUBROUTINE  user_post_process
    IMPLICIT NONE
  END SUBROUTINE user_post_process

  FUNCTION liquid_fraction (enthalpy, is_D)
    IMPLICIT NONE
    REAL (pr),         INTENT(IN) :: enthalpy(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: is_D
    REAL (pr) :: liquid_fraction(SIZE(enthalpy))

    IF (smoothing_method.EQ.0) THEN         ! C^0
      liquid_fraction = lf_piecewise(enthalpy, is_D)
    ELSE IF (smoothing_method.EQ.1) THEN    ! C^1
      liquid_fraction = lf_cubic_splines(enthalpy, is_D)
    ELSE IF (smoothing_method.EQ.2) THEN    ! C^\infty
      liquid_fraction = lf_exponent(enthalpy, is_D)
    END IF
  END FUNCTION liquid_fraction

  FUNCTION lf_piecewise (enthalpy, is_D)
    IMPLICIT NONE
    REAL (pr),         INTENT(IN) :: enthalpy(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: is_D
    REAL (pr) :: lf_piecewise(SIZE(enthalpy))
    REAL (pr) :: lf_der

    lf_der = 1.0_pr/(enthalpy_L - enthalpy_S)           ! first derivative of liquid fraction on enthalpy

    IF (.NOT.PRESENT(is_D).OR.(.NOT.is_D)) THEN
      lf_piecewise = MAX(0.0_pr, MIN(1.0_pr, (enthalpy - enthalpy_S)*lf_der))
    ELSE
      WHERE (enthalpy_S.LT.enthalpy .AND. enthalpy.LT.enthalpy_L)
        lf_piecewise = lf_der
      ELSEWHERE
        lf_piecewise = 0.0_pr
      END WHERE
    END IF
  END FUNCTION lf_piecewise

  FUNCTION lf_exponent (enthalpy, is_D)
    IMPLICIT NONE
    REAL (pr),         INTENT(IN) :: enthalpy(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: is_D
    REAL (pr) :: lf_exponent(SIZE(enthalpy))
    REAL (pr) :: enthalpy_delta, halfpoint, lf_der

    enthalpy_delta = fusion_delta + fusion_heat         ! change of enthalpy for fusion
    halfpoint = 1.0_pr + enthalpy_delta/2.0_pr          ! a middle (symmetric) point of fusion in terms of enthalpy
    lf_der = 1.0_pr/enthalpy_delta                      ! first derivative of liquid fraction on enthalpy

    lf_exponent = EXP(-4*lf_der*(enthalpy - halfpoint))

    IF (.NOT.PRESENT(is_D).OR.(.NOT.is_D)) THEN
      lf_exponent = 1.0_pr/(1.0_pr + lf_exponent)
    ELSE
      lf_exponent = 4*lf_der*lf_exponent/(1.0_pr + lf_exponent)**2
    END IF
  END FUNCTION lf_exponent

  FUNCTION lf_cubic_splines (enthalpy, is_D)
    IMPLICIT NONE
    REAL (pr),         INTENT(IN) :: enthalpy(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: is_D
    REAL (pr) :: lf_cubic_splines(SIZE(enthalpy))
    REAL (pr), DIMENSION(4) :: coeff
    REAL (pr) :: enthalpy_delta, halfpoint, delta, lf_der
    REAL (pr) :: enthalpy_Sm, enthalpy_Sp, enthalpy_Lm, enthalpy_Lp

    enthalpy_delta = fusion_delta + fusion_heat         ! change of enthalpy for fusion
    halfpoint = 1.0_pr + enthalpy_delta/2.0_pr          ! a middle (symmetric) point of fusion in terms of enthalpy
    delta = enthalpy_delta*smoothing_width              ! width of the cubic splines
    lf_der = 1.0_pr/enthalpy_delta                      ! first derivative of liquid fraction on enthalpy

    enthalpy_Sm = enthalpy_S - delta
    enthalpy_Sp = enthalpy_S + delta
    enthalpy_Lm = enthalpy_L - delta
    enthalpy_Lp = enthalpy_L + delta

    coeff = Spline_cubic(enthalpy_Sm, enthalpy_Sp, 0.0_pr, smoothing_width, 0.0_pr, lf_der)
    lf_cubic_splines = lf_piecewise(enthalpy, is_D)

    IF (.NOT.PRESENT(is_D).OR.(.NOT.is_D)) THEN
      WHERE (enthalpy.GT.enthalpy_Sm .AND. enthalpy.LE.enthalpy_Sp)
        lf_cubic_splines = coeff(2)*enthalpy**2 + coeff(3)*enthalpy + coeff(4)
      ELSEWHERE (enthalpy.GT.enthalpy_Lm .AND. enthalpy.LE.enthalpy_Lp)
        lf_cubic_splines = 1.0_pr - coeff(2)*(2*halfpoint - enthalpy)**2 - coeff(3)*(2*halfpoint-enthalpy) - coeff(4)
      END WHERE
    ELSE
      WHERE (enthalpy.GT.enthalpy_Sm .AND. enthalpy.LE.enthalpy_Sp)
        lf_cubic_splines = coeff(2)*enthalpy*2 + coeff(3)
      ELSEWHERE (enthalpy.GT.enthalpy_Lm .AND. enthalpy.LE.enthalpy_Lp)
        lf_cubic_splines = coeff(2)*(2*halfpoint - enthalpy)*2 + coeff(3)
      END WHERE
    END IF
  END FUNCTION lf_cubic_splines

  FUNCTION porosity (enthalpy)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy(:)
    REAL (pr) :: porosity(SIZE(enthalpy))
    porosity = MAX(0.0_pr, MIN(u(:,n_var_porosity), initial_porosity*(1.0_pr - liquid_fraction(enthalpy))))
  END FUNCTION porosity

  FUNCTION porosity_term ()
    IMPLICIT NONE
    REAL (pr) :: porosity_term(nwlt)
    porosity_term = (1.0_pr - u(:,n_var_porosity)) / (1.0_pr - initial_porosity)
  END FUNCTION porosity_term

  FUNCTION conductivity (temperature, liquid_fraction, is_D)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: temperature(:), liquid_fraction(:)
    REAL (pr) :: conductivity(SIZE(temperature))
    INTEGER, OPTIONAL, INTENT(IN) :: is_D

    conductivity = three_parameter_model(temperature, liquid_fraction, &
      Dconductivity_solid, Dconductivity_liquid, conductivity_fusion, is_D)
  END FUNCTION conductivity

  FUNCTION capacity (temperature, liquid_fraction, is_D)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: temperature(:), liquid_fraction(:)
    REAL (pr) :: capacity(SIZE(temperature))
    INTEGER, OPTIONAL, INTENT(IN) :: is_D

    capacity = three_parameter_model(temperature, liquid_fraction, &
      Dcapacity_solid, Dcapacity_liquid, capacity_fusion, is_D)
  END FUNCTION capacity

  ! is_D = 1 for partial derivative over \phi and 2 for derivative over T
  FUNCTION three_parameter_model (T, phi, Dsolid, Dliquid, fusion_jump, is_D)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) ::  T(:), phi(:)
    REAL (pr) :: three_parameter_model(SIZE(T))
    REAL (pr) :: Dsolid, Dliquid, fusion_jump
    INTEGER, OPTIONAL, INTENT(IN) :: is_D

    IF (.NOT.PRESENT(is_D)) THEN
      three_parameter_model = 1.0_pr + Dsolid*T + &
        (fusion_jump + (Dliquid - Dsolid)*(T - 1.0_pr))*phi
    ELSEIF (is_D .EQ. 1) THEN
      three_parameter_model = fusion_jump + (Dliquid - Dsolid)*(T - 1.0_pr)
    ELSEIF (is_D .EQ. 2) THEN
      three_parameter_model = Dsolid + (Dliquid - Dsolid)*phi
    END IF
  END FUNCTION three_parameter_model

  FUNCTION temperature (enthalpy)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy(:)
    REAL (pr) :: temperature(SIZE(enthalpy))
    REAL (pr) :: g1(SIZE(enthalpy)), g2(SIZE(enthalpy)), g3(SIZE(enthalpy)), phi(SIZE(enthalpy))

    phi = liquid_fraction(enthalpy)
    g2 = Dcapacity_solid*(1.0_pr - phi) + Dcapacity_liquid*phi
    WHERE (g2.GE.eps_zero)
      g1 = 1.0_pr*(1.0_pr - phi) + (1.0_pr + Dcapacity_solid - Dcapacity_liquid + capacity_fusion)*phi
      g3 = enthalpy_fusion(enthalpy) + ((Dcapacity_solid - Dcapacity_liquid)/2 + capacity_fusion)*phi
      temperature = (SQRT(g1**2 + 2*g2*g3) - g1)/g2
    ELSEWHERE
      temperature = enthalpy_fusion(enthalpy)
    END WHERE
  END FUNCTION temperature

  FUNCTION Dtemperature (enthalpy, temperature)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) ::  enthalpy(:), temperature(:)
    REAL (pr) :: Dtemperature(SIZE(enthalpy))
    REAL (pr) :: Dg2(SIZE(enthalpy)), Dg1(SIZE(enthalpy)), g1(SIZE(enthalpy)), g2(SIZE(enthalpy))
    REAL (pr) :: Dg3(SIZE(enthalpy)), phi(SIZE(enthalpy)), Dphi(SIZE(enthalpy))

    phi = liquid_fraction(enthalpy)
    g2 = Dcapacity_solid*(1.0_pr - phi) + Dcapacity_liquid*phi
    WHERE (g2.GE.eps_zero)
      Dphi = liquid_fraction(enthalpy, .TRUE.)
      g1 = 1.0_pr + (Dcapacity_solid - Dcapacity_liquid + capacity_fusion)*phi
      Dg1 = (Dcapacity_solid - Dcapacity_liquid + capacity_fusion)*Dphi
      Dg2 = (Dcapacity_liquid - Dcapacity_solid)*Dphi
      Dg3 = enthalpy_fusion(enthalpy, .TRUE.) + ((Dcapacity_solid - Dcapacity_liquid)/2 + capacity_fusion)*Dphi
      Dtemperature = (Dg3 - Dg1*temperature - Dg2*(temperature**2)/2) / (g2*temperature + g1)
    ELSEWHERE
      Dtemperature = enthalpy_fusion(enthalpy, .TRUE.)
    END WHERE
  END FUNCTION Dtemperature

  ! is_D = 0 for RHS, = 1 for DRHS and DRHS_diag
  FUNCTION enthalpy_fusion (enthalpy, is_D)
    IMPLICIT NONE
    REAL (pr),         INTENT(IN) :: enthalpy(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: is_D
    REAL (pr) :: enthalpy_fusion(SIZE(enthalpy))

    IF (.NOT.PRESENT(is_D).OR.(.NOT.is_D)) THEN
        enthalpy_fusion = enthalpy - fusion_heat*liquid_fraction(enthalpy)
    ELSE
        enthalpy_fusion = 1.0_pr - fusion_heat*liquid_fraction(enthalpy, .TRUE.)
    END IF
  END FUNCTION enthalpy_fusion

  FUNCTION Neumann_bc (enthalpy)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy(:)
    REAL (pr) :: Neumann_bc(SIZE(enthalpy))
    REAL (pr) :: T(SIZE(enthalpy)), phi(SIZE(enthalpy))

    T = temperature(enthalpy)
    phi = liquid_fraction(enthalpy)
    Neumann_bc = (1.0_pr - porosity(enthalpy)) * conductivity(T, phi) / &
      capacity(T, phi) * (1.0_pr - fusion_heat*liquid_fraction(enthalpy, .TRUE.))
  END FUNCTION Neumann_bc

  FUNCTION Spline_cubic (r_p, l_p, fr_p, fl_p, dfr_p, dfl_p)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: r_p, l_p, fr_p, fl_p, dfr_p, dfl_p
    REAL (pr), DIMENSION(4) :: Spline_cubic
    !Spline_cubic(1,2,3,4) are coefficient of polynomial of 3rd degree [ax^3+bx^2+cx+d]
      !where a = Spline_cubic(1) etc.
    Spline_cubic(1)= (2*fl_p - 2*fr_p - dfl_p*l_p - dfr_p*l_p + dfl_p*r_p + dfr_p*r_p)/(r_p-l_p)**3
    Spline_cubic(2) = (3*fr_p*l_p - 3*fl_p*l_p - 3*fl_p*r_p + 3*fr_p*r_p + dfl_p*l_p**2 &
      + 2*dfr_p*l_p**2 - 2*dfl_p*r_p**2 - dfr_p*r_p**2 + dfl_p*l_p*r_p - dfr_p*l_p*r_p)/(r_p-l_p)**3
    Spline_cubic(3) = (dfl_p*r_p**3 - dfr_p*l_p**3 + dfl_p*l_p*r_p**2 - 2*dfl_p*l_p**2*r_p &
      + 2*dfr_p*l_p*r_p**2 - dfr_p*l_p**2*r_p + 6*fl_p*l_p*r_p - 6*fr_p*l_p*r_p)/(r_p-l_p)**3
    Spline_cubic(4) = (r_p**3*(fl_p - dfl_p*l_p) + r_p*(dfr_p*l_p**3 + 3*fr_p*l_p**2) &
      - fr_p*l_p**3 - r_p**2*(3*fl_p*l_p - dfl_p*l_p**2 + dfr_p*l_p**2))/(r_p-l_p)**3
  END FUNCTION Spline_cubic

  FUNCTION Dirichlet_bc (enthalpy)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy(:)
    REAL (pr) :: Dirichlet_bc(SIZE(enthalpy))
    REAL (pr) :: T(SIZE(enthalpy)), DT(SIZE(enthalpy))

    T = temperature(enthalpy)
    DT = Dtemperature(enthalpy, T)
    Dirichlet_bc = F_heat_flux(T, .TRUE.)*DT
  END FUNCTION Dirichlet_bc

  FUNCTION F_heat_flux (temperature, is_D)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: temperature(:)
    REAL (pr) :: F_heat_flux(SIZE(temperature))
    LOGICAl, OPTIONAL, INTENT(IN) :: is_D

    IF (.NOT.PRESENT(is_D).OR.(.NOT.is_D)) THEN
      F_heat_flux = convective_transfer*T + emissivity*radiative_transfer*T**(dim+1)
    ELSE
      F_heat_flux = convective_transfer + emissivity*radiative_transfer*T**dim*(dim+1)
    END IF
  END FUNCTION F_heat_flux

  FUNCTION enthalpy (temperature, liquid_fraction)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: temperature(:), liquid_fraction(:)
    REAL (pr) :: enthalpy(SIZE(temperature))

    enthalpy = temperature + Dcapacity_solid/2*temperature**2 + liquid_fraction* &
      (fusion_heat + (Dcapacity_solid - Dcapacity_liquid + capacity_fusion)*(temperature - 1.0_pr) + &
      (Dcapacity_liquid - Dcapacity_solid)/2*(temperature**2 - 1.0_pr))
  END FUNCTION enthalpy
END MODULE user_case
