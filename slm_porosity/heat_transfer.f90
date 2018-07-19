MODULE user_case
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
  INTEGER n_var_enthalpy
  INTEGER n_var_porosity
  INTEGER n_var_h_star
  INTEGER n_var_lfrac
  INTEGER n_var_temp
  INTEGER n_var_pressure

  ! quantities used in RHS, DRHS, and algebraic BC
  REAL(pr), DIMENSION(:), ALLOCATABLE :: diffusivity_prev, Ddiffusivity_prev, Dh_star_prev
  REAL(pr), DIMENSION(:,:), ALLOCATABLE :: diff_h_star_prev
  REAL(pr), DIMENSION(:), ALLOCATABLE :: enthalpy_prev, temp_prev, Dtemp_prev, psi_prev
  INTEGER, DIMENSION(:), ALLOCATABLE :: i_p_face

  ! thermophysical properties
  REAL(pr) :: Dconductivity_solid
  REAL(pr) :: Dconductivity_liquid
  REAL(pr) :: conductivity_fusion
  REAL(pr) :: Dcapacity_solid
  REAL(pr) :: Dcapacity_liquid
  REAL(pr) :: capacity_fusion
  REAL(pr) :: fusion_delta
  REAL(pr) :: fusion_heat

  ! bed parameters
  REAL(pr) :: laser_power
  REAL(pr) :: scanning_speed
  REAL(pr) :: initial_porosity
  REAL(pr) :: convective_transfer
  REAL(pr) :: radiative_transfer
  REAL(pr) :: absolute_temperature

  ! macroscopic model parameters
  REAL(pr) :: absorptivity
  REAL(pr) :: emissivity

  ! initial conditions
  REAL(pr) :: initial_temp
  REAL(pr) :: initial_pool_radius
  REAL(pr), DIMENSION(3) :: initial_laser_position

  ! numerics-specific parameters
  INTEGER smoothing_method
  REAL(pr) :: smoothing_factor
  REAL(pr) :: eps_zero
  REAL(pr) :: porosity_scale
  REAL(pr) :: power_factor_2d

  ! derived quantities
  REAL(pr) :: enthalpy_S             ! temp=solidus,  phi=0
  REAL(pr) :: enthalpy_one           ! temp=1,        phi=1/2
  REAL(pr) :: enthalpy_L             ! temp=liquidus, phi=1
  REAL(pr) :: Dphi_one               ! maximum derivative of liquid fraction on enthalpy
CONTAINS

  !
  ! The following variables must be setup in this routine:
  !
  ! n_integrated     ! first n_integrated eqns will be acted on for time integration
  ! n_var_additional ! interpolated variables (adapted and saved are automatically included)
  ! n_var
  !
  SUBROUTINE user_setup_pde (verb)
    USE variable_mapping
    IMPLICIT NONE
    LOGICAL, OPTIONAL :: verb
    LOGICAL, PARAMETER :: f = .FALSE., t = .TRUE.
    LOGICAL, DIMENSION(2), PARAMETER :: &
      ff = (/.FALSE.,.FALSE./), ft = (/.FALSE.,.TRUE./), tf = (/.TRUE.,.FALSE./), tt = (/.TRUE.,.TRUE./)
    INTEGER :: i

    CALL register_var('enthalpy',         integrated=t, adapt=ff, saved=t, interpolate=ff, exact=ff, req_restart=t)
    CALL register_var('porosity',         integrated=f, adapt=tt, saved=t, interpolate=ff, exact=ff, req_restart=f)
    CALL register_var('enthalpy_star',    integrated=f, adapt=tt, saved=f, interpolate=ff, exact=ff, req_restart=f)
    CALL register_var('liquid_fraction',  integrated=f, adapt=ff, saved=t, interpolate=ff, exact=ff, req_restart=f)
    CALL register_var('temperature',      integrated=f, adapt=ff, saved=t, interpolate=ff, exact=ff, req_restart=f)
    CALL register_var('pressure',         integrated=f, adapt=ff, saved=f, interpolate=ff, exact=ff, req_restart=f)

    CALL setup_mapping()
    CALL print_variable_registery(FULL=.TRUE.)

    n_var_enthalpy  = get_index('enthalpy')
    n_var_porosity  = get_index('porosity')
    n_var_h_star    = get_index('enthalpy_star')
    n_var_lfrac     = get_index('liquid_fraction')
    n_var_temp      = get_index('temperature')
    n_var_pressure  = get_index('pressure')

    ALLOCATE(i_p_face(0:dim))
    i_p_face(0) = 1
    DO i = 1, dim
      i_p_face(i) = i_p_face(i-1)*3
    END DO

    ALLOCATE(Umn(n_var))
    Umn = 0.0_pr !set up here if mean quantities are not zero and used in scales or equation
    scaleCoeff = 1.
    scaleCoeff(n_var_porosity) = porosity_scale

    IF (verb_level.GT.0) THEN
      PRINT *, 'n_integrated = ', n_integrated
      PRINT *, 'n_var = ', n_var
      PRINT *, 'n_var_exact = ', n_var_exact
      PRINT *, '*******************Variable Names*******************'
      DO i = 1, n_var
        WRITE (*, u_variable_names_fmt) u_variable_names(i)
      END DO
      PRINT *, '****************************************************'
    END IF
  END SUBROUTINE user_setup_pde

  !
  ! Set the exact solution for comparison to the simulated solution
  !
  ! u          - array to fill in the exact solution
  ! nlocal       - number of active wavelets
  ! ne_local        - total number of equations
  ! t          - time of current time step
  ! l_n_var_exact_soln_index - index into the elements of u for which we need to
  !                            find the exact solution
  SUBROUTINE user_exact_soln (u, nlocal, t_local, l_n_var_exact_soln)
    IMPLICIT NONE
    REAL(pr), INTENT(INOUT) :: u(nlocal,n_var_exact)
    INTEGER,  INTENT(IN) :: nlocal
    REAL(pr), INTENT(IN) :: t_local
    LOGICAL,  INTENT(IN) :: l_n_var_exact_soln(n_var)
  END SUBROUTINE user_exact_soln

  SUBROUTINE user_initial_conditions (u, nlocal, ne_local, t_local, scl, scl_fltwt, iter)
    IMPLICIT NONE
    REAL(pr), INTENT(INOUT) :: u(nlocal,ne_local)
    INTEGER,  INTENT(IN) :: nlocal, ne_local
    REAL(pr), INTENT(IN) :: t_local, scl(1:n_var), scl_fltwt
    INTEGER,  INTENT(INOUT) :: iter       ! iteration of call while adapting initial grid
    REAL(pr), DIMENSION(nlocal) :: depth, lambda, temp, phi, psi, k_0
    REAL(pr) :: x_center(nlocal,dim)      ! coordinates of the center of the laser beam
    REAL(pr) :: x_surface(nlocal,dim)     ! coordinates of the surface of laser energy absorption

    IF ( IC_restart_mode.EQ.0) THEN
      x_center = TRANSPOSE(SPREAD(laser_position(t), 2, nlocal))
      x_surface(:,:) = x(:,:)
      x_surface(:,dim) = x_center(:,dim)
      depth = x_center(:,dim) - x(:,dim)
      temp = initial_temp*EXP(-SUM((x_surface-x_center)**2, 2)/initial_pool_radius**2)
      phi = lf_from_temperature(temp)
      psi = porosity(phi, SPREAD(initial_porosity, 1, nlocal))
      k_0 = conductivity(temp, phi)
      lambda = laser_heat_flux(x_surface, nlocal, 1.0_pr - initial_pool_radius**-2) / &
        (initial_temp * k_0 * (1.0_pr - psi))
      u(:,n_var_temp) = temp * EXP(-lambda*depth - (depth/initial_pool_radius)**2)
    END IF
  END SUBROUTINE user_initial_conditions

  !
  ! u_in      - fields on the adaptive grid
  ! nlocal    - number of active points
  ! ne_local  - number of equations
  !
  SUBROUTINE user_algebraic_BC (Lu, u_in, nlocal, ne_local, jlev, meth)
    IMPLICIT NONE
    REAL(pr), INTENT(INOUT) :: Lu(nlocal*ne_local)
    REAL(pr), INTENT(IN) :: u_in(nlocal*ne_local)
    INTEGER,  INTENT(IN) :: nlocal, ne_local, jlev, meth
    INTEGER :: i, ie, shift, face_type, nloc
    REAL(pr), DIMENSION(ne_local,nlocal,dim) :: du, d2u
    INTEGER :: face(dim), iloc(nwlt)

    CALL c_diff_fast (u_in, du, d2u, jlev, nlocal, meth, 10, ne_local, 1, ne_local)

    DO ie = 1, ne_local
      shift = nlocal*(ie-1)
      DO face_type = 0, 3**dim - 1
        face = INT(MOD(face_type, i_p_face(1:dim))/i_p_face(0:dim-1))-1
        IF (ANY( face(1:dim) /= 0) .AND. ie == n_var_enthalpy) THEN
          CALL get_all_indices_by_face (face_type, jlev, nloc, iloc)
          IF (nloc > 0) THEN
            IF (face(dim) > 0) THEN
              Lu(shift+iloc(1:nloc)) = &
                Neumann_bc(enthalpy_prev(iloc(1:nloc)), psi_prev(iloc(1:nloc))) * du(ie, iloc(1:nloc), dim) + &
                Dirichlet_bc(enthalpy_prev(iloc(1:nloc))) * u_in(shift+iloc(1:nloc))
            ELSEIF (dim == 3 .AND. face(2) < 0) THEN
              Lu(shift+iloc(1:nloc)) = du(ie, iloc(1:nloc), 2)
            ELSE
              Lu(shift+iloc(1:nloc)) = u_in(shift+iloc(1:nloc))
            END IF
          END IF
        END IF
      END DO
    END DO
  END SUBROUTINE user_algebraic_BC

  SUBROUTINE user_algebraic_BC_diag (Lu_diag, nlocal, ne_local, jlev, meth)
    IMPLICIT NONE
    REAL(pr), INTENT(INOUT) :: Lu_diag(nlocal*ne_local)
    INTEGER,  INTENT(IN) :: nlocal, ne_local, jlev, meth
    INTEGER :: i, ie, shift, face_type, nloc
    REAL(pr), DIMENSION(nlocal,dim) :: du, d2u
    INTEGER :: face(dim), iloc(nwlt)

    CALL c_diff_diag (du, d2u, jlev, nlocal, meth, meth, 10)

    DO ie = 1, ne_local
      shift = nlocal*(ie-1)
      DO face_type = 0, 3**dim - 1
        face = INT(MOD(face_type, i_p_face(1:dim))/i_p_face(0:dim-1))-1
        IF (ANY( face(1:dim) /= 0) .AND. ie == n_var_enthalpy) THEN
          CALL get_all_indices_by_face (face_type, jlev, nloc, iloc)
          IF (nloc > 0) THEN
            IF (face(dim) > 0) THEN
              Lu_diag(shift+iloc(1:nloc)) = &
                Neumann_bc(enthalpy_prev(iloc(1:nloc)), psi_prev(iloc(1:nloc))) * du(iloc(1:nloc), dim) + &
                Dirichlet_bc(enthalpy_prev(iloc(1:nloc)))
            ELSEIF (dim == 3 .AND. face(2) < 0) THEN
              Lu_diag(shift+iloc(1:nloc)) = du(iloc(1:nloc), 2)
            ELSE
              Lu_diag(shift+iloc(1:nloc)) = 1.0_pr
            END IF
          END IF
        END IF
      END DO
    END DO
  END SUBROUTINE user_algebraic_BC_diag

  SUBROUTINE user_algebraic_BC_rhs (rhs, ne_local, nlocal, jlev)
    IMPLICIT NONE
    REAL(pr), INTENT(INOUT) :: rhs(nlocal*ne_local)
    INTEGER,  INTENT(IN) :: ne_local, nlocal, jlev
    INTEGER :: i, ie, shift, face_type, nloc
    INTEGER :: face(dim), iloc(nwlt)

    DO ie = 1, ne_local
      shift = nlocal*(ie-1)
      DO face_type = 0, 3**dim - 1
        face = INT(MOD(face_type, i_p_face(1:dim))/i_p_face(0:dim-1))-1
        IF (ANY( face(1:dim) /= 0) .AND. ie == n_var_enthalpy) THEN
          CALL get_all_indices_by_face (face_type, jlev, nloc, iloc)
          IF (nloc > 0) THEN
            IF (face(dim) > 0) THEN
              rhs(shift+iloc(1:nloc)) = &
                laser_heat_flux(x(iloc(1:nloc),:), nloc) - other_heat_flux(temp_prev(iloc(1:nloc))) + &
                other_heat_flux(temp_prev(iloc(1:nloc)), 1) * Dtemp_prev(iloc(1:nloc))*enthalpy_prev(iloc(1:nloc))
            ELSE
              rhs(shift+iloc(1:nloc)) = 0
            END IF
          END IF
        END IF
      END DO
    END DO
  END SUBROUTINE user_algebraic_BC_rhs

  SUBROUTINE user_project (u, p, nlocal, meth)
    !--Makes u divergence free
    IMPLICIT NONE
    REAL(pr), INTENT(INOUT) :: u(nlocal,n_integrated), p(nlocal)
    INTEGER,  INTENT(IN) :: nlocal, meth
  END SUBROUTINE user_project

  FUNCTION user_rhs (u_in, p)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: u_in(ng,ne)
    REAL(pr), INTENT(IN) :: p(ng)
    REAL(pr) :: user_rhs(n)
    INTEGER :: ie, shift, i
    INTEGER, PARAMETER :: meth = 1
    REAL(pr), DIMENSION(ng,dim) :: for_du
    REAL(pr), DIMENSION(dim,ng,dim):: d2u, d2u_dummy

    ie = n_var_enthalpy
    shift = ng*(ie-1)

    IF (IMEXswitch.LE.0) THEN
      user_rhs(shift+1:shift+ng) = 0.0_pr
    END IF

    IF (IMEXswitch.GE.0) THEN
      for_du = diff_h_star_prev * SPREAD(diffusivity_prev, 2, dim)
      CALL c_diff_fast(for_du, d2u, d2u_dummy, j_lev, ng, meth, 10, dim, 1, dim)
      DO i = 1, dim
        user_rhs(shift+1:shift+ng) = user_rhs(shift+1:shift+ng) + d2u(i,:,i)
      END DO
    END IF
  END FUNCTION user_rhs

  FUNCTION user_Drhs (pert_u, u_prev, meth)
    IMPLICIT NONE
    REAL(pr), INTENT(IN), DIMENSION(ng,ne) :: pert_u, u_prev
    INTEGER,  INTENT(IN) :: meth
    REAL(pr) :: user_Drhs(n)
    INTEGER :: ie, shift, i
    INTEGER, SAVE :: k = 0
    REAL(pr), DIMENSION(ng) :: pert_h_star
    REAL(pr), DIMENSION(ng,dim) :: for_du, diff_pert_h_star
    REAL(pr), DIMENSION(ne,ng,dim) :: du, du_dummy
    REAL(pr), DIMENSION(dim,ng,dim) :: d2u, d2u_dummy

    ie = n_var_enthalpy
    shift = ng*(ie-1)

    IF (IMEXswitch.LE.0) THEN
      user_Drhs(shift+1:shift+ng) = 0.0_pr
    END IF

    IF (IMEXswitch.GE.0) THEN
      pert_h_star = Dh_star_prev*pert_u(:,ie)
      CALL c_diff_fast(pert_h_star, du, du_dummy, j_lev, ng, meth, 10, ne, 1, ne)
      diff_pert_h_star = du(ie,:,:)
      for_du = &
        diff_h_star_prev * SPREAD(Ddiffusivity_prev * pert_h_star, 2, dim) + &
        diff_pert_h_star * SPREAD(diffusivity_prev, 2, dim)
      CALL c_diff_fast(for_du, d2u, d2u_dummy, j_lev, ng, meth, 10, dim, 1, dim)
      DO i = 1, dim
        user_Drhs(shift+1:shift+ng) = user_Drhs(shift+1:shift+ng) + d2u(i,:,i)
      END DO
    END IF
    IF (user_Drhs(1).NE.user_Drhs(1)) THEN
      PRINT *, '--- NaN in user_Drhs ---'
      CALL ABORT
    END IF
    IF (BTEST(verb_level,1)) THEN
      PRINT *, 'Drhs', k
      k = k + 1
    END IF
  END FUNCTION user_Drhs

  FUNCTION user_Drhs_diag (meth)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: meth
    REAL(pr) :: user_Drhs_diag(n)
    INTEGER :: ie, shift, i
    REAL(pr), DIMENSION(ng,dim) :: du, d2u, for_du

    ie = n_var_enthalpy
    shift = ng*(ie-1)

    IF (IMEXswitch.LE.0) THEN
      user_Drhs_diag(shift+1:shift+ng) = 1.0_pr
    END IF

    IF (IMEXswitch.GE.0) THEN
      CALL c_diff_diag(du, d2u, j_lev, ng, meth, meth, -11)
      for_du = &
        du * diff_h_star_prev * SPREAD(Ddiffusivity_prev, 2, dim) + &
        d2u * SPREAD(diffusivity_prev * Dh_star_prev, 2, dim)
      user_Drhs_diag(shift+1:shift+ng) = SUM(for_du, 2)
    END IF
  END FUNCTION user_Drhs_diag

  FUNCTION user_chi (nlocal, t_local)
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: nlocal
    REAL(pr), INTENT(IN) :: t_local
    REAL(pr), DIMENSION(nlocal) :: user_chi

    user_chi = 0.0_pr
  END FUNCTION user_chi


  FUNCTION user_mapping (xlocal, nlocal, t_local)
    USE curvilinear
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: xlocal(nlocal,dim)
    INTEGER,  INTENT(IN) :: nlocal
    REAL(pr), INTENT(IN) :: t_local
    REAL(pr) :: user_mapping(nlocal,dim)

    user_mapping(:,1:dim) = xlocal(:,1:dim)
  END FUNCTION user_mapping

  !
  ! Calculate any statitics
  !
  ! startup_flag: 0 - when adapting to IC, 1 - in the main integration loop
  !
  SUBROUTINE user_stats (u, j_mn, startup_flag)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: u(nwlt,n_var)
    INTEGER,  INTENT(IN) :: startup_flag, j_mn
    REAL(pr) :: max_temp               ! maximum of temperature
    REAL(pr) :: volume                 ! volume of the melting pool
    REAL(pr) :: width, length, depth   ! dimensions of the melting pool
    REAL(pr) :: phi(nwlt), h_arr(dim,nwlt), vars(dim+3)
    REAL(pr) :: xmin, xmax, ymax, zmin, h_max
    REAL(pr), PARAMETER :: threshold = 0.5_pr
    CHARACTER(LEN=16), PARAMETER :: file_name = 'melting_pool.txt'
    CHARACTER(LEN=100) :: header, columns

    IF (startup_flag.EQ.0) THEN
      IF (par_rank.EQ.0) THEN
        header = '# time      max_temp    volume      length      depth'
        IF (dim.EQ.3) header = TRIM(header) // '       width'
        OPEN(555, FILE=file_name, STATUS='replace', ACTION='write')
        WRITE(555, *) TRIM(header)
        CLOSE(555)
      END IF
    ELSE
      CALL get_all_local_h (h_arr)
      h_max = MAXVAL(h_arr)
      phi = u(:,n_var_lfrac)
      volume = SUM(phi*dA)
      max_temp = MAXVAL(temp_prev)
      CALL parallel_global_sum(REAL=volume)
      CALL parallel_global_sum(REALMAXVAL=max_temp)

      xmin = domain_bound(phi, threshold, 1,   -1, h_max)
      xmax = domain_bound(phi, threshold, 1,    1, h_max)
      ymax = domain_bound(phi, threshold, 2,    1, h_max)
      zmin = domain_bound(phi, threshold, dim, -1, h_max)

      length = xmax - xmin
      width = 2*ymax
      depth = xyzlimits(2,dim) - zmin

      IF (par_rank.EQ.0) THEN
        columns = '(' // REPEAT('ES10.3, 2x', SIZE(vars)) // ')'
        vars = (/ t, max_temp, volume, length, depth /)
        IF (dim.EQ.3) vars(SIZE(vars)) = width
        OPEN(555, FILE=file_name, STATUS='old', POSITION='append', ACTION='write')
        WRITE(555, TRIM(columns)) vars
        CLOSE(555)
      END IF
    END IF
  END SUBROUTINE user_stats

  SUBROUTINE user_cal_force (u, n, t_local, force, drag, lift)
    !--Calculates drag and lift on obstacle using penalization formula
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: u(n,dim)
    INTEGER,  INTENT(IN) :: n
    REAL(pr), INTENT(IN) :: t_local
    REAL(pr), INTENT(INOUT) :: force(dim)
    REAL(pr), INTENT(OUT) :: drag, lift
    drag = 0.0_pr
    lift = 0.0_pr
    !
    ! There is no obstacle in flow
    !
  END SUBROUTINE user_cal_force

  !
  ! Read input from "case_name"_pde.inp file
  ! case_name is in string file_name read from command line
  ! in read_command_line_input()
  !
  SUBROUTINE user_read_input ()
    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: fmt_real = '(A30,F10.7)'
    CHARACTER(LEN=*), PARAMETER :: stars = REPEAT('*', 20)
    REAL(pr) :: thickness

    ! thermophysical properties
    call input_real ('Dconductivity_solid', Dconductivity_solid, 'stop')
    call input_real ('Dconductivity_liquid', Dconductivity_liquid, 'stop')
    call input_real ('conductivity_fusion', conductivity_fusion, 'stop')
    call input_real ('Dcapacity_solid', Dcapacity_solid, 'stop')
    call input_real ('Dcapacity_liquid', Dcapacity_liquid, 'stop')
    call input_real ('capacity_fusion', capacity_fusion, 'stop')
    call input_real ('fusion_delta', fusion_delta, 'stop')
    call input_real ('fusion_heat', fusion_heat, 'stop')

    ! bed parameters
    call input_real ('laser_power', laser_power, 'stop')
    call input_real ('scanning_speed', scanning_speed, 'stop')
    call input_real ('initial_porosity', initial_porosity, 'stop')
    call input_real ('convective_transfer', convective_transfer, 'stop')
    call input_real ('radiative_transfer', radiative_transfer, 'stop')
    call input_real ('absolute_temperature', absolute_temperature, 'stop')

    ! macroscopic model parameters
    call input_real ('absorptivity', absorptivity, 'stop')
    call input_real ('emissivity', emissivity, 'stop')

    ! initial conditions
    call input_real ('initial_temp', initial_temp, 'stop')
    call input_real ('initial_pool_radius', initial_pool_radius, 'stop')
    call input_real_vector ('initial_laser_position', initial_laser_position, 3, 'stop')

    ! numerics-specific parameters
    call input_integer ('smoothing_method', smoothing_method, 'stop')
    call input_real ('smoothing_factor', smoothing_factor, 'stop')
    call input_real ('eps_zero', eps_zero, 'stop')
    call input_real ('porosity_scale', porosity_scale, 'stop')
    call input_real ('power_factor_2d', power_factor_2d, 'stop')

    ! calculate the real quantities
    enthalpy_S = enthalpy(1.0_pr - fusion_delta/2, 0.0_pr)
    enthalpy_L = enthalpy(1.0_pr + fusion_delta/2, 1.0_pr)
    thickness = (enthalpy_L - enthalpy_S - fusion_heat)/capacity(1.0_pr, 0.5_pr)

    IF (par_rank.EQ.0) THEN
      PRINT *, stars, ' Real quantities ', stars
      PRINT fmt_real, 'enthalpy_S =', enthalpy_S
      PRINT fmt_real, 'enthalpy_L =', enthalpy_L
      PRINT fmt_real, 'interface thickness =', thickness
    END IF

    ! smooth solid-liquid interface
    fusion_delta = fusion_delta*smoothing_factor

    ! calculate the derived quantities
    enthalpy_S = enthalpy(1.0_pr - fusion_delta/2, 0.0_pr)
    enthalpy_one = enthalpy(1.0_pr, 0.5_pr)
    enthalpy_L = enthalpy(1.0_pr + fusion_delta/2, 1.0_pr)
    Dphi_one = 1.0_pr/(enthalpy_L - enthalpy_S)
    thickness = (enthalpy_L - enthalpy_S - fusion_heat)/capacity(1.0_pr, 0.5_pr)

    IF (par_rank.EQ.0) THEN
      PRINT *, stars, ' Computational quantities ', stars
      PRINT fmt_real, 'enthalpy_S =', enthalpy_S
      PRINT fmt_real, 'enthalpy_L =', enthalpy_L
      PRINT fmt_real, 'enthalpy(1) =', enthalpy_one
      PRINT fmt_real, 'Dphi(1) =', Dphi_one
      PRINT fmt_real, 'capacity(1) =', capacity(1.0_pr, 0.5_pr)
      PRINT fmt_real, 'conductivity(1) =', conductivity(1.0_pr, 0.5_pr)
      PRINT fmt_real, 'interface thickness =', thickness
    END IF
  END SUBROUTINE user_read_input

  !
  ! calculate any additional variables
  !
  ! arg
  ! flag - 0 calledwhile adapting to IC, 1 - called in main time integration loop
  !
  ! These additional variables are calculated and left in real space.
  !
  SUBROUTINE user_additional_vars (t_local, flag)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: t_local
    INTEGER,  INTENT(IN) :: flag ! 0- called during adaption to IC, 1 called during main integration loop

    IF (.NOT.flag) THEN
      u(:,n_var_porosity) = initial_porosity
      u(:,n_var_enthalpy) = enthalpy(u(:,n_var_temp), lf_from_temperature(u(:,n_var_temp)))
    END IF
    ! calculate the interpolated variables only
    u(:,n_var_temp) = temperature(u(:,n_var_enthalpy))
    u(:,n_var_lfrac) = liquid_fraction(u(:,n_var_enthalpy))
    u(:,n_var_porosity) = porosity(u(:,n_var_lfrac), u(:,n_var_porosity))
    u(:,n_var_h_star) = u(:,n_var_enthalpy) - fusion_heat*u(:,n_var_lfrac)
  END SUBROUTINE user_additional_vars

  !
  ! calculate any additional scalar variables
  !
  ! arg
  ! flag - 0 calledwhile adapting to IC, 1 - called in main time integration loop
  !
  SUBROUTINE user_scalar_vars (flag)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: flag ! 0- called during adaption to IC, 1 called during main integration loop
  END SUBROUTINE user_scalar_vars

  !
  !************ Calculating Scales ***************************
  !
  ! Note the order of the components in the scl array
  ! correspond to u_tn, v_tn, w_tn, u_tn-1, v_tn-1, w_tn-1, u_tn-2, v_tn-2, w_tn-2
  !
  SUBROUTINE user_scales (flag, use_default, u, nlocal, ne_local, l_n_var_adapt, l_n_var_adapt_index, scl, scl_fltwt)
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: flag ! 0 during initial adaptation, then 1 during time advancement
    LOGICAL,  INTENT(INOUT) :: use_default
    REAL(pr), INTENT(IN) :: u(nlocal,ne_local)
    INTEGER,  INTENT(IN) :: nlocal, ne_local
    LOGICAL,  INTENT(IN) :: l_n_var_adapt(ne_local)
    INTEGER,  INTENT(IN) :: l_n_var_adapt_index(1:ne_local)
    REAL(pr), INTENT(INOUT) :: scl(ne_local)
    REAL(pr), INTENT(IN) :: scl_fltwt !weight for temporal filter on scl
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
    LOGICAL,  INTENT(INOUT) :: use_default
    REAL(pr), INTENT(IN)    :: u(nwlt,n_integrated)
    REAL(pr), INTENT(INOUT) :: cfl_out
    REAL(pr) :: h_arr(dim,nwlt)

    use_default = .FALSE.
    CALL get_all_local_h (h_arr)
    cfl_out = MAXVAL(dt/h_arr(1,:)*scanning_speed)
    CALL parallel_global_sum(REALMAXVAL=cfl_out)
    IF (u(1,1).NE.u(1,1)) THEN
      PRINT *, '--- INFINITE ---'
      CALL ABORT
    END IF
  END SUBROUTINE user_cal_cfl

  SUBROUTINE user_init_sgs_model
    IMPLICIT NONE
  END SUBROUTINE user_init_sgs_model

  SUBROUTINE user_sgs_force (u_loc, nlocal)
    IMPLICIT NONE
    REAL(pr), INTENT(INOUT) :: u_loc(nlocal,n_integrated)
    INTEGER,  INTENT(IN) :: nlocal
  END SUBROUTINE user_sgs_force

  FUNCTION user_sound_speed (u, neq, nwlt)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: u(nwlt,neq)
    INTEGER,  INTENT(IN) :: nwlt, neq
    REAL(pr) :: user_sound_speed(nwlt)

    user_sound_speed(:) = 0.0_pr
  END FUNCTION user_sound_speed

  SUBROUTINE reallocate_scalar (var)
    IMPLICIT NONE
    REAL(pr), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: var
    IF (ALLOCATED(var)) DEALLOCATE(var)
    ALLOCATE(var(nwlt))
  END SUBROUTINE reallocate_scalar

  SUBROUTINE reallocate_vector (var)
    IMPLICIT NONE
    REAL(pr), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: var
    IF (ALLOCATED(var)) DEALLOCATE(var)
    ALLOCATE(var(nwlt,dim))
  END SUBROUTINE reallocate_vector

  SUBROUTINE user_pre_process
    IMPLICIT NONE
    REAL(pr), DIMENSION(nwlt) :: h_star, phi, psi, temp, k_0, c_p
    REAL(pr), DIMENSION(nwlt) :: Dphi, Dtemp, Dk_0, Dc_p
    REAL(pr), DIMENSION(nwlt,ne) :: for_du
    REAL(pr), DIMENSION(ne,nwlt,dim) :: du, du_dummy
    INTEGER, PARAMETER :: meth = 1

    h_star = enthalpy_star(u(:,n_var_enthalpy))
    phi = liquid_fraction(u(:,n_var_enthalpy))
    psi = porosity(phi, u(:,n_var_porosity))
    temp = temperature(u(:,n_var_enthalpy))
    k_0 = conductivity(temp, phi)
    c_p = capacity(temp, phi)

    Dphi = liquid_fraction(u(:,n_var_enthalpy), 1)
    Dtemp = temperature(u(:,n_var_enthalpy), temp)
    Dk_0 = conductivity(temp, phi, 1)*Dphi + conductivity(temp, phi, 2)*Dtemp
    Dc_p = capacity(temp, phi, 1)*Dphi + capacity(temp, phi, 2)*Dtemp

    ! quantities used in RHS and DRHS
    CALL reallocate_scalar(diffusivity_prev); diffusivity_prev = k_0 / c_p * porosity_term(psi)
    CALL reallocate_scalar(Ddiffusivity_prev); Ddiffusivity_prev = (Dk_0*c_p - Dc_p*k_0) / c_p**2 * porosity_term(psi)
    CALL reallocate_scalar(Dh_star_prev); Dh_star_prev = enthalpy_star(u(:,n_var_enthalpy), 1)
    ! quantities used in algebraic BC
    CALL reallocate_scalar(enthalpy_prev); enthalpy_prev = u(:,n_var_enthalpy)
    CALL reallocate_scalar(temp_prev); temp_prev = temp
    CALL reallocate_scalar(Dtemp_prev); Dtemp_prev = Dtemp
    CALL reallocate_scalar(psi_prev); psi_prev = psi

    ! NB: here, ng is not equal to nwlt
    for_du = RESHAPE((/ h_star /), SHAPE(for_du))
    CALL c_diff_fast(for_du, du, du_dummy, j_lev, nwlt, meth, 10, ne, 1, ne)
    CALL reallocate_vector(diff_h_star_prev); diff_h_star_prev = du(1,:,:)
  END SUBROUTINE user_pre_process

  SUBROUTINE user_post_process
    IMPLICIT NONE
  END SUBROUTINE user_post_process

  ELEMENTAL FUNCTION liquid_fraction (enthalpy, is_D)
    IMPLICIT NONE
    REAL(pr),          INTENT(IN) :: enthalpy
    INTEGER, OPTIONAL, INTENT(IN) :: is_D
    REAL(pr) :: liquid_fraction

    IF (smoothing_method.EQ.0) THEN         ! C^0
      liquid_fraction = lf_piecewise(enthalpy, is_D)
    ELSE IF (smoothing_method.EQ.1) THEN    ! C^\infty
      liquid_fraction = lf_exponent(enthalpy, is_D)
    END IF
  END FUNCTION liquid_fraction

  ELEMENTAL FUNCTION lf_from_temperature (temp)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: temp
    REAL(pr) :: lf_from_temperature

    lf_from_temperature = 1.0_pr/(1.0_pr + EXP(-4*(temp - 1.0_pr)/fusion_delta))
  END FUNCTION lf_from_temperature

  ELEMENTAL FUNCTION lf_piecewise (enthalpy, is_D)
    IMPLICIT NONE
    REAL(pr),          INTENT(IN) :: enthalpy
    INTEGER, OPTIONAL, INTENT(IN) :: is_D
    REAL(pr) :: lf_piecewise

    IF (.NOT.PRESENT(is_D)) THEN
      lf_piecewise = MAX(0.0_pr, MIN(1.0_pr, (enthalpy - enthalpy_S)*Dphi_one))
    ELSE
      IF (enthalpy_S.LT.enthalpy .AND. enthalpy.LT.enthalpy_L) THEN
        lf_piecewise = Dphi_one
      ELSE
        lf_piecewise = 0.0_pr
      END IF
    END IF
  END FUNCTION lf_piecewise

  ELEMENTAL FUNCTION lf_exponent (enthalpy, is_D)
    IMPLICIT NONE
    REAL(pr),          INTENT(IN) :: enthalpy
    INTEGER, OPTIONAL, INTENT(IN) :: is_D
    REAL(pr) :: lf_exponent

    lf_exponent = EXP(-4*Dphi_one*(enthalpy - enthalpy_one))

    IF (.NOT.PRESENT(is_D)) THEN
      lf_exponent = 1.0_pr/(1.0_pr + lf_exponent)
    ELSE
      lf_exponent = 4*Dphi_one*lf_exponent/(1.0_pr + lf_exponent)**2
    END IF
  END FUNCTION lf_exponent

  ELEMENTAL FUNCTION porosity (liquid_fraction, previous)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: liquid_fraction, previous
    REAL(pr) :: porosity
    porosity = MAX(0.0_pr, MIN(previous, initial_porosity*(1.0_pr - liquid_fraction)))
  END FUNCTION porosity

  ELEMENTAL FUNCTION porosity_term (porosity)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: porosity
    REAL(pr) :: porosity_term
    porosity_term = (1.0_pr - porosity) / (1.0_pr - initial_porosity)
  END FUNCTION porosity_term

  ELEMENTAL FUNCTION conductivity (temp, phi, is_D)
    IMPLICIT NONE
    REAL(pr),          INTENT(IN) :: temp, phi
    INTEGER, OPTIONAL, INTENT(IN) :: is_D
    REAL(pr) :: conductivity

    conductivity = three_parameter_model(temp, phi, &
      Dconductivity_solid, Dconductivity_liquid, conductivity_fusion, is_D)
  END FUNCTION conductivity

  ELEMENTAL FUNCTION capacity (temp, phi, is_D)
    IMPLICIT NONE
    REAL(pr),          INTENT(IN) :: temp, phi
    INTEGER, OPTIONAL, INTENT(IN) :: is_D
    REAL(pr) :: capacity

    capacity = three_parameter_model(temp, phi, &
      Dcapacity_solid, Dcapacity_liquid, capacity_fusion, is_D)
  END FUNCTION capacity

  ! is_D: 1 - partial derivative of liquid_fraction, 2 - partial derivative of temperature
  ELEMENTAL FUNCTION three_parameter_model (temp, phi, Dsolid, Dliquid, jump, is_D)
    IMPLICIT NONE
    REAL(pr),          INTENT(IN) :: temp, phi, Dsolid, Dliquid, jump
    INTEGER, OPTIONAL, INTENT(IN) :: is_D
    REAL(pr) :: three_parameter_model

    IF (.NOT.PRESENT(is_D)) THEN
      three_parameter_model = 1.0_pr + Dsolid*temp + (jump + (Dliquid - Dsolid)*(temp - 1.0_pr))*phi
    ELSE IF (is_D.EQ.1) THEN
      three_parameter_model = jump + (Dliquid - Dsolid)*(temp - 1.0_pr)
    ELSE IF (is_D.EQ.2) THEN
      three_parameter_model = Dsolid + (Dliquid - Dsolid)*phi
    END IF
  END FUNCTION three_parameter_model

  ! return Dtemperature if the second argument is provided
  ELEMENTAL FUNCTION temperature (enthalpy, temp)
    IMPLICIT NONE
    REAL(pr),           INTENT(IN) :: enthalpy
    REAL(pr), OPTIONAL, INTENT(IN) :: temp
    REAL(pr) :: temperature, h_star, phi, g1, g2, g3
    REAL(pr) :: Dh_star, Dphi, Dg1, Dg2, Dg3

    h_star = enthalpy_star(enthalpy)
    phi = liquid_fraction(enthalpy)
    g1 = 1.0_pr + (Dcapacity_solid - Dcapacity_liquid + capacity_fusion)*phi
    g2 = Dcapacity_solid*(1.0_pr - phi) + Dcapacity_liquid*phi
    g3 = h_star + ((Dcapacity_solid - Dcapacity_liquid)/2 + capacity_fusion)*phi

    IF (.NOT.PRESENT(temp)) THEN
      IF (g2.GE.eps_zero) THEN
        temperature = (SQRT(g1**2 + 2*g2*g3) - g1)/g2
      ELSE
        temperature = h_star
      END IF
    ELSE
      Dh_star = enthalpy_star(enthalpy, 1)
      Dphi = liquid_fraction(enthalpy, 1)
      Dg1 = (Dcapacity_solid - Dcapacity_liquid + capacity_fusion)*Dphi
      Dg2 = (Dcapacity_liquid - Dcapacity_solid)*Dphi
      Dg3 = Dh_star + ((Dcapacity_solid - Dcapacity_liquid)/2 + capacity_fusion)*Dphi
      IF (g2.GE.eps_zero) THEN
        temperature = (Dg3 - Dg1*temp - Dg2*(temp**2)/2) / (g2*temp + g1)
      ELSE
        temperature = Dh_star
      END IF
    END IF
  END FUNCTION temperature

  ELEMENTAL FUNCTION enthalpy (temp, phi)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: temp, phi
    REAL(pr) :: enthalpy

    enthalpy = temp + Dcapacity_solid/2*temp**2 + phi*(fusion_heat + &
      (Dcapacity_solid - Dcapacity_liquid + capacity_fusion)*(temp - 1.0_pr) + &
      (Dcapacity_liquid - Dcapacity_solid)/2*(temp**2 - 1.0_pr))
  END FUNCTION enthalpy

  ELEMENTAL FUNCTION enthalpy_star (enthalpy, is_D)
    IMPLICIT NONE
    REAL(pr),          INTENT(IN) :: enthalpy
    INTEGER, OPTIONAL, INTENT(IN) :: is_D
    REAL(pr) :: enthalpy_star

    IF (.NOT.PRESENT(is_D)) THEN
      enthalpy_star = enthalpy - fusion_heat*liquid_fraction(enthalpy)
    ELSE
      enthalpy_star = 1.0_pr - fusion_heat*liquid_fraction(enthalpy, 1)
    END IF
  END FUNCTION enthalpy_star

  ELEMENTAL FUNCTION Neumann_bc (enthalpy, porosity)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: enthalpy, porosity
    REAL(pr) :: Neumann_bc, temp, phi, Dh_star

    temp = temperature(enthalpy)
    phi = liquid_fraction(enthalpy)
    Dh_star = enthalpy_star(enthalpy, 1)
    Neumann_bc = (1.0_pr - porosity) * conductivity(temp, phi) / capacity(temp, phi) * Dh_star
  END FUNCTION Neumann_bc

  ELEMENTAL FUNCTION Dirichlet_bc (enthalpy)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: enthalpy
    REAL(pr) :: Dirichlet_bc, temp, Dtemp

    temp = temperature(enthalpy)
    Dtemp = temperature(enthalpy, temp)
    Dirichlet_bc = other_heat_flux(temp, 1)*Dtemp
  END FUNCTION Dirichlet_bc

  ELEMENTAL FUNCTION other_heat_flux (temp, is_D)
    IMPLICIT NONE
    REAL(pr),          INTENT(IN) :: temp
    INTEGER, OPTIONAL, INTENT(IN) :: is_D
    REAL(pr) :: other_heat_flux

    IF (.NOT.PRESENT(is_D)) THEN
      other_heat_flux = convective_transfer*temp + &
        emissivity*radiative_transfer*(temp + absolute_temperature)**(dim+1)
    ELSE
      other_heat_flux = convective_transfer + &
        emissivity*radiative_transfer*(temp + absolute_temperature)**dim*(dim+1)
    END IF
  END FUNCTION other_heat_flux

  FUNCTION laser_heat_flux (x, nlocal, factor)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: x(nlocal,dim)
    INTEGER,  INTENT(IN) :: nlocal
    REAL(pr), INTENT(IN), OPTIONAL :: factor
    REAL(pr) :: laser_heat_flux(nlocal), x_center(nlocal,dim), factor_

    x_center = TRANSPOSE(SPREAD(laser_position(t), 2, nlocal))
    factor_ = 1.0_pr; IF (PRESENT(factor)) factor_ = factor
    laser_heat_flux = absorptivity*laser_power/pi**(.5*(dim-1))*EXP(-SUM(factor_*(x-x_center)**2, 2))
    IF (dim.EQ.2) THEN
      laser_heat_flux = laser_heat_flux*power_factor_2d
    END IF
  END FUNCTION laser_heat_flux

  FUNCTION laser_position (time)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: time
    REAL(pr) :: laser_position(dim)
    ! NB: (:) cannot be removed
    laser_position(:) = initial_laser_position(:dim)
    laser_position(dim) = xyzlimits(2, dim)
    laser_position(1) = laser_position(1) + scanning_speed*time
  END FUNCTION laser_position

  ELEMENTAL FUNCTION linear_interpolation (x, x1, x2, y1, y2)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: x, x1, x2, y1, y2
    REAL(pr) :: linear_interpolation
    IF (ABS(x2-x1).LE.eps_zero) THEN
      linear_interpolation = (y1 + y2)/2
    ELSE
      linear_interpolation = (y2 - y1)/(x2 - x1)*(x - x1) + y1
    END IF
  END FUNCTION linear_interpolation

  FUNCTION domain_bound (field, threshold, axis, sgn, min_distance)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: field(nwlt), threshold, min_distance
    INTEGER,  INTENT(IN) :: axis, sgn
    REAL(pr) :: domain_bound, x1, x2
    INTEGER x1_loc(1), x2_loc(1), i
    LOGICAL, DIMENSION(nwlt) :: domain, mask

    domain = field.GE.threshold
    x1 = MAXVAL(sgn*x(:,axis), MASK=domain)
    x1_loc = MAXLOC(sgn*x(:,axis), MASK=domain)
    mask = .NOT.domain .AND. sgn*x(:,axis).GT.sgn*x(x1_loc(1),axis)
    DO i = 1, dim
      IF (i.NE.axis) mask = mask .AND. x(:,i).EQ.x(x1_loc(1),i)
    END DO
    x2 = MINVAL(sgn*x(:,axis), MASK=mask)
    x2_loc = MINLOC(sgn*x(:,axis), MASK=mask)
    domain_bound = linear_interpolation(threshold, field(x1_loc(1)), field(x2_loc(1)), x1, x2)
    IF (COUNT(mask).EQ.0 .OR. ABS(x2-x1).GT.min_distance) domain_bound = -1.0_pr/0.0_pr
    CALL parallel_global_sum(REALMAXVAL=domain_bound)
    domain_bound = sgn*domain_bound
  END FUNCTION domain_bound

END MODULE user_case
