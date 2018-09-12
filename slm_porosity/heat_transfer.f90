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
  ! Style guide for the following user case:
  !   1) use long names for global variables (two or more words) and short names for local ones or of user types
  !   2) use ELEMENTAL and PURE function wherever possible
  !   3) avoid duplicating code as much as possible
  !   4) write 'IMPLICIT NONE' in all functions and subroutines
  !   5) avoid lines more than 120 symbols
  !   6) align column and comments
  !   7) prefer obvious names, otherwise add comments
  !   8) prefer intermediate local variables and functions to multiline commands
  !   9) use ' instead of "
  !

  !
  ! Case specific types
  !
  TYPE model3
    REAL(pr) :: Dsolid, Dliquid, jump
  END TYPE

  TYPE model5
    REAL(pr) :: Dsolid, D2solid, Dliquid, D2liquid, jump
  END TYPE

  ABSTRACT INTERFACE
    ! NB: fortran prohibits creating pointers to elemental functions (only to pure),
    !     but is possible to call a pure function from the elemental one
    !     and it will be considered as elemental. Fortran's magic :)
    PURE FUNCTION func_with_derivative(x, arg, is_D)
      USE precision
      IMPLICIT NONE
      REAL(pr), INTENT (IN) :: x, arg
      INTEGER,  INTENT (IN), OPTIONAL :: is_D
      REAL(pr) :: func_with_derivative
    END FUNCTION func_with_derivative
  END INTERFACE

  !
  ! Case specific variables
  !
  ! indices for variables
  INTEGER n_var_enthalpy
  INTEGER n_var_porosity
  INTEGER n_var_diffus
  INTEGER n_var_lfrac
  INTEGER n_var_temp
  INTEGER n_var_pressure

  ! quantities used in DRHS and algebraic BC
  REAL(pr), DIMENSION(:),   ALLOCATABLE :: diffusivity_prev, Ddiffusivity_prev, Dh_star_prev
  REAL(pr), DIMENSION(:,:), ALLOCATABLE :: grad_h_star_prev
  REAL(pr), DIMENSION(:),   ALLOCATABLE :: enthalpy_prev, temp_prev, Dtemp_prev, psi_prev
  INTEGER,  DIMENSION(:),   ALLOCATABLE :: i_p_face

  ! thermophysical properties
  TYPE(model3) :: conductivity3
  TYPE(model3) :: capacity3
  TYPE(model5) :: enthalpy5
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
  REAL(pr) :: initial_laser_position(3)
  REAL(pr) :: powder_depth

  ! numerics-specific parameters
  INTEGER  :: smoothing_method
  REAL(pr) :: fusion_smoothing
  REAL(pr) :: powder_smoothing
  REAL(pr) :: eps_zero
  REAL(pr) :: lfrac_scale
  REAL(pr) :: power_factor_2d
  REAL(pr) :: tol_newton
  INTEGER  :: max_iter_newton
  INTEGER  :: thickness_points
  LOGICAL  :: check_resolution

  ! derived quantities
  REAL(pr) :: enthalpy_S             ! temp=solidus,  phi=0
  REAL(pr) :: enthalpy_one           ! temp=1,        phi=1/2
  REAL(pr) :: enthalpy_L             ! temp=liquidus, phi=1
  REAL(pr) :: Dphi_one               ! maximum derivative of liquid fraction on enthalpy
  REAL(pr) :: interface_thickness    ! the solid--liquid interface thickness

  PROCEDURE(func_with_derivative), POINTER :: func_newton => null()

CONTAINS

  !
  ! The following variables must be setup in this routine:
  !
  ! n_integrated     - first n_integrated eqns will be acted on for time integration
  ! n_var_additional - interpolated variables (adapted and saved are automatically included)
  ! n_var
  !
  SUBROUTINE user_setup_pde (verb)
    USE variable_mapping
    IMPLICIT NONE
    LOGICAL, OPTIONAL :: verb
    LOGICAL, PARAMETER :: f_ = .FALSE., t_ = .TRUE.
    LOGICAL, DIMENSION(2), PARAMETER :: &
      ff = (/.FALSE.,.FALSE./), ft = (/.FALSE.,.TRUE./), tf = (/.TRUE.,.FALSE./), tt = (/.TRUE.,.TRUE./)
    INTEGER :: i

    ! Adapt mesh to
    ! 1) diffusivity as the driven nonlinearity,
    ! 2) liquid_fraction for the BC,
    ! 3) temperature for the case of constant thermophysical properties.
    CALL register_var('enthalpy',         integrated=t_, adapt=ff, saved=t_, interpolate=ff, exact=ff, req_restart=t_)
    CALL register_var('porosity',         integrated=f_, adapt=ff, saved=t_, interpolate=ff, exact=ff, req_restart=f_)
    CALL register_var('diffusivity',      integrated=f_, adapt=tt, saved=t_, interpolate=ff, exact=ff, req_restart=f_)
    CALL register_var('liquid_fraction',  integrated=f_, adapt=tt, saved=t_, interpolate=ff, exact=ff, req_restart=f_)
    CALL register_var('temperature',      integrated=f_, adapt=ft, saved=t_, interpolate=ff, exact=ff, req_restart=f_)

    CALL setup_mapping()
    CALL print_variable_registery(FULL=.TRUE.)

    n_var_enthalpy  = get_index('enthalpy')
    n_var_porosity  = get_index('porosity')
    n_var_diffus    = get_index('diffusivity')
    n_var_lfrac     = get_index('liquid_fraction')
    n_var_temp      = get_index('temperature')
    n_var_pressure  = n_var_enthalpy  ! for compatibility

    ALLOCATE(i_p_face(0:dim))
    i_p_face(0) = 1
    DO i = 1, dim
      i_p_face(i) = i_p_face(i-1)*3
    END DO

    ALLOCATE(Umn(n_var))
    Umn = 0.0_pr ! set up here if mean quantities are not zero and used in scales or equation
    scaleCoeff = 1.0_pr
    scaleCoeff(n_var_lfrac) = lfrac_scale

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
  ! u                        - array to fill in the exact solution
  ! nlocal                   - number of active wavelets
  ! ne_local                 - total number of equations
  ! t                        - time of current time step
  ! l_n_var_exact_soln_index - index into the elements of u for which we need to find the exact solution
  !
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
      ! 1. Introduce the initial melting pool in terms of the temperature field
      x_center = TRANSPOSE(SPREAD(laser_position(t), 2, nlocal))
      x_surface(:,:) = x(:,:)
      x_surface(:,dim) = x_center(:,dim)
      depth = x_center(:,dim) - x(:,dim)
      temp = initial_temp*EXP(-SUM((x_surface-x_center)**2, 2)/initial_pool_radius**2)
      phi = lf_from_temperature(temp)
      psi = porosity(SPREAD(initial_porosity, 1, nlocal), phi)
      k_0 = conductivity(temp, phi)
      lambda = laser_heat_flux(x_surface, nlocal, 1.0_pr - initial_pool_radius**-2) / &
        (initial_temp * k_0 * (1.0_pr - psi))
      temp = temp * EXP(-lambda*depth - (depth/initial_pool_radius)**2)

      ! 2. Calculate the enthalpy field from the temperature one
      !   a) find the approximate initial condition
      u(:,n_var_enthalpy) = enthalpy(temp, lf_from_temperature(temp))
      !   b) update it by solving the corresponding equation
      func_newton => enthalpy_equation
      u(:,n_var_enthalpy) = find_root(u(:,n_var_enthalpy), temp)
    END IF
  END SUBROUTINE user_initial_conditions

  !
  ! u_in     - fields on the adaptive grid
  ! nlocal   - number of active points
  ! ne_local - number of equations
  !
  SUBROUTINE user_algebraic_BC (Lu, u_in, nlocal, ne_local, jlev, meth)
    IMPLICIT NONE
    REAL(pr), INTENT(INOUT) :: Lu(nlocal*ne_local)
    REAL(pr), INTENT(IN) :: u_in(nlocal*ne_local)
    INTEGER,  INTENT(IN) :: nlocal, ne_local, jlev, meth
    INTEGER :: i, ie, shift, face_type, nloc
    REAL(pr), DIMENSION(ne_local,nlocal,dim) :: du, d2u
    INTEGER :: face(dim), iloc(nwlt)

    CALL c_diff_fast(u_in, du, d2u, jlev, nlocal, grad_meth(meth), 10, ne_local, 1, ne_local)

    DO ie = 1, ne_local
      shift = nlocal*(ie-1)
      DO face_type = 0, 3**dim - 1
        face = INT(MOD(face_type, i_p_face(1:dim))/i_p_face(0:dim-1))-1
        IF (ANY( face(1:dim) /= 0) .AND. ie == n_var_enthalpy) THEN
          CALL get_all_indices_by_face(face_type, jlev, nloc, iloc)
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

    CALL c_diff_diag(du, d2u, jlev, nlocal, grad_meth(meth), div_meth(meth), -10)

    DO ie = 1, ne_local
      shift = nlocal*(ie-1)
      DO face_type = 0, 3**dim - 1
        face = INT(MOD(face_type, i_p_face(1:dim))/i_p_face(0:dim-1))-1
        IF (ANY( face(1:dim) /= 0) .AND. ie == n_var_enthalpy) THEN
          CALL get_all_indices_by_face(face_type, jlev, nloc, iloc)
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
          CALL get_all_indices_by_face(face_type, jlev, nloc, iloc)
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

  !
  ! To make u divergence free
  !
  SUBROUTINE user_project (u, p, nlocal, meth)
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
    REAL(pr), DIMENSION(ng)         :: h, h_star, phi, psi
    REAL(pr), DIMENSION(ng,ne)      :: for_du
    REAL(pr), DIMENSION(ng,dim)     :: for_d2u
    REAL(pr), DIMENSION(ne,ng,dim)  :: du, du_dummy
    REAL(pr), DIMENSION(dim,ng,dim) :: d2u, d2u_dummy

    ie = n_var_enthalpy
    shift = ng*(ie-1)

    IF (IMEXswitch.LE.0) THEN
      user_rhs(shift+1:shift+ng) = 0.0_pr
    END IF

    IF (IMEXswitch.GE.0) THEN
      h = u_in(:,n_var_enthalpy)
      h_star = enthalpy_star(h)
      phi = liquid_fraction(h)
      psi = porosity(psi_prev, phi)
      for_du = RESHAPE(h_star, SHAPE(for_du))
      CALL c_diff_fast(for_du, du, du_dummy, j_lev, nwlt, grad_meth(meth), 10, ne, 1, ne)
      for_d2u = du(1,:,:) * SPREAD(diffusivity(h, psi), 2, dim)
      CALL c_diff_fast(for_d2u, d2u, d2u_dummy, j_lev, ng, div_meth(meth), 10, dim, 1, dim)
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
    REAL(pr), DIMENSION(ng)         :: pert_h_star
    REAL(pr), DIMENSION(ng,dim)     :: for_d2u, diff_pert_h_star
    REAL(pr), DIMENSION(ne,ng,dim)  :: du, du_dummy
    REAL(pr), DIMENSION(dim,ng,dim) :: d2u, d2u_dummy

    ie = n_var_enthalpy
    shift = ng*(ie-1)

    IF (IMEXswitch.LE.0) THEN
      user_Drhs(shift+1:shift+ng) = 0.0_pr
    END IF

    IF (IMEXswitch.GE.0) THEN
      pert_h_star = Dh_star_prev*pert_u(:,ie)
      CALL c_diff_fast(pert_h_star, du, du_dummy, j_lev, ng, grad_meth(meth), 10, ne, 1, ne)
      diff_pert_h_star = du(ie,:,:)
      for_d2u = &
        grad_h_star_prev * SPREAD(Ddiffusivity_prev * pert_u(:,ie), 2, dim) + &
        diff_pert_h_star * SPREAD(diffusivity_prev, 2, dim)
      CALL c_diff_fast(for_d2u, d2u, d2u_dummy, j_lev, ng, div_meth(meth), 10, dim, 1, dim)
      DO i = 1, dim
        user_Drhs(shift+1:shift+ng) = user_Drhs(shift+1:shift+ng) + d2u(i,:,i)
      END DO
    END IF
    IF (user_Drhs(1).NE.user_Drhs(1)) STOP '--- NaN in user_Drhs ---'
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
    REAL(pr), DIMENSION(ng,dim) :: du, d2u, for_d2u

    ie = n_var_enthalpy
    shift = ng*(ie-1)

    IF (IMEXswitch.LE.0) THEN
      user_Drhs_diag(shift+1:shift+ng) = 1.0_pr
    END IF

    IF (IMEXswitch.GE.0) THEN
      CALL c_diff_diag(du, d2u, j_lev, ng, grad_meth(meth), div_meth(meth), -11)
      for_d2u = &
        du * grad_h_star_prev * SPREAD(Ddiffusivity_prev, 2, dim) + &
        d2u * SPREAD(diffusivity_prev * Dh_star_prev, 2, dim)
      user_Drhs_diag(shift+1:shift+ng) = SUM(for_d2u, 2)
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
    CHARACTER(LEN=200) :: melting_pool, header, columns

    melting_pool = TRIM(res_path)//'melting_pool.txt'
    IF (startup_flag.EQ.0) THEN
      IF (par_rank.EQ.0) THEN
        header = '# time      max_temp    volume      length      depth'
        IF (dim.EQ.3) header = TRIM(header) // '       width'
        OPEN(555, FILE=melting_pool, STATUS='replace', ACTION='write')
        WRITE(555, *) TRIM(header)
        CLOSE(555)
      END IF
    ELSE
      CALL get_all_local_h(h_arr)
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
        columns = '(10(ES10.3, 2x))'
        vars = (/ t, max_temp, volume, length, depth /)
        IF (dim.EQ.3) vars(SIZE(vars)) = width
        OPEN(555, FILE=melting_pool, STATUS='old', POSITION='append', ACTION='write')
        WRITE(555, TRIM(columns)) vars
        CLOSE(555)
      END IF
    END IF
  END SUBROUTINE user_stats

  !
  ! Calculate drag and lift on obstacle using penalization formula
  !
  SUBROUTINE user_cal_force (u, n, t_local, force, drag, lift)
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
  ! Read input from "case_name".inp file
  ! case_name is in string file_name read from command line
  ! in read_command_line_input()
  !
  SUBROUTINE user_read_input ()
    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: fmt_real = '(A30,F10.7)'
    CHARACTER(LEN=*), PARAMETER :: stars = REPEAT('*', 20)
    REAL(pr) :: Denthalpy_L

    ! thermophysical properties
    call input_real('Dconductivity_solid', conductivity3%Dsolid, 'stop')
    call input_real('Dconductivity_liquid', conductivity3%Dliquid, 'stop')
    call input_real('conductivity_jump', conductivity3%jump, 'stop')
    call input_real('Dcapacity_solid', capacity3%Dsolid, 'stop')
    call input_real('Dcapacity_liquid', capacity3%Dliquid, 'stop')
    call input_real('capacity_jump', capacity3%jump, 'stop')
    call input_real('fusion_delta', fusion_delta, 'stop')
    call input_real('fusion_heat', fusion_heat, 'stop')
    Denthalpy_L = 1.0_pr + capacity3%Dsolid - capacity3%Dliquid + capacity3%jump
    enthalpy5 = model5(1.0_pr, capacity3%Dsolid, Denthalpy_L, capacity3%Dliquid, fusion_heat)

    ! bed parameters
    call input_real('laser_power', laser_power, 'stop')
    call input_real('scanning_speed', scanning_speed, 'stop')
    call input_real('initial_porosity', initial_porosity, 'stop')
    call input_real('convective_transfer', convective_transfer, 'stop')
    call input_real('radiative_transfer', radiative_transfer, 'stop')
    call input_real('absolute_temperature', absolute_temperature, 'stop')

    ! macroscopic model parameters
    call input_real('absorptivity', absorptivity, 'stop')
    call input_real('emissivity', emissivity, 'stop')

    ! initial conditions
    call input_real('initial_temp', initial_temp, 'stop')
    call input_real('initial_pool_radius', initial_pool_radius, 'stop')
    call input_real_vector('initial_laser_position', initial_laser_position, 3, 'stop')
    call input_real('powder_depth', powder_depth, 'stop')

    ! numerics-specific parameters
    call input_integer('smoothing_method', smoothing_method, 'stop')
    call input_real('fusion_smoothing', fusion_smoothing, 'stop')
    call input_real('powder_smoothing', powder_smoothing, 'stop')
    call input_real('eps_zero', eps_zero, 'stop')
    call input_real('lfrac_scale', lfrac_scale, 'stop')
    call input_real('power_factor_2d', power_factor_2d, 'stop')
    call input_real('tol_newton', tol_newton, 'stop')
    call input_integer('max_iter_newton', max_iter_newton, 'stop')
    call input_integer('thickness_points', thickness_points, 'stop')
    call input_logical('check_resolution', check_resolution, 'stop')

    ! calculate the real quantities
    enthalpy_S = enthalpy(1.0_pr - fusion_delta/2, 0.0_pr)
    enthalpy_L = enthalpy(1.0_pr + fusion_delta/2, 1.0_pr)
    interface_thickness = (enthalpy_L - enthalpy_S - fusion_heat)/capacity(1.0_pr, 0.5_pr)

    IF (par_rank.EQ.0) THEN
      PRINT *, stars, ' Real quantities ', stars
      PRINT fmt_real, 'enthalpy_S =', enthalpy_S
      PRINT fmt_real, 'enthalpy_L =', enthalpy_L
      PRINT fmt_real, 'interface thickness =', interface_thickness
    END IF

    ! smooth the solid--liquid interface
    fusion_delta = fusion_delta*fusion_smoothing

    ! calculate the derived quantities
    enthalpy_S = enthalpy(1.0_pr - fusion_delta/2, 0.0_pr)
    enthalpy_one = enthalpy(1.0_pr, 0.5_pr)
    enthalpy_L = enthalpy(1.0_pr + fusion_delta/2, 1.0_pr)
    Dphi_one = 1.0_pr/(enthalpy_L - enthalpy_S)
    interface_thickness = (enthalpy_L - enthalpy_S - fusion_heat)/capacity(1.0_pr, 0.5_pr)

    IF (par_rank.EQ.0) THEN
      PRINT *, stars, ' Computational quantities ', stars
      PRINT fmt_real, 'enthalpy_S =', enthalpy_S
      PRINT fmt_real, 'enthalpy_L =', enthalpy_L
      PRINT fmt_real, 'enthalpy(1) =', enthalpy_one
      PRINT fmt_real, 'Dphi(1) =', Dphi_one
      PRINT fmt_real, 'capacity(1) =', capacity(1.0_pr, 0.5_pr)
      PRINT fmt_real, 'conductivity(1) =', conductivity(1.0_pr, 0.5_pr)
      PRINT fmt_real, 'interface thickness =', interface_thickness
    END IF
  END SUBROUTINE user_read_input

  !
  ! Calculate any additional variables
  !
  ! flag: 0 - called while adapting to IC, 1 - called in main time integration loop
  !
  ! These additional variables are calculated and left in real space.
  !
  SUBROUTINE user_additional_vars (t_local, flag)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: t_local
    INTEGER,  INTENT(IN) :: flag
    REAL(pr) :: depth(nwlt), x_center(nwlt,dim)

    IF (flag.EQ.0) THEN
      ! the IC for porosity is calculated here since it is not the integrated variable
      x_center = TRANSPOSE(SPREAD(laser_position(t), 2, nwlt))
      depth = x_center(:,dim) - x(:,dim)
      u(:,n_var_porosity) = initial_porosity/(1.0_pr + EXP(-4*(powder_depth - depth)/powder_smoothing))
    END IF

    ! calculate the interpolated variables only
    u(:,n_var_temp) = temperature(u(:,n_var_enthalpy))
    u(:,n_var_lfrac) = liquid_fraction(u(:,n_var_enthalpy))
    u(:,n_var_porosity) = porosity(u(:,n_var_porosity), u(:,n_var_lfrac))
    u(:,n_var_diffus) = diffusivity(u(:,n_var_enthalpy), u(:,n_var_porosity))
  END SUBROUTINE user_additional_vars

  !
  ! Calculate any additional scalar variables
  !
  ! flag: 0 - called while adapting to IC, 1 - called in main time integration loop
  !
  SUBROUTINE user_scalar_vars (flag)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: flag
  END SUBROUTINE user_scalar_vars

  !
  ! Calculating Scales
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
    REAL(pr) :: h_arr(dim,nwlt), cfl_laser, cfl_diffusive

    use_default = .FALSE.
    CALL get_all_local_h(h_arr)
    cfl_laser = MAXVAL(dt/h_arr(1,:)*scanning_speed)
    cfl_diffusive = MAXVAL(dt/SUM(h_arr**2, 1)/diffusivity_prev) / (1.0_pr - fusion_heat*Dphi_one)
    CALL parallel_global_sum(REALMAXVAL=cfl_laser)
    CALL parallel_global_sum(REALMAXVAL=cfl_diffusive)
    cfl_out = cfl_laser
    IF (time_integration_method == TIME_INT_RK) cfl_out = MAX(cfl_out, cfl_diffusive)
    IF (par_rank.EQ.0) THEN
      PRINT *, 'CFL_laser =', cfl_laser, 'CFL_diffusive =', cfl_diffusive
    END IF
    IF (u(1,1).NE.u(1,1)) STOP '--- NaN in user_cal_cfl ---'
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

  FUNCTION grad_meth (meth)
    IMPLICIT NONE
    INTEGER :: meth, grad_meth
    grad_meth = meth + 2
  END FUNCTION grad_meth

  FUNCTION div_meth (meth)
    IMPLICIT NONE
    INTEGER :: meth, div_meth
    div_meth = meth + 4
  END FUNCTION div_meth

  SUBROUTINE user_pre_process
    IMPLICIT NONE
    INTEGER, PARAMETER :: meth = 1, nvar = 2
    REAL(pr), DIMENSION(nwlt) :: h, h_star, phi, psi, temp
    REAL(pr), DIMENSION(nwlt) :: Dphi, Dpsi, Dtemp
    REAL(pr), DIMENSION(nwlt,nvar)     :: for_du
    REAL(pr), DIMENSION(nvar,nwlt,dim) :: du, du_dummy
    REAL(pr), DIMENSION(dim) :: max_grad_temp, thickness, resolution

    ! 1. Use the calculated (from the enthalpy) values instead of the interpolated ones
    h = u(:,n_var_enthalpy)
    h_star = enthalpy_star(h)
    phi = liquid_fraction(h)
    psi = porosity(u(:,n_var_porosity), phi)
    temp = temperature(h)

    Dphi = liquid_fraction(h, 1)
    Dpsi = porosity(u(:,n_var_porosity), phi, Dphi)

    ! NB: 1) ng is not equal to nwlt after mesh adaptation
    !     2) ne is not equal to n_integrated until first call of time integration method
    for_du = RESHAPE((/ h_star, temp /), SHAPE(for_du))
    CALL c_diff_fast(for_du, du, du_dummy, j_lev, nwlt, grad_meth(meth), 10, nvar, 1, nvar)

    ! 2. Store the calculated variables as global variables with suffix "_prev"
    ! quantities used in DRHS
    CALL reallocate_scalar(diffusivity_prev); diffusivity_prev = diffusivity(h, psi)
    CALL reallocate_scalar(Ddiffusivity_prev); Ddiffusivity_prev = diffusivity(h, psi, Dpsi)
    CALL reallocate_scalar(Dh_star_prev); Dh_star_prev = enthalpy_star(h, 1)
    CALL reallocate_vector(grad_h_star_prev); grad_h_star_prev = du(1,:,:)

    ! quantities used in algebraic BC
    CALL reallocate_scalar(enthalpy_prev); enthalpy_prev = h
    CALL reallocate_scalar(temp_prev); temp_prev = temp
    CALL reallocate_scalar(Dtemp_prev); Dtemp_prev = temperature(h, temp)
    CALL reallocate_scalar(psi_prev); psi_prev = psi


    ! 3. Check if J_MX provides a sufficient resolution for solid--liquid interface thickness
    IF (check_resolution) THEN
      max_grad_temp = MAXVAL(ABS(du(2,:,:)), 1)
      CALL parallel_vector_sum(REALMAXVAL=max_grad_temp, LENGTH=dim)
      thickness = interface_thickness/max_grad_temp
      resolution = (xyzlimits(2,:) - xyzlimits(1,:)) / (mxyz*2**j_mx)
      IF (par_rank.EQ.0) THEN
        WRITE(*, '(A35, 3ES10.2)') 'The minimum interface thickness =', thickness
        WRITE(*, '(A35, 3ES10.2)') 'The minimum mesh step =', resolution
        WRITE(*, '(A35, 3ES10.2)') 'Their ratio =', thickness/resolution
      END IF
      IF (ANY(thickness.LT.resolution*thickness_points)) STOP '--- Insufficient resolution ---'
    END IF
  END SUBROUTINE user_pre_process

  SUBROUTINE user_post_process
    IMPLICIT NONE

    IF (ISNAN(SUM(u))) STOP '--- NaN in user_post_process ---'
    IF (scanning_speed*t.GE.(xyzlimits(1,1) + xyzlimits(2,1))) STOP '--- Finished ---'
  END SUBROUTINE user_post_process

  ELEMENTAL FUNCTION liquid_fraction (enthalpy, is_D)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: enthalpy
    INTEGER,  INTENT(IN), OPTIONAL :: is_D
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
    REAL(pr), INTENT(IN) :: enthalpy
    INTEGER,  INTENT(IN), OPTIONAL :: is_D
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
    REAL(pr), INTENT(IN) :: enthalpy
    INTEGER,  INTENT(IN), OPTIONAL :: is_D
    REAL(pr) :: lf_exponent

    lf_exponent = EXP(-4*Dphi_one*(enthalpy - enthalpy_one))

    IF (.NOT.PRESENT(is_D)) THEN
      lf_exponent = 1.0_pr/(1.0_pr + lf_exponent)
    ELSE
      lf_exponent = 4*Dphi_one*lf_exponent/(1.0_pr + lf_exponent)**2
    END IF
  END FUNCTION lf_exponent

  ! return Dporosity if the third argument is provided
  ELEMENTAL FUNCTION porosity (previous, phi, Dphi)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: previous, phi
    REAL(pr), INTENT(IN), OPTIONAL :: Dphi
    REAL(pr) :: porosity

    porosity = MAX(0.0_pr, MIN(previous, initial_porosity*(1.0_pr - phi)))
    IF (PRESENT(Dphi)) THEN
      IF (porosity.LT.previous) THEN
        porosity = -initial_porosity*Dphi
      ELSE
        porosity = 0.0_pr
      END IF
    END IF
  END FUNCTION porosity

  ELEMENTAL FUNCTION diffusivity (enthalpy, psi, Dpsi)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: enthalpy, psi
    REAL(pr), INTENT(IN), OPTIONAL :: Dpsi
    REAL(pr) :: diffusivity, phi, temp, k_0, c_p, pterm
    REAL(pr) :: Dphi, Dtemp, Dk_0, Dc_p, Dpterm

    phi = liquid_fraction(enthalpy)
    temp = temperature(enthalpy)
    k_0 = conductivity(temp, phi)
    c_p = capacity(temp, phi)

    Dphi = liquid_fraction(enthalpy, 1)
    Dtemp = temperature(enthalpy, temp)
    Dk_0 = conductivity(temp, phi, 1)*Dtemp + conductivity(temp, phi, 2)*Dphi
    Dc_p = capacity(temp, phi, 1)*Dtemp + capacity(temp, phi, 2)*Dphi
    pterm = (1.0_pr - psi) / (1.0_pr - initial_porosity)

    IF (.NOT.PRESENT(Dpsi)) THEN
      diffusivity = k_0 / c_p * pterm
    ELSE
      Dpterm = -Dpsi / (1.0_pr - initial_porosity)
      diffusivity = (Dk_0*c_p - Dc_p*k_0) / c_p**2 * pterm + k_0 / c_p * Dpterm
    END IF
  END FUNCTION diffusivity

  ELEMENTAL FUNCTION conductivity (temp, phi, is_D)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: temp, phi
    INTEGER,  INTENT(IN), OPTIONAL :: is_D
    REAL(pr) :: conductivity

    conductivity = three_parameter_model(temp, phi, conductivity3, is_D)
  END FUNCTION conductivity

  ELEMENTAL FUNCTION capacity (temp, phi, is_D)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: temp, phi
    INTEGER,  INTENT(IN), OPTIONAL :: is_D
    REAL(pr) :: capacity

    capacity = three_parameter_model(temp, phi, capacity3, is_D)
  END FUNCTION capacity

  ! if is_D present, return the partial derivative with respect to: 1 - temperature, 2 - liquid_fraction
  ELEMENTAL FUNCTION three_parameter_model (temp, phi, m, is_D)
    IMPLICIT NONE
    REAL(pr),     INTENT(IN) :: temp, phi
    TYPE(model3), INTENT(IN) :: m
    INTEGER,      INTENT(IN), OPTIONAL :: is_D
    REAL(pr) :: three_parameter_model, delta

    delta = m%Dliquid - m%Dsolid
    IF (.NOT.PRESENT(is_D)) THEN
      three_parameter_model = 1.0_pr + m%Dsolid*temp + (m%jump + delta*(temp - 1.0_pr))*phi
    ELSE IF (is_D.EQ.1) THEN
      three_parameter_model = m%Dsolid + delta*phi
    ELSE IF (is_D.EQ.2) THEN
      three_parameter_model = m%jump + delta*(temp - 1.0_pr)
    END IF
  END FUNCTION three_parameter_model

  ! if is_D present, return the partial derivative with respect to: 1 - temperature, 2 - liquid_fraction
  ELEMENTAL FUNCTION five_parameter_model (temp, phi, m, is_D)
    IMPLICIT NONE
    REAL(pr),     INTENT(IN) :: temp, phi
    TYPE(model5), INTENT(IN) :: m
    INTEGER,      INTENT(IN), OPTIONAL :: is_D
    REAL(pr) :: five_parameter_model, delta, delta2, before_phi

    delta = m%Dliquid - m%Dsolid
    delta2 = m%D2liquid - m%D2solid
    before_phi = m%jump + delta*(temp - 1.0_pr) + delta2*(temp**2 - 1.0_pr)/2

    IF (.NOT.PRESENT(is_D)) THEN
      five_parameter_model = m%Dsolid*temp + m%D2solid*temp**2/2 + before_phi*phi
    ELSE IF (is_D.EQ.1) THEN
      five_parameter_model = m%Dsolid + m%D2solid*temp + (delta + delta2*temp)*phi
    ELSE IF (is_D.EQ.2) THEN
      five_parameter_model = before_phi
    END IF
  END FUNCTION five_parameter_model

  ! return Dtemperature if the second argument is provided
  ELEMENTAL FUNCTION temperature (enthalpy, temp)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: enthalpy
    REAL(pr), INTENT(IN), OPTIONAL :: temp
    REAL(pr) :: temperature, h_star, phi, g1, g2, g3
    REAL(pr) :: Dh_star, Dphi, Dg1, Dg2, Dg3
    REAL(pr) :: delta, delta2

    h_star = enthalpy_star(enthalpy)
    phi = liquid_fraction(enthalpy)
    delta = enthalpy5%Dliquid - enthalpy5%Dsolid; delta2 = enthalpy5%D2liquid - enthalpy5%D2solid
    g1 = enthalpy5%Dsolid + delta*phi; g2 = enthalpy5%D2solid + delta2*phi; g3 = h_star + (delta + delta2/2)*phi

    IF (.NOT.PRESENT(temp)) THEN
      IF (g2.GE.eps_zero) THEN
        temperature = (SQRT(g1**2 + 2*g2*g3) - g1)/g2
      ELSE
        temperature = h_star
      END IF
    ELSE
      Dh_star = enthalpy_star(enthalpy, 1)
      Dphi = liquid_fraction(enthalpy, 1)
      Dg1 = delta*Dphi; Dg2 = delta2*Dphi; Dg3 = Dh_star + (delta + delta2/2)*Dphi
      IF (g2.GE.eps_zero) THEN
        temperature = (Dg3 - Dg1*temp - Dg2*(temp**2)/2) / (g2*temp + g1)
      ELSE
        temperature = Dh_star
      END IF
    END IF
  END FUNCTION temperature

  ! if is_D present, return the partial derivative with respect to: 1 - temperature, 2 - liquid_fraction
  ELEMENTAL FUNCTION enthalpy (temp, phi, is_D)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: temp, phi
    INTEGER,  INTENT(IN), OPTIONAL :: is_D
    REAL(pr) :: enthalpy

    enthalpy = five_parameter_model(temp, phi, enthalpy5, is_D)
  END FUNCTION enthalpy

  PURE FUNCTION enthalpy_equation (h, temp, is_D)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: h, temp
    INTEGER,  INTENT(IN), OPTIONAL :: is_D
    REAL(pr) :: enthalpy_equation, phi, Dphi

    phi = liquid_fraction(h)
    Dphi = liquid_fraction(h, 1)
    IF (.NOT.PRESENT(is_D)) THEN
      enthalpy_equation = h - enthalpy(temp, phi)
    ELSE
      enthalpy_equation = 1.0_pr - enthalpy(temp, phi, 2)*Dphi
    END IF
  END FUNCTION enthalpy_equation

  ELEMENTAL FUNCTION enthalpy_star (enthalpy, is_D)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: enthalpy
    INTEGER,  INTENT(IN), OPTIONAL :: is_D
    REAL(pr) :: enthalpy_star

    IF (.NOT.PRESENT(is_D)) THEN
      enthalpy_star = enthalpy - fusion_heat*liquid_fraction(enthalpy)
    ELSE
      enthalpy_star = 1.0_pr - fusion_heat*liquid_fraction(enthalpy, 1)
    END IF
  END FUNCTION enthalpy_star

  ELEMENTAL FUNCTION Neumann_bc (enthalpy, psi)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: enthalpy, psi
    REAL(pr) :: Neumann_bc, temp, phi, Dh_star

    temp = temperature(enthalpy)
    phi = liquid_fraction(enthalpy)
    Dh_star = enthalpy_star(enthalpy, 1)
    Neumann_bc = (1.0_pr - psi) * conductivity(temp, phi) / capacity(temp, phi) * Dh_star
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
    REAL(pr), INTENT(IN) :: temp
    INTEGER,  INTENT(IN), OPTIONAL :: is_D
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
    IF (COUNT(mask).EQ.0 .OR. ABS(x2-x1).GT.min_distance) domain_bound = -HUGE(REAL(pr))
    CALL parallel_global_sum(REALMAXVAL=domain_bound)
    domain_bound = sgn*domain_bound
  END FUNCTION domain_bound

  !
  ! Find root of equation using the Newton algorithm
  !
  ! NB: Dummy procedures are not allowed as arguments in elemental functions,
  !     but can be transfered to them as pointers to pure functions!
  !
  ! func_newton - global pointer to function with its derivative
  ! x0          - initial guess
  ! arg         - additional argument for func_newton
  ELEMENTAL FUNCTION find_root (x0, arg) RESULT(x)
    IMPLICIT NONE
    REAL(pr), INTENT(IN) :: x0, arg
    REAL(pr) :: x, f
    INTEGER  :: i
    i = 0
    x = x0
    f = func_newton(x0, arg)
    DO WHILE(ABS(f).GT.tol_newton .AND. i.LT.max_iter_newton)
      x = x - f/func_newton(x, arg, 1)
      f = func_newton(x, arg)
      i = i + 1
    END DO
  END FUNCTION find_root

END MODULE user_case
