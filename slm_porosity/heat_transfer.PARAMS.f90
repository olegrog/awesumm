MODULE user_case_db

  !
  ! n_var_db must be set to
  !   Count( Union (n_var_adapt, n_var_interpolate))
  ! This is the total number of variable that either are used for adaptation
  ! criteria or need to be interpolated to new grid at each time step
  ! or that are needed for vector derivatives.
  !
  INTEGER , PARAMETER :: n_var_db  = 5
  INTEGER , PARAMETER :: n_ivar_db = 0
  INTEGER , PARAMETER :: max_dim   = 3

  !
  ! NOTE: db_tree and db_tree1 will not use the above parameters, although
  !       db_lines and db_tree_f will use them for database initialization
  !       and the program might be terminated for inappropriate parameters.
  !
  !       n_var_db - maximum number of REAL(pr) variables in the database
  !       n_ivar_db - maximum number of INTEGER variables in the database
  !       max_dim - maximum dimension in the database
  !

END MODULE user_case_db
