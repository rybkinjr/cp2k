!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE optimize_dfet_mat_potential
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_maxabsval,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_multiply,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_trace,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_add, &
                                              dbcsr_dot, &
                                              dbcsr_copy, &
                                              dbcsr_set, &
                                              dbcsr_release_p
   USE embed_types,                     ONLY: opt_dfet_mat_pot_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_dmfet_potential'

   PUBLIC :: prepare_dfet_mat_opt, release_dfet_mat_opt,  check_dfet_mat, dfet_mat_opt_step, &
             release_dfet_mat_term, dfet_mat_ks_matrix, add_dfet_mat_to_w_func

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param opt_dmfet ...
!> \param opt_dmfet_section ...
! **************************************************************************************************
   SUBROUTINE read_opt_dfet_mat_section(opt_dfet_mat, opt_dfet_mat_section)
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      TYPE(section_vals_type), POINTER                   :: opt_dfet_mat_section

      ! Read keywords

      CALL section_vals_val_get(opt_dfet_mat_section, "N_ITER", &
                                i_val=opt_dfet_mat%n_iter)

      CALL section_vals_val_get(opt_dfet_mat_section, "TRUST_RAD", &
                                r_val=opt_dfet_mat%trust_rad)

   END SUBROUTINE read_opt_dfet_mat_section

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param opt_dmfet ...
!> \param opt_dmfet_section ...
! **************************************************************************************************
   SUBROUTINE prepare_dfet_mat_opt(qs_env, opt_dfet_mat, opt_dfet_mat_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      TYPE(section_vals_type), POINTER                   :: opt_dfet_mat_section
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      INTEGER                                            :: nspins, ispin, diff_size

      ! Read the input
      CALL read_opt_dfet_mat_section(opt_dfet_mat, opt_dfet_mat_section)

      nspins = 1
      IF (opt_dfet_mat%open_shell_embed) THEN
         nspins = 2
      ENDIF

      ! Get the orbital coefficients
      CALL get_qs_env(qs_env, matrix_s=matrix_s)

      NULLIFY (opt_dfet_mat%p_total)
      CALL dbcsr_allocate_matrix_set(opt_dfet_mat%p_diff, nspins)
      CALL dbcsr_allocate_matrix_set(opt_dfet_mat%v_mat, nspins)
      DO ispin = 1, nspins
         ALLOCATE (opt_dfet_mat%p_diff(ispin)%matrix)
         CALL dbcsr_init_p(opt_dfet_mat%p_diff(ispin)%matrix)
         CALL dbcsr_create(opt_dfet_mat%p_diff(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(opt_dfet_mat%p_diff(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(opt_dfet_mat%p_diff(ispin)%matrix, 0.0_dp)
         !
         ALLOCATE (opt_dfet_mat%v_mat(ispin)%matrix)
         CALL dbcsr_init_p(opt_dfet_mat%v_mat(ispin)%matrix)
         CALL dbcsr_create(opt_dfet_mat%v_mat(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(opt_dfet_mat%v_mat(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(opt_dfet_mat%v_mat(ispin)%matrix, 0.0_dp)
         !
      ENDDO

      ! Array to store functional values
      ALLOCATE (opt_dfet_mat%w_func(opt_dfet_mat%n_iter))
      opt_dfet_mat%w_func = 0.0_dp

      ! Allocate max_diff and int_diff
      diff_size = 1
      IF (opt_dfet_mat%open_shell_embed) diff_size = 2
      ALLOCATE (opt_dfet_mat%max_diff(diff_size))
      ALLOCATE (opt_dfet_mat%int_diff(diff_size))

   END SUBROUTINE prepare_dfet_mat_opt

! **************************************************************************************************
!> \brief ...
!> \param opt_dmfet ...
! **************************************************************************************************
   SUBROUTINE release_dfet_mat_opt(opt_dfet_mat)
      TYPE(opt_dfet_mat_pot_type)                           :: opt_dfet_mat
      INTEGER                                            :: ispin, nspins

      nspins = size(opt_dfet_mat%p_diff)

      DO ispin = 1, nspins
         CALL dbcsr_release_p(opt_dfet_mat%p_diff(ispin)%matrix)
         CALL dbcsr_release_p(opt_dfet_mat%v_mat(ispin)%matrix)
      ENDDO
      DEALLOCATE (opt_dfet_mat%p_diff)
      DEALLOCATE (opt_dfet_mat%v_mat)

   END SUBROUTINE release_dfet_mat_opt

! **************************************************************************************************
!> \brief ...
!> \param opt_dmfet ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE dfet_mat_opt_step(opt_dfet_mat)
      TYPE(opt_dfet_mat_pot_type)                           :: opt_dfet_mat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dfet_opt_step', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins

      CALL timeset(routineN, handle)
      nspins = 1
      IF (opt_dfet_mat%open_shell_embed) nspins = 2

      DO ispin = 1, nspins
         CALL dbcsr_add(opt_dfet_mat%v_mat(ispin)%matrix, opt_dfet_mat%p_diff(ispin)%matrix,&
               alpha_scalar=1.0_dp, beta_scalar=-opt_dfet_mat%trust_rad)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE dfet_mat_opt_step

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param opt_dmfet ...
!> \param i_force_eval ...
! **************************************************************************************************
   SUBROUTINE dfet_mat_ks_matrix(qs_env, opt_dfet_mat, i_force_eval)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      INTEGER                                            :: i_force_eval

      INTEGER                                            :: ispin, nspins
      LOGICAL                                            :: my_change_spin
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      my_change_spin = .FALSE.
      IF ((opt_dfet_mat%change_spin) .AND. (i_force_eval .EQ. 2)) my_change_spin = .TRUE.

      ! Get Kohn-Sham environment, Kohn-Sham matrix, overlap ...
      NULLIFY (ks_env)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, dft_control=dft_control, ks_env=ks_env)
      nspins = dft_control%nspins

      DO ispin = 1, nspins
         ALLOCATE (opt_dfet_mat%dfet_mat_term(ispin)%matrix)
         CALL dbcsr_init_p(opt_dfet_mat%dfet_mat_term(ispin)%matrix)
         CALL dbcsr_create(opt_dfet_mat%dfet_mat_term(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(opt_dfet_mat%dfet_mat_term(ispin)%matrix, 0.0_dp)
      ENDDO

      DO ispin = 1, nspins

         IF (.NOT. my_change_spin) THEN
            CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(ispin)%matrix, opt_dfet_mat%v_mat(ispin)%matrix)
         ELSE
            IF (ispin .EQ. 1) THEN
               CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(ispin)%matrix, opt_dfet_mat%v_mat(2)%matrix)
            ENDIF
            IF (ispin .EQ. 2) THEN
               CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(ispin)%matrix, opt_dfet_mat%v_mat(1)%matrix)
            ENDIF
         ENDIF

         ! Special case of a closed-shell subsystem in an open-shell embedding: matrix potential
         ! is here a sum of alpha- and beta- counterparts.
         IF ((nspins .EQ. 1) .AND. (opt_dfet_mat%open_shell_embed)) THEN
            CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(1)%matrix, opt_dfet_mat%v_mat(1)%matrix)
            ! Add alpha- and beta- parts
            CALL dbcsr_add(opt_dfet_mat%dfet_mat_term(1)%matrix, opt_dfet_mat%v_mat(2)%matrix,&
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         ENDIF

      ENDDO

      ! Set embedding potentials in the Kphn-Sham environment
      CALL set_ks_env(ks_env=ks_env, dfet_mat_term=opt_dfet_mat%dfet_mat_term)

   END SUBROUTINE dfet_mat_ks_matrix



   FUNCTION add_dfet_mat_to_w_func(opt_dfet_mat, qs_env) RESULT(dfet_mat_corr)
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: trace_tmp, dfet_mat_corr
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_rho_type), POINTER                         :: rho
      CHARACTER(LEN=*), PARAMETER :: routineN = 'add_dfet_mat_to_w_func', routineP = moduleN//':'//routineN
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, rho=rho, dft_control=dft_control, ks_env=ks_env)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      nspins = dft_control%nspins

      dfet_mat_corr = 0.0_dp

      ! If reference system is closed-shell but the embedding is open-shell
      IF ((nspins .EQ. 1) .AND. (opt_dfet_mat%open_shell_embed)) THEN
         CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(1)%matrix, opt_dfet_mat%v_mat(1)%matrix)
         ! Add alpha- and beta- parts
         CALL dbcsr_add(opt_dfet_mat%dfet_mat_term(1)%matrix, opt_dfet_mat%v_mat(2)%matrix,&
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         CALL dbcsr_dot(opt_dfet_mat%dfet_mat_term(1)%matrix, rho_ao(1)%matrix,  trace_tmp)
         dfet_mat_corr = trace_tmp
      ELSE
         DO ispin = 1, nspins
            CALL dbcsr_dot(rho_ao(ispin)%matrix, opt_dfet_mat%dfet_mat_term(ispin)%matrix,  trace_tmp)
            dfet_mat_corr = dfet_mat_corr+trace_tmp
         ENDDO
      ENDIF

      CALL timestop(handle)

   END FUNCTION add_dfet_mat_to_w_func


! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param opt_dmfet ...
! **************************************************************************************************
   SUBROUTINE release_dfet_mat_term(opt_dfet_mat)
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat

      INTEGER                                            :: ispin, nspins

      ! Release dbcsr objects
      nspins = size(opt_dfet_mat%dfet_mat_term)
      DO ispin = 1, nspins
         CALL dbcsr_release(opt_dfet_mat%dfet_mat_term(ispin)%matrix)
         DEALLOCATE (opt_dfet_mat%dfet_mat_term(ispin)%matrix)
      ENDDO
      DEALLOCATE (opt_dfet_mat%dfet_mat_term)

   END SUBROUTINE release_dfet_mat_term

! **************************************************************************************************
!> \brief ...
!> \param opt_dmfet ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE check_dfet_mat(opt_dfet_mat, qs_env)
      TYPE(opt_dfet_mat_pot_type)                           :: opt_dfet_mat
      TYPE(qs_environment_type), POINTER                 :: qs_env

      REAL(KIND=dp)                                      :: diff, diff_beta
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s

      CALL get_qs_env(qs_env, matrix_s=matrix_s)


      CALL dbcsr_dot(opt_dfet_mat%p_diff(1)%matrix, matrix_s(1)%matrix,  diff)
      write(*,*) 'diff', diff

      IF (opt_dfet_mat%open_shell_embed) THEN
         CALL dbcsr_dot(opt_dfet_mat%p_diff(1)%matrix, matrix_s(1)%matrix,  diff_beta)
         write(*,*) 'diff', diff_beta
      ENDIF

   END SUBROUTINE check_dfet_mat

END MODULE optimize_dfet_mat_potential
