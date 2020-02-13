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
                                              section_vals_val_get,&
                                              section_vals_get_subs_vals
   USE kinds,                           ONLY: dp,&
                                              default_string_length
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE input_constants,                 ONLY: do_potential_id
   USE particle_types,                  ONLY: particle_type
   USE input_cp2k_hfx,                  ONLY: ri_pmat
   USE hfx_types,                       ONLY: hfx_type,&
                                              hfx_ri_type,&
                                              hfx_ri_init_read_input
   USE hfx_ri,                          ONLY: hfx_ri_update_ks
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
   SUBROUTINE read_opt_dfet_mat_section(opt_dfet_mat, opt_dfet_mat_section, ri_section)
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      TYPE(section_vals_type), POINTER                   :: opt_dfet_mat_section
      TYPE(section_vals_type), POINTER                   :: ri_section

      ! Read keywords

      CALL section_vals_val_get(opt_dfet_mat_section, "N_ITER", &
                                i_val=opt_dfet_mat%n_iter)

      CALL section_vals_val_get(opt_dfet_mat_section, "TRUST_RAD", &
                                r_val=opt_dfet_mat%trust_rad)

      NULLIFY(ri_section)
      ri_section => section_vals_get_subs_vals(opt_dfet_mat_section, "RI")

   END SUBROUTINE read_opt_dfet_mat_section




   SUBROUTINE overlap_ri_init(ri_data, ri_section, qs_kind_set, &
                            particle_set, dft_control, para_env, nelectron_total)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      !TYPE(hfx_type)                       :: x_data
      TYPE(section_vals_type), POINTER                   :: ri_section
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_para_env_type)                             :: para_env
      INTEGER, INTENT(IN)                                :: nelectron_total

      CHARACTER(LEN=*), PARAMETER :: routineN = 'overlap_ri_init', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: orb_basis_type, ri_basis_type
      INTEGER                                            :: handle, unit_nr, unit_nr_dbcsr
      TYPE(section_vals_type), POINTER                   :: hf_sub_section

      CALL timeset(routineN, handle)

      NULLIFY (hf_sub_section)

      ! Temprorary here
      ri_data%hfx_pot%omega = 0.0_dp
      ri_data%hfx_pot%cutoff_radius = 0.0_dp
      ri_data%eps_schwarz = 0.0_dp

      ASSOCIATE (hfx_pot=>ri_data%hfx_pot)
         ! Set by force: potential is identity
         hfx_pot%potential_type = do_potential_id
         !hfx_pot%omega = x_data%potential_parameter%omega
         !hfx_pot%cutoff_radius = x_data%potential_parameter%cutoff_radius
      END ASSOCIATE
      ri_data%ri_section => ri_section
      !ri_data%eps_schwarz = x_data%screening_parameter%eps_schwarz

      unit_nr_dbcsr = 0

      unit_nr = -1

      ! Set by force
      orb_basis_type = "ORB"
      ri_basis_type = "RI_HFX"

      CALL hfx_ri_init_read_input(ri_data, ri_section, qs_kind_set, particle_set, &
                                  orb_basis_type, ri_basis_type, dft_control%auto_basis_ri_hfx, &
                                  para_env, unit_nr, unit_nr_dbcsr, nelectron_total)

      ! Set by force
      !ri_data%flavor = ri_pmat


      CALL timestop(handle)

   END SUBROUTINE overlap_ri_init

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
      TYPE(qs_kind_type), DIMENSION(:), &
         POINTER                                         :: qs_kind_set
        TYPE(particle_type), DIMENSION(:), &
         POINTER                                         :: particle_set
      TYPE(dft_control_type), POINTER          :: dft_control
      TYPE(cp_para_env_type), POINTER                             :: para_env
      INTEGER                                   :: nelectron_total
      TYPE(section_vals_type), POINTER                   :: ri_section
      
      ! Read the input
      CALL read_opt_dfet_mat_section(opt_dfet_mat, opt_dfet_mat_section, ri_section)

      ! Prepare the ri ovelap calculations
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, dft_control=dft_control,&
                      para_env=para_env, nelectron_total=nelectron_total)
      CALL  overlap_ri_init(opt_dfet_mat%ri_data, ri_section, qs_kind_set, &
                            particle_set, dft_control, para_env, nelectron_total)

      nspins = 1
      IF (opt_dfet_mat%open_shell_embed) THEN
         nspins = 2
      ENDIF

      ! Get the orbital coefficients
      CALL get_qs_env(qs_env, matrix_s=matrix_s)

      NULLIFY (opt_dfet_mat%p_total)
      NULLIFY (opt_dfet_mat%p_diff)
      NULLIFY (opt_dfet_mat%v_mat)
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
   SUBROUTINE dfet_mat_opt_step(opt_dfet_mat, qs_env)
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dfet_opt_step', routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: hf_fraction, ehfx, step
      LOGICAL                                            :: geometry_did_change
      INTEGER                                            :: handle, ispin, nspins
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao_2d, mat_2d
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gradient

      CALL timeset(routineN, handle)
      nspins = 1
      IF (opt_dfet_mat%open_shell_embed) nspins = 2

      ! DM difference serve as rho_ao
      !rho_ao => opt_dfet_mat%p_diff
      rho_ao_2d(1:nspins, 1:1) => opt_dfet_mat%p_diff

      NULLIFY(gradient)
      CALL dbcsr_allocate_matrix_set(gradient, nspins)
      DO ispin = 1, nspins
         ALLOCATE (gradient(ispin)%matrix)
         CALL dbcsr_init_p(gradient(ispin)%matrix)
         CALL dbcsr_create(gradient(ispin)%matrix, template=opt_dfet_mat%v_mat(1)%matrix)
      ENDDO
      mat_2d(1:nspins, 1:1) => gradient(1:nspins)

      ! Some settings
      geometry_did_change = .TRUE.
      hf_fraction = 1.0_dp

      ! Compute the gradient in HFX-like update
      CALL hfx_ri_update_ks(qs_env=qs_env, ri_data=opt_dfet_mat%ri_data, ks_matrix=mat_2d, ehfx=ehfx, &
                            rho_ao=rho_ao_2d, geometry_did_change = geometry_did_change, nspins=nspins, &
                            hf_fraction=hf_fraction)

      step = opt_dfet_mat%trust_rad
      DO ispin = 1, nspins
         CALL dbcsr_add(opt_dfet_mat%v_mat(ispin)%matrix, gradient(ispin)%matrix,alpha_scalar=1.0_dp, beta_scalar=step)
      ENDDO

      ! Deallocate gradient
      DO ispin = 1, nspins
         CALL dbcsr_release(gradient(ispin)%matrix)
         DEALLOCATE (gradient(ispin)%matrix)
      ENDDO
      DEALLOCATE (gradient)

      CALL timestop(handle)

   END SUBROUTINE dfet_mat_opt_step

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param opt_dmfet ...
!> \param i_force_eval ...
! **************************************************************************************************
   SUBROUTINE dfet_mat_ks_matrix(qs_env, qs_env_ref, opt_dfet_mat, i_force_eval)
      TYPE(qs_environment_type), POINTER                 :: qs_env, qs_env_ref
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      INTEGER                                            :: i_force_eval

      INTEGER                                            :: ispin, nspins
      LOGICAL                                            :: my_change_spin, geometry_did_change
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, term_local
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      REAL(KIND=dp)                                      :: trace, ehfx, hf_fraction
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao_2d, mat_2d

      my_change_spin = .FALSE.
      IF ((opt_dfet_mat%change_spin) .AND. (i_force_eval .EQ. 2)) my_change_spin = .TRUE.

      ! Get Kohn-Sham environment, Kohn-Sham matrix, overlap ...
      NULLIFY (ks_env)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, ks_env=ks_env)
      nspins = 1
      IF (opt_dfet_mat%open_shell_embed) nspins = 2

      NULLIFY (opt_dfet_mat%dfet_mat_term)
      CALL dbcsr_allocate_matrix_set(opt_dfet_mat%dfet_mat_term, nspins)
      NULLIFY(term_local)
      CALL dbcsr_allocate_matrix_set(term_local, nspins)
      DO ispin = 1, nspins
         ALLOCATE (opt_dfet_mat%dfet_mat_term(ispin)%matrix)
         CALL dbcsr_init_p(opt_dfet_mat%dfet_mat_term(ispin)%matrix)
         !CALL dbcsr_create(opt_dfet_mat%dfet_mat_term(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_create(opt_dfet_mat%dfet_mat_term(ispin)%matrix, template=opt_dfet_mat%v_mat(1)%matrix)
         !CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(ispin)%matrix, matrix_s(1)%matrix)
         !CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(ispin)%matrix, opt_dfet_mat%v_mat(1)%matrix)
         !CALL dbcsr_set(opt_dfet_mat%dfet_mat_term(ispin)%matrix, 0.0_dp)
         ALLOCATE (term_local(ispin)%matrix)
         CALL dbcsr_init_p(term_local(ispin)%matrix)
         !CALL dbcsr_create(term_local(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_create(term_local(ispin)%matrix, template=opt_dfet_mat%v_mat(1)%matrix)
         !CALL dbcsr_copy(term_local(ispin)%matrix, matrix_s(1)%matrix)
         !CALL dbcsr_copy(term_local(ispin)%matrix, opt_dfet_mat%v_mat(1)%matrix)
         !CALL dbcsr_set(term_local(ispin)%matrix, 0.0_dp)
      ENDDO

      !
      ! Apply HFX-like update to compte the dfet_mat_term
      !
      mat_2d(1:nspins, 1:1) => term_local
      ! Matrix potential serve as rho_ao
      rho_ao_2d(1:nspins, 1:1) => opt_dfet_mat%v_mat
      ! Some settings
      geometry_did_change = .TRUE.
      hf_fraction = 1.0_dp

      ! Compute the gradient in HFX-like update
      CALL hfx_ri_update_ks(qs_env=qs_env_ref, ri_data=opt_dfet_mat%ri_data, ks_matrix=mat_2d, ehfx=ehfx, &
                            rho_ao=rho_ao_2d, geometry_did_change = geometry_did_change, nspins=nspins, &
                            hf_fraction=hf_fraction)

      DO ispin = 1, nspins

         IF (.NOT. my_change_spin) THEN
            CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(ispin)%matrix, term_local(ispin)%matrix)
         ELSE
            IF (ispin .EQ. 1) THEN
               CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(ispin)%matrix, term_local(2)%matrix)
            ENDIF
            IF (ispin .EQ. 2) THEN
               CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(ispin)%matrix, term_local(1)%matrix)
            ENDIF
         ENDIF

         ! Special case of a closed-shell subsystem in an open-shell embedding: matrix potential
         ! is here a sum of alpha- and beta- counterparts.
         IF ((nspins .EQ. 1) .AND. (opt_dfet_mat%open_shell_embed)) THEN
            CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(1)%matrix, term_local(1)%matrix)
            ! Add alpha- and beta- parts
            CALL dbcsr_add(opt_dfet_mat%dfet_mat_term(1)%matrix, term_local(2)%matrix,&
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         ENDIF

      ENDDO
      CALL dbcsr_trace(opt_dfet_mat%dfet_mat_term(1)%matrix, trace)
      write(*,*) 'dfet_mat_Term', trace

      ! Set embedding potentials in the Kohn-Sham environment
      CALL set_ks_env(ks_env=ks_env, dfet_mat_term=opt_dfet_mat%dfet_mat_term)

      ! Release the local term
      DO ispin = 1, nspins
         CALL dbcsr_release_p(term_local(ispin)%matrix)
         CALL dbcsr_release_p(term_local(ispin)%matrix)
      ENDDO
      DEALLOCATE (term_local)

   END SUBROUTINE dfet_mat_ks_matrix



   FUNCTION add_dfet_mat_to_w_func(opt_dfet_mat, qs_env) RESULT(dfet_mat_corr)
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: handle, ispin, nspins
      TYPE(dft_control_type), POINTER                    :: dft_control
      REAL(KIND=dp)                                      :: trace_tmp, dfet_mat_corr
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(qs_rho_type), POINTER                         :: rho
      CHARACTER(LEN=*), PARAMETER :: routineN = 'add_dfet_mat_to_w_func', routineP = moduleN//':'//routineN
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, rho=rho, ks_env=ks_env)

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

      REAL(KIND=dp)                                      :: diff, diff_beta, nel, nel_tot, trace, s_sq, dp_sq
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s

      CALL get_qs_env(qs_env, matrix_s=matrix_s)


      CALL dbcsr_dot(opt_dfet_mat%p_diff(1)%matrix, matrix_s(1)%matrix,  diff)
      CALL dbcsr_dot(opt_dfet_mat%p_diff(1)%matrix, opt_dfet_mat%p_diff(1)%matrix,  dp_sq)
      CALL dbcsr_dot(opt_dfet_mat%p_subsys(1)%matrix, matrix_s(1)%matrix,  nel)
      CALL dbcsr_dot(opt_dfet_mat%p_total(1)%matrix, matrix_s(1)%matrix,  nel_tot)
      CALL dbcsr_dot(matrix_s(1)%matrix, matrix_s(1)%matrix,  s_sq)
      CALL dbcsr_trace(opt_dfet_mat%p_diff(1)%matrix, trace)
      write(*,*) 'diff', diff, 'n el', nel, nel_tot, 'diff trace', trace, 's_s_trace', s_sq, 'deltap*deltap trace',dp_sq
      CALL dbcsr_trace(matrix_s(1)%matrix, trace)
      write(*,*) 'overlap trace', trace

      IF (opt_dfet_mat%open_shell_embed) THEN
         CALL dbcsr_dot(opt_dfet_mat%p_diff(1)%matrix, matrix_s(1)%matrix,  diff_beta)
         write(*,*) 'diff', diff_beta
      ENDIF

   END SUBROUTINE check_dfet_mat

END MODULE optimize_dfet_mat_potential
