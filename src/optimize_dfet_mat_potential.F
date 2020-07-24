!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE optimize_dfet_mat_potential
   USE memory_utilities,                ONLY: reallocate
   USE dfet_mat_utils,                  ONLY: mulliken_pot_partitioning
   USE pw_methods,                      ONLY: &
        pw_axpy, pw_copy, pw_derive, pw_dr2, pw_integral_ab, pw_integrate_function, pw_scale, &
        pw_transfer, pw_zero
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              copy_dbcsr_to_fm,&
                                              dbcsr_allocate_matrix_set,&
                                              copy_dbcsr_to_fm
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
                                              dbcsr_release_p, &
                                              dbcsr_get_block_p, &
                                              dbcsr_get_info, &
                                              dbcsr_iterator_blocks_left, &
                                              dbcsr_iterator_next_block, &
                                              dbcsr_iterator_start, &
                                              dbcsr_iterator_stop, &
                                              dbcsr_iterator_type, &
                                              dbcsr_frobenius_norm, &
                                              dbcsr_maxabs, &
                                              dbcsr_scale, &
                                              dbcsr_get_info
   USE embed_types,                     ONLY: opt_dfet_mat_pot_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_get_subs_vals
   USE kinds,                           ONLY: dp,&
                                              default_string_length,&
                                              default_path_length
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE input_constants,                 ONLY: do_potential_id, &
                                              dfet_mat_steep_desc, &
                                              dfet_mat_conj_grad
   USE particle_types,                  ONLY: particle_type
   USE input_cp2k_hfx,                  ONLY: ri_pmat
   USE hfx_types,                       ONLY: hfx_type,&
                                              hfx_ri_type,&
                                              hfx_ri_init_read_input,&
                                              hfx_ri_release
   USE hfx_ri,                          ONLY: hfx_ri_update_ks
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace
   USE qs_collocate_density,            ONLY: calculate_drho_elec,&
                                              calculate_rho_elec
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_release,&
                                              pw_type
   USE optimize_embedding_potential,     ONLY: max_dens_diff
   USE atomic_kind_types,                ONLY: atomic_kind_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, &
        cp_fm_release, cp_fm_set_element, cp_fm_to_fm_submat, cp_fm_type, &
        cp_fm_read_unformatted, cp_fm_write_unformatted
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_files,                        ONLY: close_file,&
                                              open_file
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_dmfet_potential'

   PUBLIC :: prepare_dfet_mat_opt, release_dfet_mat_opt,  check_dfet_mat, dfet_mat_opt_step, &
             release_dfet_mat_term, dfet_mat_ks_matrix, print_dfet_mat_opt_info,&
             print_dfet_mat, given_dfet_mat_pot

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
      INTEGER                                            :: dfet_mat_optimizer, i, n_rep
      INTEGER, DIMENSION(:), POINTER                     :: tmp_v

      ! Read keywords

      CALL section_vals_val_get(opt_dfet_mat_section, "N_ITER", &
                                i_val=opt_dfet_mat%n_iter)

      CALL section_vals_val_get(opt_dfet_mat_section, "TRUST_RAD", &
                                r_val=opt_dfet_mat%trust_rad)

      CALL section_vals_val_get(opt_dfet_mat_section, "DENS_CONV_MAX", &
                                r_val=opt_dfet_mat%conv_max)

      CALL section_vals_val_get(opt_dfet_mat_section, "DENS_CONV_INT", &
                                r_val=opt_dfet_mat%conv_int)

      CALL section_vals_val_get(opt_dfet_mat_section, "SPIN_DENS_CONV_MAX", &
                                r_val=opt_dfet_mat%conv_max_spin)

      CALL section_vals_val_get(opt_dfet_mat_section, "SPIN_DENS_CONV_INT", &
                                r_val=opt_dfet_mat%conv_int_spin)

      CALL section_vals_val_get(opt_dfet_mat_section, "READ_EMBED_POT_MAT", &
                                l_val=opt_dfet_mat%read_embed_pot_mat)

      CALL section_vals_val_get(opt_dfet_mat_section, "OPTIMIZER", i_val=dfet_mat_optimizer)
      SELECT CASE (dfet_mat_optimizer)
      CASE (dfet_mat_steep_desc)
         opt_dfet_mat%sd = .TRUE.
         opt_dfet_mat%cg = .FALSE.
      CASE (dfet_mat_conj_grad)
         opt_dfet_mat%sd = .FALSE.
         opt_dfet_mat%cg = .TRUE.
      CASE DEFAULT
         opt_dfet_mat%sd = .TRUE.
      END SELECT

      CALL section_vals_val_get(opt_dfet_mat_section, "REDUCE_DFET_MAT", &
                                l_val=opt_dfet_mat%reduce_dfet_mat)

      CALL section_vals_val_get(opt_dfet_mat_section, "REDUCE_THRESH", &
                                r_val=opt_dfet_mat%reduce_thresh)

      ALLOCATE (opt_dfet_mat%real_atoms(0))
      NULLIFY (tmp_v)
      CALL section_vals_val_get(opt_dfet_mat_section, "REAL_ATOMS", n_rep_val=n_rep)
      DO i = 1, n_rep
         CALL section_vals_val_get(opt_dfet_mat_section, "REAL_ATOMS", i_vals=tmp_v, i_rep_val=i)
         CALL reallocate(opt_dfet_mat%real_atoms, 1, SIZE(opt_dfet_mat%real_atoms) + SIZE(tmp_v))
         opt_dfet_mat%real_atoms(SIZE(opt_dfet_mat%real_atoms) - SIZE(tmp_v) + 1:SIZE(opt_dfet_mat%real_atoms)) &
         = tmp_v
      END DO
      IF (ASSOCIATED(tmp_v)) THEN
         DEALLOCATE (tmp_v)
      END IF

      NULLIFY(ri_section)
      ri_section => section_vals_get_subs_vals(opt_dfet_mat_section, "RI")

   END SUBROUTINE read_opt_dfet_mat_section


   SUBROUTINE overlap_ri_init(ri_data, ri_section, qs_kind_set, &
                            particle_set, atomic_kind_set, dft_control, para_env, nelectron_total, unit_nr)
      TYPE(hfx_ri_type), POINTER                         :: ri_data
 !     TYPE(hfx_type)                       :: x_data
      TYPE(section_vals_type), POINTER                   :: ri_section
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_para_env_type)                             :: para_env
      INTEGER, INTENT(IN)                                :: nelectron_total

      CHARACTER(LEN=*), PARAMETER :: routineN = 'overlap_ri_init', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: orb_basis_type, ri_basis_type
      INTEGER                                            :: handle, unit_nr, unit_nr_dbcsr
      !TYPE(section_vals_type), POINTER                   :: hf_sub_section
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! Temprorary here
      NULLIFY (ri_data)
      ALLOCATE (ri_data)
      ri_data%hfx_pot%omega = 0.0_dp
      ri_data%hfx_pot%cutoff_radius = 0.0_dp
      ri_data%eps_schwarz = 0.0_dp

      ri_data%ri_section => ri_section
      !ri_data%hfx_section => hfx_section
      !ri_data%eps_schwarz = x_data%screening_parameter%eps_schwarz

      logger => cp_get_default_logger()
      unit_nr_dbcsr = cp_print_key_unit_nr(logger, ri_data%ri_section, "RI_INFO", &
                                           extension=".dbcsrLog")
      ASSOCIATE (hfx_pot=>ri_data%hfx_pot)
         ! Set by force: potential is identity
         hfx_pot%potential_type = do_potential_id
         !hfx_pot%omega = x_data%potential_parameter%omega
         !hfx_pot%cutoff_radius = x_data%potential_parameter%cutoff_radius
      END ASSOCIATE
      !ri_data%eps_schwarz = x_data%screening_parameter%eps_schwarz

      unit_nr_dbcsr = 0

      !unit_nr = -1

      ! Set by force
      orb_basis_type = "ORB"
      ri_basis_type = "RI_HFX"

      CALL hfx_ri_init_read_input(ri_data, ri_section, qs_kind_set, particle_set, atomic_kind_set, &
                                  orb_basis_type, ri_basis_type, dft_control%auto_basis_ri_hfx, &
                                  para_env, unit_nr, unit_nr_dbcsr, nelectron_total)

      CALL timestop(handle)

   END SUBROUTINE overlap_ri_init

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param opt_dmfet ...
!> \param opt_dmfet_section ...
! **************************************************************************************************
   SUBROUTINE prepare_dfet_mat_opt(qs_env, opt_dfet_mat, opt_dfet_mat_section, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      TYPE(section_vals_type), POINTER                   :: opt_dfet_mat_section
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      INTEGER                                            :: nspins, ispin, diff_size
      TYPE(qs_kind_type), DIMENSION(:), &
         POINTER                                         :: qs_kind_set
        TYPE(particle_type), DIMENSION(:), &
         POINTER                                         :: particle_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: nelectron_total
      TYPE(section_vals_type), POINTER                   :: ri_section
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      INTEGER                                            :: unit_nr

      ! Read the input
      CALL read_opt_dfet_mat_section(opt_dfet_mat, opt_dfet_mat_section, ri_section)

      ! Prepare the ri ovelap calculations
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, atomic_kind_set=atomic_kind_set, &
                      dft_control=dft_control, para_env=para_env, nelectron_total=nelectron_total)
      CALL  overlap_ri_init(opt_dfet_mat%ri_data, ri_section, qs_kind_set, &
                            particle_set, atomic_kind_set, dft_control, para_env, nelectron_total, unit_nr)

      nspins = 1
      IF (opt_dfet_mat%open_shell_embed) nspins = 2

      ! Get the overlaps
      CALL get_qs_env(qs_env, matrix_s=matrix_s)

      !NULLIFY (opt_dfet_mat%p_total)
      NULLIFY (opt_dfet_mat%p_diff)
      NULLIFY (opt_dfet_mat%v_mat)
      CALL dbcsr_allocate_matrix_set(opt_dfet_mat%p_diff, nspins)
      CALL dbcsr_allocate_matrix_set(opt_dfet_mat%v_mat, nspins)
      DO ispin = 1, nspins
         ALLOCATE (opt_dfet_mat%p_diff(ispin)%matrix)
         CALL dbcsr_init_p(opt_dfet_mat%p_diff(ispin)%matrix)
         !CALL dbcsr_create(opt_dfet_mat%p_diff(ispin)%matrix, template=matrix_s(1)%matrix)
         !CALL dbcsr_copy(opt_dfet_mat%p_diff(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_create(opt_dfet_mat%p_diff(ispin)%matrix, template=opt_dfet_mat%p_total(1)%matrix)
         CALL dbcsr_copy(opt_dfet_mat%p_diff(ispin)%matrix, opt_dfet_mat%p_total(1)%matrix)
         CALL dbcsr_set(opt_dfet_mat%p_diff(ispin)%matrix, 0.0_dp)
         !
         ALLOCATE (opt_dfet_mat%v_mat(ispin)%matrix)
         CALL dbcsr_init_p(opt_dfet_mat%v_mat(ispin)%matrix)
         CALL dbcsr_create(opt_dfet_mat%v_mat(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(opt_dfet_mat%v_mat(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(opt_dfet_mat%v_mat(ispin)%matrix, 0.0_dp)
         !
      ENDDO
      ! Quantities for conjugate gradient
      If (opt_dfet_mat%cg) then
         NULLIFY (opt_dfet_mat%prev_step)
         NULLIFY (opt_dfet_mat%prev_dir)
         NULLIFY (opt_dfet_mat%dir)
         NULLIFY (opt_dfet_mat%prev_v_mat)
         NULLIFY (opt_dfet_mat%prev_grad)
         CALL dbcsr_allocate_matrix_set(opt_dfet_mat%prev_step, nspins)
         CALL dbcsr_allocate_matrix_set(opt_dfet_mat%prev_dir, nspins)
         CALL dbcsr_allocate_matrix_set(opt_dfet_mat%dir, nspins)
         CALL dbcsr_allocate_matrix_set(opt_dfet_mat%prev_v_mat, nspins)
         CALL dbcsr_allocate_matrix_set(opt_dfet_mat%prev_grad, nspins)
         DO ispin = 1, nspins
            ALLOCATE (opt_dfet_mat%prev_step(ispin)%matrix)
            CALL dbcsr_init_p(opt_dfet_mat%prev_step(ispin)%matrix)
            CALL dbcsr_create(opt_dfet_mat%prev_step(ispin)%matrix, &
                             template=opt_dfet_mat%p_total(1)%matrix)
            CALL dbcsr_set(opt_dfet_mat%prev_step(ispin)%matrix, 0.0_dp)
            ALLOCATE (opt_dfet_mat%prev_dir(ispin)%matrix)
            CALL dbcsr_init_p(opt_dfet_mat%prev_dir(ispin)%matrix)
            CALL dbcsr_create(opt_dfet_mat%prev_dir(ispin)%matrix, &
                             template=opt_dfet_mat%p_total(1)%matrix)
            CALL dbcsr_set(opt_dfet_mat%prev_dir(ispin)%matrix, 0.0_dp)
            ALLOCATE (opt_dfet_mat%dir(ispin)%matrix)
            CALL dbcsr_init_p(opt_dfet_mat%dir(ispin)%matrix)
            CALL dbcsr_create(opt_dfet_mat%dir(ispin)%matrix, &
                             template=opt_dfet_mat%p_total(1)%matrix)
            CALL dbcsr_set(opt_dfet_mat%dir(ispin)%matrix, 0.0_dp)
            ALLOCATE (opt_dfet_mat%prev_v_mat(ispin)%matrix)
            CALL dbcsr_init_p(opt_dfet_mat%prev_v_mat(ispin)%matrix)
            CALL dbcsr_create(opt_dfet_mat%prev_v_mat(ispin)%matrix, &
                             template=opt_dfet_mat%v_mat(1)%matrix)
            CALL dbcsr_set(opt_dfet_mat%prev_v_mat(ispin)%matrix, 0.0_dp)
            ALLOCATE (opt_dfet_mat%prev_grad(ispin)%matrix)
            CALL dbcsr_init_p(opt_dfet_mat%prev_grad(ispin)%matrix)
            CALL dbcsr_create(opt_dfet_mat%prev_grad(ispin)%matrix, &
                             template=opt_dfet_mat%v_mat(1)%matrix)
            CALL dbcsr_set(opt_dfet_mat%prev_grad(ispin)%matrix, 0.0_dp)
         ENDDO
         ! Restrictions on the step size
         opt_dfet_mat%max_step_len = 3.0_dp*opt_dfet_mat%trust_rad
         opt_dfet_mat%min_step_len = 1.0*opt_dfet_mat%trust_rad
      Endif
      ! Array to store functional values
      ALLOCATE (opt_dfet_mat%w_func(opt_dfet_mat%n_iter))
      opt_dfet_mat%w_func = 0.0_dp
      opt_dfet_mat%dfet_mat_corr = 0.0_dp

      ! Allocate max_diff and int_diff
      diff_size = 1
      IF (opt_dfet_mat%open_shell_embed) diff_size = 2
      ALLOCATE (opt_dfet_mat%max_grad(diff_size))
      ALLOCATE (opt_dfet_mat%int_grad(diff_size))

      ! Read potential in matrix form is required
      IF (opt_dfet_mat%read_embed_pot_mat) &
         CALL read_embed_pot_mat(qs_env,opt_dfet_mat_section, opt_dfet_mat%open_shell_embed,&
             opt_dfet_mat%v_mat, .FALSE.)

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
      DEALLOCATE (opt_dfet_mat%max_grad)
      DEALLOCATE (opt_dfet_mat%int_grad)

      ! Release RI data
      CALL hfx_ri_release(opt_dfet_mat%ri_data)
      
      If (opt_dfet_mat%cg) then
         DO ispin = 1, nspins
            CALL dbcsr_release_p(opt_dfet_mat%prev_step(ispin)%matrix)
            CALL dbcsr_release_p(opt_dfet_mat%prev_dir(ispin)%matrix)
            CALL dbcsr_release_p(opt_dfet_mat%dir(ispin)%matrix)
            CALL dbcsr_release_p(opt_dfet_mat%prev_v_mat(ispin)%matrix)
            CALL dbcsr_release_p(opt_dfet_mat%prev_grad(ispin)%matrix)
         ENDDO
         DEALLOCATE (opt_dfet_mat%prev_step)
         DEALLOCATE (opt_dfet_mat%prev_dir)
         DEALLOCATE (opt_dfet_mat%dir)
         DEALLOCATE (opt_dfet_mat%prev_v_mat)
         DEALLOCATE (opt_dfet_mat%prev_grad)
      Endif
 
      DEALLOCATE (opt_dfet_mat%real_atoms)
      DEALLOCATE (opt_dfet_mat%w_func)

   END SUBROUTINE release_dfet_mat_opt

! **************************************************************************************************
!> \brief ...
!> \param opt_dmfet ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE dfet_mat_opt_step(opt_dfet_mat, qs_env, step)
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dfet_opt_step', routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: hf_fraction, ehfx, opt_step
      LOGICAL                                            :: geometry_did_change, calc_grad
      INTEGER                                            :: handle, ispin, nspins, step
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao_2d, mat_2d
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gradient

      CALL timeset(routineN, handle)
      nspins = 1
      IF (opt_dfet_mat%open_shell_embed) nspins = 2

      calc_grad = .TRUE.
      ! If the optimizer is CG and the energy decreases, we trace back and don't recalculate gradient
      If (step .GT. 1) then
         If ((opt_dfet_mat%cg) &
            .AND. ( (opt_dfet_mat%w_func(step) - opt_dfet_mat%w_func(step-1)) .LT. 0.0_dp )) &
             calc_grad = .FALSE.
      Endif

      If (.NOT. calc_grad) then ! CG backtracking
          CALL cg_backtrack(opt_dfet_mat)
      Else ! Everything else: normal CG and steepest descend
         NULLIFY(gradient)
         CALL dbcsr_allocate_matrix_set(gradient, nspins)
         DO ispin = 1, nspins
            ALLOCATE (gradient(ispin)%matrix)
            CALL dbcsr_init_p(gradient(ispin)%matrix)
            CALL dbcsr_create(gradient(ispin)%matrix, template=opt_dfet_mat%v_mat(1)%matrix)
         ENDDO

         ! Make auxiliary matrices
         NULLIFY (mat_2d)
         CALL dbcsr_allocate_matrix_set(mat_2d, nspins, 1)
         NULLIFY (rho_ao_2d)
         CALL dbcsr_allocate_matrix_set(rho_ao_2d, nspins, 1)
         Do ispin = 1, nspins
            ALLOCATE (mat_2d(ispin, 1)%matrix)
            CALL dbcsr_init_p(mat_2d(ispin, 1)%matrix)
            CALL dbcsr_create(mat_2d(ispin, 1)%matrix, template=gradient(ispin)%matrix)
            CALL dbcsr_copy(mat_2d(ispin, 1)%matrix, gradient(ispin)%matrix)
            ! DM differences serve as rho_ao
            ALLOCATE (rho_ao_2d(ispin, 1)%matrix)
            CALL dbcsr_init_p(rho_ao_2d(ispin, 1)%matrix)
            CALL dbcsr_create(rho_ao_2d(ispin, 1)%matrix, template=opt_dfet_mat%p_diff(ispin)%matrix)
            CALL dbcsr_copy(rho_ao_2d(ispin, 1)%matrix, opt_dfet_mat%p_diff(ispin)%matrix)
         Enddo
      
         ! Some settings
         geometry_did_change = .FALSE.
         ! If it's first step and we have restarted, no need to recompute integrals
         If ((step .EQ. 1) .AND. (.NOT. opt_dfet_mat%read_embed_pot_mat)) THEN
             geometry_did_change = .TRUE.
         Endif
         hf_fraction = 1.0_dp
         write(*,*) 'Geometry did change?', geometry_did_change
      
         ! Compute the gradient in HFX-like update
         write(*,*) 'rho_ao_2d', associated(rho_ao_2d)
         write(*,*) dbcsr_frobenius_norm(rho_ao_2d(1,1)%matrix)
         write(*,*) 'mat_2d', associated(mat_2d)
         write(*,*) dbcsr_frobenius_norm(mat_2d(1,1)%matrix)
         write(*,*) 'before update, step'
         ehfx = 0.0_dp
         CALL hfx_ri_update_ks(qs_env=qs_env, ri_data=opt_dfet_mat%ri_data, ks_matrix=mat_2d, ehfx=ehfx, &
                               rho_ao=rho_ao_2d, geometry_did_change = geometry_did_change, nspins=nspins, &
                               hf_fraction=hf_fraction)
         write(*,*) 'after update, step'
         ! Copy from auxiliary structure to gradient
         DO ispin = 1, nspins
            CALL dbcsr_copy(gradient(ispin)%matrix, mat_2d(ispin, 1)%matrix)
         ENDDO
         ! Do conjugate gradient if requested
         If (opt_dfet_mat%cg) then
            CALL cg_step(opt_dfet_mat, gradient)
         Else ! Gradient ascent
            opt_step = -opt_dfet_mat%trust_rad
            DO ispin = 1, nspins
               CALL dbcsr_add(opt_dfet_mat%v_mat(ispin)%matrix, gradient(ispin)%matrix,alpha_scalar=1.0_dp, &
                              beta_scalar=opt_step)
            ENDDO
         Endif
      
         ! Compute norms and deallocate gradient
         DO ispin = 1, nspins
            opt_dfet_mat%int_grad(ispin) = dbcsr_frobenius_norm(gradient(ispin)%matrix)
            opt_dfet_mat%max_grad(ispin) = dbcsr_maxabs(gradient(ispin)%matrix)
            CALL dbcsr_release(gradient(ispin)%matrix)
            DEALLOCATE (gradient(ispin)%matrix)
         ENDDO
         write(*,*) 'step', step
         write(*,*) 'gradient norm', opt_dfet_mat%int_grad
         DEALLOCATE (gradient)
   
         ! Deallocate auxiliary matrices
         DO ispin = 1, nspins
            CALL dbcsr_release(mat_2d(ispin, 1)%matrix)
            DEALLOCATE (mat_2d(ispin, 1)%matrix)
            CALL dbcsr_release(rho_ao_2d(ispin, 1)%matrix)
            DEALLOCATE (rho_ao_2d(ispin, 1)%matrix)
         ENDDO
         DEALLOCATE (mat_2d)
         DEALLOCATE (rho_ao_2d)

      Endif

      CALL timestop(handle)

   END SUBROUTINE dfet_mat_opt_step


   SUBROUTINE cg_step(opt_dfet_mat, gradient)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gradient
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      INTEGER                                            :: ispin, nspins
      REAL(KIND=dp)                                      :: opt_step
      REAL(KIND=dp), DIMENSION(:),ALLOCATABLE            :: beta
      REAL(KIND=dp)                                      :: b, alpha

      nspins = 1
      IF (opt_dfet_mat%open_shell_embed) nspins = 2
         ALLOCATE(beta(nspins))

      If (opt_dfet_mat%i_iter .eq. 1) then ! First step: take steepest ascent step
         opt_step = -opt_dfet_mat%trust_rad
         DO ispin = 1, nspins
            CALL dbcsr_add(opt_dfet_mat%v_mat(ispin)%matrix, gradient(ispin)%matrix,alpha_scalar=1.0_dp, &
                           beta_scalar=opt_step)
            ! Save the step
            CALL dbcsr_add(opt_dfet_mat%prev_step(ispin)%matrix, gradient(ispin)%matrix,alpha_scalar=1.0_dp, &
                              beta_scalar=opt_step)
            ! Save the direction (gradient) and step
            CALL dbcsr_copy(opt_dfet_mat%prev_dir(ispin)%matrix, gradient(ispin)%matrix)
            CALL dbcsr_scale(opt_dfet_mat%prev_dir(ispin)%matrix, -1.0_dp)
            CALL dbcsr_copy(opt_dfet_mat%prev_step(ispin)%matrix, opt_dfet_mat%prev_dir(ispin)%matrix)
            CALL dbcsr_scale(opt_dfet_mat%prev_step(ispin)%matrix, opt_step)
         ENDDO
      Else
         ! Find the CG direction
         CALL Polak_Ribiere(gradient, opt_dfet_mat%prev_dir, opt_dfet_mat%dir, beta, nspins)
         write(*,*) 'beta, cg', beta
         DO ispin = 1, nspins
            b = beta(ispin)
            CALL dbcsr_copy(opt_dfet_mat%dir(ispin)%matrix, gradient(ispin)%matrix)
            CALL dbcsr_add(opt_dfet_mat%dir(ispin)%matrix, opt_dfet_mat%prev_dir(ispin)%matrix,alpha_scalar=-1.0_dp, &
                              beta_scalar=b)
         ENDDO
         ! Update the potential
         CALL line_search(opt_dfet_mat, gradient, nspins, alpha)
         DO ispin = 1, nspins
            ! Save potential
            CALL dbcsr_copy(opt_dfet_mat%prev_v_mat(ispin)%matrix, opt_dfet_mat%v_mat(ispin)%matrix)
            ! Step size is the same for both spins
            CALL dbcsr_add(opt_dfet_mat%v_mat(ispin)%matrix, opt_dfet_mat%dir(ispin)%matrix, &
                          alpha_scalar=1.0_dp, beta_scalar=alpha)
            ! Save quantities for the next iteration
            CALL dbcsr_copy(opt_dfet_mat%prev_dir(ispin)%matrix, opt_dfet_mat%dir(ispin)%matrix)
            CALL dbcsr_set(opt_dfet_mat%prev_step(ispin)%matrix, 0.0_dp)
            CALL dbcsr_add(opt_dfet_mat%prev_step(ispin)%matrix, opt_dfet_mat%dir(ispin)%matrix, &
                          alpha_scalar=1.0_dp, beta_scalar=alpha)
            CALL dbcsr_copy(opt_dfet_mat%prev_grad(ispin)%matrix, gradient(ispin)%matrix)
         ENDDO
 
         DEALLOCATE(beta)
      Endif
 
   END SUBROUTINE cg_step

   SUBROUTINE cg_backtrack(opt_dfet_mat)
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      REAL(kind=dp)                                      :: opt_step
      INTEGER                                            :: ispin, nspins
   
      nspins = 1
      If (opt_dfet_mat%open_shell_embed) nspins = 2

      opt_step = -opt_dfet_mat%trust_rad
      DO ispin = 1, nspins
         CALL dbcsr_copy(opt_dfet_mat%v_mat(ispin)%matrix, opt_dfet_mat%prev_v_mat(ispin)%matrix)
         CALL dbcsr_add(opt_dfet_mat%v_mat(ispin)%matrix, opt_dfet_mat%prev_grad(ispin)%matrix, &
                          alpha_scalar=1.0_dp, beta_scalar=opt_step)
         ! Save direction and step
         CALL dbcsr_copy(opt_dfet_mat%prev_dir(ispin)%matrix, opt_dfet_mat%prev_grad(ispin)%matrix)
         CALL dbcsr_scale(opt_dfet_mat%prev_dir(ispin)%matrix, -1.0_dp)
         CALL dbcsr_copy(opt_dfet_mat%prev_step(ispin)%matrix, opt_dfet_mat%prev_dir(ispin)%matrix)
         CALL dbcsr_scale(opt_dfet_mat%prev_step(ispin)%matrix, opt_step)
      ENDDO

   END SUBROUTINE cg_backtrack


   SUBROUTINE Polak_Ribiere(gradient, prev_dir, dir, beta, nspins)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gradient, prev_dir, dir
      REAL(KIND=dp), DIMENSION(:),ALLOCATABLE            :: beta
      INTEGER                                            :: ispin, nspins
      REAL(kind=dp)                                      :: denom, numer

      DO ispin = 1, nspins
         CALL dbcsr_copy(dir(ispin)%matrix, gradient(ispin)%matrix)
         CALL dbcsr_scale(dir(ispin)%matrix, -1.0_dp)
         CALL dbcsr_add(dir(ispin)%matrix, prev_dir(ispin)%matrix,alpha_scalar=1.0_dp, &
                           beta_scalar=-1.0_dp)
         CALL dbcsr_dot(gradient(ispin)%matrix, dir(ispin)%matrix,  denom)
         CALL dbcsr_dot(prev_dir(ispin)%matrix, prev_dir(ispin)%matrix,  numer)
         beta(ispin) = denom/numer
      ENDDO

   END SUBROUTINE Polak_Ribiere

   SUBROUTINE line_search(opt_dfet_mat, gradient, nspins, alpha)
      REAL(kind=dp)                                      :: a, b, E0, E1, g1, dx, g, x, alpha
      INTEGER                                            :: iter, nspins, ispin
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gradient

      ! Previous and current functional values
      iter = opt_dfet_mat%i_iter
      E0 = opt_dfet_mat%w_func(iter-1) 
      E1 = opt_dfet_mat%w_func(iter)
      ! Norms of the gradient and step 
      g1 = 0.0_dp
      dx = 0.0_dp
      DO ispin = 1, nspins
         g = dbcsr_frobenius_norm(gradient(ispin)%matrix)
         g1 = g1 + g
         x = dbcsr_frobenius_norm(gradient(ispin)%matrix)
         dx = dx + x
      ENDDO
      ! Compute coefficients of quadratic fitting function: y = ax**2+bx+c
      ! c is not needed for finding the extremum
      a = ((E0-E1) + g1*dx)/dx**2
      b = g1 - 2.0_dp*a*dx
      ! Compute the step size alpha
      If (a .NE. 0.0_dp) then
         alpha = - b/(2.0_dp*a)
      Else
         alpha = opt_dfet_mat%min_step_len
      Endif
      If (alpha .GT. opt_dfet_mat%max_step_len) alpha = opt_dfet_mat%max_step_len
      If (alpha .LT. opt_dfet_mat%min_step_len) alpha = opt_dfet_mat%min_step_len
      write(*,*) 'step size', alpha, - b/(2.0_dp*a), 'b=',b, 'a=',a, 'g=',g,'dx=',dx

   END SUBROUTINE line_search

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param opt_dmfet ...
!> \param i_force_eval ...
! **************************************************************************************************
   SUBROUTINE dfet_mat_ks_matrix(qs_env_ref, qs_env, opt_dfet_mat, i_force_eval, step, output_unit)
      TYPE(qs_environment_type), POINTER                 :: qs_env, qs_env_ref
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      INTEGER                                            :: i_force_eval, step

      INTEGER                                            :: ispin, nspins, nspins_subsys, nspins_ref
      LOGICAL                                            :: my_change_spin, geometry_did_change, zero_pot
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: term_local
      TYPE(qs_ks_env_type), POINTER                      :: ks_env, ks_env_ref
      REAL(KIND=dp)                                      :: trace_tmp, ehfx, hf_fraction
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao_2d, mat_2d
      TYPE(dft_control_type), POINTER                    :: dft_control, dft_control_ref
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao_ref
      TYPE(qs_rho_type), POINTER                         :: rho_ref
      INTEGER                                            :: output_unit

      my_change_spin = .FALSE.
      IF ((opt_dfet_mat%change_spin) .AND. (i_force_eval .EQ. 2)) my_change_spin = .TRUE.
      !
      zero_pot = .FALSE.
      geometry_did_change = .FALSE.
      write(*,*) 'read_embed_pot_mat', opt_dfet_mat%read_embed_pot_mat
      If (step .EQ. 1) then
         If (.NOT. opt_dfet_mat%read_embed_pot_mat) then
             zero_pot = .TRUE.
         Else
         ! If we restart, then we need to compute 3-center integrals here
            geometry_did_change = .TRUE.
         Endif
      Endif
      write(*,*) 'zero_pot=', zero_pot
      write(*,*) 'NORM V before hfx', dbcsr_frobenius_norm(opt_dfet_mat%v_mat(1)%matrix)

      ! Global spin states
      nspins = 1
      If (opt_dfet_mat%open_shell_embed) nspins = 2

      ! Subsystem spin states
      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control)
      nspins_subsys = dft_control%nspins
      ! KS embedding term 
      NULLIFY (opt_dfet_mat%dfet_mat_term)
      CALL dbcsr_allocate_matrix_set(opt_dfet_mat%dfet_mat_term, nspins_subsys)
      DO ispin = 1, nspins_subsys
         ALLOCATE (opt_dfet_mat%dfet_mat_term(ispin)%matrix)
         CALL dbcsr_init_p(opt_dfet_mat%dfet_mat_term(ispin)%matrix)
         CALL dbcsr_create(opt_dfet_mat%dfet_mat_term(ispin)%matrix, template=opt_dfet_mat%v_mat(1)%matrix)
         CALL dbcsr_set(opt_dfet_mat%dfet_mat_term(ispin)%matrix, 0.0_dp)
      ENDDO
      !
      ! Apply HFX-like update to compte the dfet_mat_term
      !
      If (.NOT. zero_pot) then

          NULLIFY(term_local)
          CALL dbcsr_allocate_matrix_set(term_local, nspins)
          DO ispin = 1, nspins
             ALLOCATE (term_local(ispin)%matrix)
             CALL dbcsr_init_p(term_local(ispin)%matrix)
             CALL dbcsr_create(term_local(ispin)%matrix, template=opt_dfet_mat%v_mat(1)%matrix)
             CALL dbcsr_set(term_local(ispin)%matrix, 0.0_dp)
          ENDDO
    
          NULLIFY (mat_2d)
          NULLIFY (rho_ao_2d)
          mat_2d(1:nspins, 1:1) => term_local(1:nspins)
          ! Matrix potential serve as rho_ao
          rho_ao_2d(1:nspins, 1:1) => opt_dfet_mat%v_mat(1:nspins)
          ! Some settings
          hf_fraction = 1.0_dp

          ! Compute the term with HFX-like update 
          write(*,*) 'Came to update?, ks_mat'
          CALL hfx_ri_update_ks(qs_env=qs_env_ref, ri_data=opt_dfet_mat%ri_data, ks_matrix=mat_2d, ehfx=ehfx, &
                                rho_ao=rho_ao_2d, geometry_did_change = geometry_did_change, nspins=nspins, &
                                hf_fraction=hf_fraction)
          write(*,*) 'Made the update?, ks_mat'
          DO ispin = 1, nspins

             ! Special case of a closed-shell subsystem in an open-shell embedding: matrix potential
             ! is here a sum of alpha- and beta- counterparts.

             IF ((nspins_subsys .EQ. 1) .AND. (opt_dfet_mat%open_shell_embed)) THEN
                CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(1)%matrix, term_local(1)%matrix)
                ! Add alpha- and beta- parts
                CALL dbcsr_add(opt_dfet_mat%dfet_mat_term(1)%matrix, term_local(2)%matrix,&
                               alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
             ELSE
                IF (.NOT. my_change_spin) THEN
                   CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(ispin)%matrix, term_local(ispin)%matrix)
                ELSE
                   IF (ispin .EQ. 1) THEN
                      CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(ispin)%matrix, term_local(2)%matrix)
                   ENDIF
                   IF (ispin .EQ. 2) THEN
                      CALL dbcsr_copy(opt_dfet_mat%dfet_mat_term(ispin)%matrix, term_local(1)%matrix)
                   ENDIF
                ENDIF ! My_change_spin
             ENDIF ! Special case


          ENDDO

          ! Since the ks matrix term is the same for all systems, we compute its trace with reference
          ! density matrix here, only with one force_eval and print it only once
          If (i_force_eval .EQ. 1) then
            ! Print the embedding KS matrix term
            CALL print_dfet_mat(qs_env, opt_dfet_mat, opt_dfet_mat%i_iter, .False., 'term')
   
            ! Likewise, partitioning should be done for one subsystem:
            ! find atomic contributions to dfet_mat_term
            CALL mulliken_pot_partitioning(qs_env_ref,&
                 opt_dfet_mat%dfet_mat_term, output_unit, 1, opt_dfet_mat%reduce_dfet_mat,&
                 opt_dfet_mat%reduce_thresh, opt_dfet_mat%i_iter, opt_dfet_mat%real_atoms)
 
            CALL get_qs_env(qs_env_ref, rho=rho_ref, ks_env=ks_env_ref, dft_control=dft_control_ref)
            CALL qs_rho_get(rho_ref, rho_ao=rho_ao_ref)
            nspins_ref = dft_control_ref%nspins
            opt_dfet_mat%dfet_mat_corr = 0.0_dp
            trace_tmp = 0.0_dp
            ! If reference system is closed-shell but the embedding is open-shell
            IF ((nspins_ref .EQ. 1) .AND. (opt_dfet_mat%open_shell_embed)) THEN
               ! Alpha- and beta-parts are traced with one density and summed
               CALL dbcsr_dot(opt_dfet_mat%dfet_mat_term(1)%matrix, rho_ao_ref(1)%matrix,  trace_tmp)
               opt_dfet_mat%dfet_mat_corr = opt_dfet_mat%dfet_mat_corr + trace_tmp
               CALL dbcsr_dot(opt_dfet_mat%dfet_mat_term(2)%matrix, rho_ao_ref(1)%matrix,  trace_tmp)
               opt_dfet_mat%dfet_mat_corr = opt_dfet_mat%dfet_mat_corr + trace_tmp
            ELSE ! Closed shell or open shell in open-shell embedding
               DO ispin = 1, nspins
                  CALL dbcsr_dot(rho_ao_ref(ispin)%matrix, opt_dfet_mat%dfet_mat_term(ispin)%matrix,&
                                 trace_tmp)
                  opt_dfet_mat%dfet_mat_corr = opt_dfet_mat%dfet_mat_corr + trace_tmp
               ENDDO
            ENDIF
          Endif ! i_force_eval

          ! Release the local term
          DO ispin = 1, nspins
             CALL dbcsr_release_p(term_local(ispin)%matrix)
          ENDDO
          DEALLOCATE (term_local)

      ENDIF

      ! Set embedding potentials in the Kohn-Sham environment
      CALL set_ks_env(ks_env=ks_env, dfet_mat_term=opt_dfet_mat%dfet_mat_term)

   END SUBROUTINE dfet_mat_ks_matrix

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
! **************************************************************************************************
   SUBROUTINE print_dfet_mat(qs_env, opt_dfet_mat, i_iter, final_one, mat_type)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      INTEGER                                            :: nao
      INTEGER                                            :: i_iter
      LOGICAL                                            :: final_one
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env

      TYPE(cp_fm_type), POINTER                          :: dfet_mat_pot, dfet_mat_pot_beta
      CHARACTER(LEN=default_path_length)                 :: filename 
      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: input
      CHARACTER(LEN=4)                                   :: mat_type
      CHARACTER(LEN=8)                                   :: headline

      NULLIFY (input)
      CALL get_qs_env(qs_env=qs_env, input=input)

      ! We print an unformatted file
         logger => cp_get_default_logger()
         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                              "DFT%QS%OPT_DFET_MAT%EMBED_POT_MAT"), cp_p_file)) THEN
 
            ! Transfer data from dbcsr to fm type to be written
            CALL get_qs_env(qs_env=qs_env, &
                            para_env=para_env)
            NULLIFY (blacs_env)
            CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)
            CALL dbcsr_get_info(opt_dfet_mat%v_mat(1)%matrix, nfullrows_total=nao)

            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=nao, ncol_global=nao)
            NULLIFY (dfet_mat_pot)
            CALL cp_fm_create(dfet_mat_pot, fm_struct, name="dfet_mat_pot")
            IF (mat_type .EQ. 'pot ') CALL copy_dbcsr_to_fm(opt_dfet_mat%v_mat(1)%matrix, dfet_mat_pot)
            IF (mat_type .EQ. 'term') CALL copy_dbcsr_to_fm(opt_dfet_mat%dfet_mat_term(1)%matrix, dfet_mat_pot)
            If (opt_dfet_mat%open_shell_embed) then
               NULLIFY (dfet_mat_pot_beta)
               CALL cp_fm_create(dfet_mat_pot_beta, fm_struct, name="dfet_mat_pot_beta")
               IF (mat_type .EQ. 'pot ') CALL copy_dbcsr_to_fm(opt_dfet_mat%v_mat(2)%matrix, dfet_mat_pot_beta)
               IF (mat_type .EQ. 'term') CALL copy_dbcsr_to_fm(opt_dfet_mat%dfet_mat_term(2)%matrix, dfet_mat_pot_beta)
            Endif
            
            ! Make the names
            IF (mat_type .EQ. 'pot ') headline = 'mat_pot_'
            IF (mat_type .EQ. 'term') headline = 'ks_term_'

            ! Now, write it
            IF (.NOT. final_one) THEN
               WRITE (filename, '(a10,I3.3)') headline, i_iter
            ELSE
               WRITE (filename, '(a10,I3.3)') headline
            ENDIF
            unit_nr = cp_print_key_unit_nr(logger, input, "DFT%QS%OPT_DFET_MAT%EMBED_POT_MAT", extension=".wfn", &
                                     file_form="UNFORMATTED", middle_name=TRIM(filename), file_position="REWIND")
            CALL cp_fm_write_unformatted(dfet_mat_pot, unit_nr)
            If (opt_dfet_mat%open_shell_embed) CALL cp_fm_write_unformatted(dfet_mat_pot_beta, unit_nr)
            IF (unit_nr > 0) THEN
               CALL close_file(unit_nr)
            ENDIF

            ! Release structures
            CALL cp_fm_release(dfet_mat_pot)
            If (opt_dfet_mat%open_shell_embed) CALL cp_fm_release(dfet_mat_pot_beta)
            CALL cp_fm_struct_release(fm_struct)
            CALL cp_blacs_env_release(blacs_env)

         ENDIF ! Should output

   END SUBROUTINE print_dfet_mat

! **************************************************************************************************
!> \brief Read the embedding potential from the binary file
! **************************************************************************************************
   SUBROUTINE read_embed_pot_mat(qs_env, section, open_shell_embed, v_mat, stand_alone)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(LEN=default_path_length)                 :: filename
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_fm_type), POINTER                          :: dfet_mat_pot, dfet_mat_pot_beta
      INTEGER                                            :: restart_unit, nao
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: v_mat
      LOGICAL                                            :: open_shell_embed, stand_alone

      ! Make fm objects
      CALL dbcsr_get_info(v_mat(1)%matrix, nfullrows_total=nao)
      CALL get_qs_env(qs_env=qs_env, &
                      para_env=para_env)
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)
      CALL dbcsr_get_info(v_mat(1)%matrix, nfullrows_total=nao)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                            nrow_global=nao, ncol_global=nao)
      NULLIFY (dfet_mat_pot)
      CALL cp_fm_create(dfet_mat_pot, fm_struct, name="dfet_mat_pot")
      CALL copy_dbcsr_to_fm(v_mat(1)%matrix, dfet_mat_pot)
      If (open_shell_embed) then
         NULLIFY (dfet_mat_pot_beta)
         CALL cp_fm_create(dfet_mat_pot_beta, fm_struct, name="dfet_mat_pot_beta")
         CALL copy_dbcsr_to_fm(v_mat(2)%matrix, dfet_mat_pot_beta)
      Endif
       
      ! Read matrices into fm objects
      IF (para_env%ionode) THEN

         ! Get the restart file name
         CALL dfet_restart_file_name(filename, section, stand_alone)
         restart_unit = -1
         CALL open_file(file_name=filename, &
                        file_action="READ", &
                        file_form="UNFORMATTED", &
                        file_status="OLD", &
                        unit_number=restart_unit)
      ENDIF
          
      CALL cp_fm_read_unformatted(dfet_mat_pot, restart_unit)
      If (open_shell_embed) CALL cp_fm_read_unformatted(dfet_mat_pot_beta, restart_unit)

      ! Close restart file
      IF (para_env%ionode) THEN
         CALL close_file(unit_number=restart_unit)
      ENDIF

      ! Copy from fm to dbcsr
      CALL copy_fm_to_dbcsr(dfet_mat_pot, v_mat(1)%matrix)
      If (open_shell_embed) CALL copy_fm_to_dbcsr(dfet_mat_pot_beta, v_mat(2)%matrix)
   
      ! Release my_embed_pot_coef
      CALL cp_fm_release(dfet_mat_pot)
      IF (open_shell_embed) CALL cp_fm_release(dfet_mat_pot_beta)
      ! Release fm structure
       CALL cp_fm_struct_release(fm_struct)
      ! Release blacs environment
      CALL cp_blacs_env_release(blacs_env)

   END SUBROUTINE read_embed_pot_mat

! **************************************************************************************************
!> \brief Find the embedding restart file name
!> \param filename ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE dfet_restart_file_name(filename, section, stand_alone)
      CHARACTER(LEN=default_path_length), INTENT(OUT)    :: filename
      TYPE(section_vals_type), POINTER                   :: section
      LOGICAL                                            :: stand_alone

      LOGICAL                                            :: exist

      exist = .FALSE.
      If (.NOT. stand_alone) then
         CALL section_vals_val_get(section, "DFET_MAT_RESTART_FILE_NAME", c_val=filename)
      Else
         CALL section_vals_val_get(section, "DFET_MAT_TERM_FILE_NAME", c_val=filename)
      Endif
      INQUIRE (FILE=filename, exist=exist)
      IF (.NOT. exist) &
         CPABORT("DFET matrix embedding restart file not found. ")

   END SUBROUTINE dfet_restart_file_name

! **************************************************************************************************
!> \brief ...
!> \param , output_unitoutput_unit ...
!> \param step_num ...
! **************************************************************************************************
   SUBROUTINE print_dfet_mat_opt_info(output_unit, step_num, opt_dfet_mat)
      INTEGER                                            :: output_unit, step_num
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat

      write(*,*) 'Are we here?', output_unit
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,8('-'),A,I5,1X,12('-'))") &
            "  Optimize matrix embedding potential info at step = ", step_num
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Functional value         = ", opt_dfet_mat%w_func(step_num)
         IF (step_num .GT. 1) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Real energy change         = ", opt_dfet_mat%w_func(step_num) - &
               opt_dfet_mat%w_func(step_num - 1)

           ! WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
           !    " Step size                  = ", opt_dfet_mat%step_len

         END IF

         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Trust radius               = ", opt_dfet_mat%trust_rad

         WRITE (UNIT=output_unit, FMT="(T2,71('-'))")
      END IF

   END SUBROUTINE print_dfet_mat_opt_info

! **************************************************************************************************
!> \brief ...
!> \param opt_dmfet ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE check_dfet_mat(opt_dfet_mat, qs_env, output_unit)
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      TYPE(qs_environment_type), POINTER                 :: qs_env

      REAL(KIND=dp)                                      :: diff_int, diff_max, tot_rho_dummy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_p_type)                                    :: rho_g, rho_r
      INTEGER                                            :: ispin, nspins
      INTEGER                                            :: output_unit
      LOGICAL                                            :: conv_int_diff, conv_max_diff
      REAL(dp), DIMENSION(:), ALLOCATABLE                :: local_max_diff, local_int_diff
      TYPE(dft_control_type), POINTER                    :: dft_control
      !
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = 1
      If (opt_dfet_mat%open_shell_embed) nspins = 2

      If (.NOT. dft_control%qs_control%gapw) then
      ! Only for GPW (no GAPW)
          ALLOCATE(local_max_diff(nspins))
          ALLOCATE(local_int_diff(nspins))
          NULLIFY (ks_env)
          NULLIFY (pw_env)
          CALL get_qs_env(qs_env, ks_env=ks_env, pw_env=pw_env)
          ! pw stuff
          NULLIFY (auxbas_pw_pool)
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
          ! get some of the grids ready
          NULLIFY (rho_r%pw, rho_g%pw)
          CALL pw_pool_create_pw(auxbas_pw_pool, rho_r%pw, &
                                 use_data=REALDATA3D, &
                                 in_space=REALSPACE)
          CALL pw_pool_create_pw(auxbas_pw_pool, rho_g%pw, &
                                 use_data=COMPLEXDATA1D, &
                                 in_space=RECIPROCALSPACE)
          DO ispin = 1, nspins
             CALL calculate_rho_elec(matrix_p=opt_dfet_mat%p_diff(ispin)%matrix, &
                                  rho=rho_r, &
                                  rho_gspace=rho_g, &
                                  total_rho=tot_rho_dummy, &
                                  ks_env=ks_env)
             local_int_diff(ispin) = pw_integrate_function(fun=rho_r%pw, oprt='ABS')
             local_max_diff(ispin) = max_dens_diff(rho_r)
          ENDDO
          ! Release stuff
          CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r%pw)
          CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)
      Endif

      ! Compute target convergence values on gradients
      diff_int = maxval(opt_dfet_mat%int_grad)
      diff_max = maxval(opt_dfet_mat%max_grad)

      !
      ! Converged?
      !

      ! Maximum value criterium
      ! Closed shell
      IF (diff_max .LE. opt_dfet_mat%conv_max) THEN
         conv_max_diff = .TRUE.
      ELSE
         conv_max_diff = .FALSE.
      ENDIF

      ! Frobenius value criterium
      IF (diff_int .LE. opt_dfet_mat%conv_int) THEN
         conv_int_diff = .TRUE.
      ELSE
         conv_int_diff = .FALSE.
      ENDIF

      ! Both creteria satisfied?
      IF ((conv_max_diff) .AND. (conv_int_diff)) THEN
         opt_dfet_mat%converged = .TRUE.
      ELSE
         opt_dfet_mat%converged = .FALSE.
      ENDIF

      ! Print the information
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
            " Convergence check :"

         ! Maximum value of density
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Maximum gradient element            = ", diff_max
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Convergence limit for max. gradient = ", opt_dfet_mat%conv_max

         IF (conv_max_diff) THEN
            WRITE (UNIT=output_unit, FMT="(T2,2A)") &
               " Convergence in max. gradient =     ", &
               "             YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,2A)") &
               " Convergence in max. gradient =     ", &
               "              NO"
         END IF

         ! Frobenius abs. value of density
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Frobenius norm of gradient difference      = ", diff_int
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Conv. limit for Frobenius norm of gradient = ", opt_dfet_mat%conv_int

         IF (conv_int_diff) THEN
            WRITE (UNIT=output_unit, FMT="(T2,2A)") &
               " Convergence in Frobenius norm of gradient =     ", &
               "             YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,2A)") &
               " Convergence in Frobenius norm of gradient =     ", &
               "              NO"
         END IF

         ! Additional information for GPW case
         If (.NOT. dft_control%qs_control%gapw) then
            ! Maximum value of density
            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Maximum density difference                = ", local_max_diff(1)
 
            IF (opt_dfet_mat%open_shell_embed) THEN
 
               WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
                  " Maximum density difference, beta spin      = ", local_max_diff(2)
 
            ENDIF
 
            ! Integrated abs. value of density
            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Integrated density difference             = ", local_int_diff(1)
            IF (opt_dfet_mat%open_shell_embed) THEN
               WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
                  " Integrated density difference, beta spin   = ", local_int_diff(2)
            ENDIF
 
         ENDIF ! GAPW

      ENDIF ! Unit

      If (.NOT. dft_control%qs_control%gapw) then
         DEALLOCATE (local_max_diff)
         DEALLOCATE (local_int_diff)
      ENDIF

      ! Frobenius value criterium
      IF (diff_int .LE. opt_dfet_mat%conv_int) THEN
         conv_int_diff = .TRUE.
      ELSE
         conv_int_diff = .FALSE.
      ENDIF

      ! Both creteria satisfied?
      IF ((conv_max_diff) .AND. (conv_int_diff)) THEN
         opt_dfet_mat%converged = .TRUE.
      ELSE
         opt_dfet_mat%converged = .FALSE.
      ENDIF

      IF ((opt_dfet_mat%converged) .AND. (output_unit > 0)) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") REPEAT("*", 79)
         WRITE (UNIT=output_unit, FMT="(T2,A,T25,A,T78,A)") &
            "***", "EMBEDDING POTENTIAL OPTIMIZATION COMPLETED", "***"
         WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("*", 79)
      END IF

   END SUBROUTINE check_dfet_mat
! **************************************************************************************************
!> \brief Creates and reads matrix potential for a stand-alone embedded calculation
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE given_dfet_mat_pot(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: open_shell_embed
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input, qs_section
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dfet_mat_term, p_total
      INTEGER                                            :: nspins, ispin
      TYPE(qs_rho_type), POINTER                         :: rho

      NULLIFY (input, dft_control,qs_section)
      CALL get_qs_env(qs_env=qs_env, &
                      input=input, &
                      dft_control=dft_control, &
                      ks_env=ks_env)
      qs_section => section_vals_get_subs_vals(input, "DFT%QS")
      open_shell_embed = .FALSE.
      IF (dft_control%nspins .EQ. 2) open_shell_embed = .TRUE.
      nspins = dft_control%nspins

      ! Get the density structure
      NULLIFY (rho, p_total)
      CALL get_qs_env(qs_env, rho=rho)
      NULLIFY (p_total)
      CALL qs_rho_get(rho, rho_ao=p_total)
      ! Create matrices
      NULLIFY (dfet_mat_term)
      CALL dbcsr_allocate_matrix_set(dfet_mat_term, nspins)
      DO ispin = 1, nspins
         ALLOCATE (dfet_mat_term(ispin)%matrix)
         CALL dbcsr_init_p(dfet_mat_term(ispin)%matrix)
         CALL dbcsr_create(dfet_mat_term(ispin)%matrix, template=p_total(1)%matrix)
         CALL dbcsr_set(dfet_mat_term(ispin)%matrix, 0.0_dp)
      ENDDO

      ! Read the matrix
      CALL read_embed_pot_mat(qs_env, qs_section, open_shell_embed, dfet_mat_term, .TRUE.)

      ! Set embedding potentials in the Kohn-Sham environment
      CALL set_ks_env(ks_env=ks_env, dfet_mat_term=dfet_mat_term, given_dfet_mat_term = .TRUE.)

   END  SUBROUTINE given_dfet_mat_pot


END MODULE optimize_dfet_mat_potential
