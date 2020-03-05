!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE optimize_dfet_mat_potential
   USE pw_methods,                      ONLY: &
        pw_axpy, pw_copy, pw_derive, pw_dr2, pw_integral_ab, pw_integrate_function, pw_scale, &
        pw_transfer, pw_zero
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              copy_dbcsr_to_fm,&
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
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_dmfet_potential'

   PUBLIC :: prepare_dfet_mat_opt, release_dfet_mat_opt,  check_dfet_mat, dfet_mat_opt_step, &
             release_dfet_mat_term, dfet_mat_ks_matrix, add_dfet_mat_to_w_func, print_dfet_mat_opt_info

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

      CALL section_vals_val_get(opt_dfet_mat_section, "DENS_CONV_MAX", &
                                r_val=opt_dfet_mat%conv_max)

      CALL section_vals_val_get(opt_dfet_mat_section, "DENS_CONV_INT", &
                                r_val=opt_dfet_mat%conv_int)

      CALL section_vals_val_get(opt_dfet_mat_section, "SPIN_DENS_CONV_MAX", &
                                r_val=opt_dfet_mat%conv_max_spin)

      CALL section_vals_val_get(opt_dfet_mat_section, "SPIN_DENS_CONV_INT", &
                                r_val=opt_dfet_mat%conv_int_spin)


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

      ! Array to store functional values
      ALLOCATE (opt_dfet_mat%w_func(opt_dfet_mat%n_iter))
      opt_dfet_mat%w_func = 0.0_dp

      ! Allocate max_diff and int_diff
      diff_size = 1
      IF (opt_dfet_mat%open_shell_embed) diff_size = 2
      ALLOCATE (opt_dfet_mat%max_dens_diff(diff_size))
      ALLOCATE (opt_dfet_mat%int_dens_diff(diff_size))

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
      DEALLOCATE (opt_dfet_mat%max_dens_diff)
      DEALLOCATE (opt_dfet_mat%int_dens_diff)

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
      LOGICAL                                            :: geometry_did_change
      INTEGER                                            :: handle, ispin, nspins, step
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao_2d, mat_2d
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gradient

      write(*,*) 'Came to dfet_mat_opt_step?'

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
      geometry_did_change = .FALSE.
      If (step .EQ. 1) THEN
          geometry_did_change = .TRUE.
      Endif
      hf_fraction = 1.0_dp

      ! Compute the gradient in HFX-like update
      CALL hfx_ri_update_ks(qs_env=qs_env, ri_data=opt_dfet_mat%ri_data, ks_matrix=mat_2d, ehfx=ehfx, &
                            rho_ao=rho_ao_2d, geometry_did_change = geometry_did_change, nspins=nspins, &
                            hf_fraction=hf_fraction)
           write(*,*) 'computed overlaps?'
      opt_step = -opt_dfet_mat%trust_rad
      DO ispin = 1, nspins
         CALL dbcsr_add(opt_dfet_mat%v_mat(ispin)%matrix, gradient(ispin)%matrix,alpha_scalar=1.0_dp, &
                        beta_scalar=opt_step)
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
   SUBROUTINE dfet_mat_ks_matrix(qs_env, opt_dfet_mat, i_force_eval, step)
      TYPE(qs_environment_type), POINTER                 :: qs_env !, qs_env_ref
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat
      INTEGER                                            :: i_force_eval, step

      INTEGER                                            :: ispin, nspins
      LOGICAL                                            :: my_change_spin, geometry_did_change, zero_pot
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, term_local
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      REAL(KIND=dp)                                      :: trace, ehfx, hf_fraction
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao_2d, mat_2d

      my_change_spin = .FALSE.
      IF ((opt_dfet_mat%change_spin) .AND. (i_force_eval .EQ. 2)) my_change_spin = .TRUE.

      !
      zero_pot = .FALSE.
      geometry_did_change = .FALSE.
      If (step .EQ. 1) THEN
          zero_pot = .TRUE.
          geometry_did_change = .TRUE.
      Endif

      ! Get Kohn-Sham environment, Kohn-Sham matrix, overlap ...
      NULLIFY (ks_env)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, ks_env=ks_env)
      nspins = 1
      IF (opt_dfet_mat%open_shell_embed) nspins = 2

      NULLIFY (opt_dfet_mat%dfet_mat_term)
      CALL dbcsr_allocate_matrix_set(opt_dfet_mat%dfet_mat_term, nspins)
      DO ispin = 1, nspins
         ALLOCATE (opt_dfet_mat%dfet_mat_term(ispin)%matrix)
         CALL dbcsr_init_p(opt_dfet_mat%dfet_mat_term(ispin)%matrix)
         CALL dbcsr_create(opt_dfet_mat%dfet_mat_term(ispin)%matrix, template=opt_dfet_mat%v_mat(1)%matrix)
         CALL dbcsr_set(opt_dfet_mat%dfet_mat_term(ispin)%matrix, 0.0_dp)
         !ALLOCATE (term_local(ispin)%matrix)
         !CALL dbcsr_init_p(term_local(ispin)%matrix)
         !CALL dbcsr_create(term_local(ispin)%matrix, template=opt_dfet_mat%v_mat(1)%matrix)
         !CALL dbcsr_set(term_local(ispin)%matrix, 0.0_dp)
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

          mat_2d(1:nspins, 1:1) => term_local
          ! Matrix potential serve as rho_ao
          rho_ao_2d(1:nspins, 1:1) => opt_dfet_mat%v_mat
          ! Some settings
          hf_fraction = 1.0_dp

          ! Compute the gradient in HFX-like update
          CALL hfx_ri_update_ks(qs_env=qs_env, ri_data=opt_dfet_mat%ri_data, ks_matrix=mat_2d, ehfx=ehfx, &
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

          ! Release the local term
          DO ispin = 1, nspins
             CALL dbcsr_release_p(term_local(ispin)%matrix)
             CALL dbcsr_release_p(term_local(ispin)%matrix)
          ENDDO
          DEALLOCATE (term_local)

      ENDIF

      ! Set embedding potentials in the Kohn-Sham environment
      CALL set_ks_env(ks_env=ks_env, dfet_mat_term=opt_dfet_mat%dfet_mat_term)


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

      CALL get_qs_env(qs_env, rho=rho, ks_env=ks_env, dft_control=dft_control)

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
!> \param output_unit ...
!> \param step_num ...
! **************************************************************************************************
   SUBROUTINE print_dfet_mat_opt_info(output_unit, step_num, opt_dfet_mat)
      INTEGER                                            :: output_unit, step_num
      TYPE(opt_dfet_mat_pot_type)                        :: opt_dfet_mat

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

      REAL(KIND=dp)                                      :: diff_int, diff_spin_int, diff_max, &
                                                            diff_spin_max,  tot_rho_dummy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_p_type)                                    :: rho_g, rho_r
      INTEGER                                            :: ispin, nspins
      INTEGER                                            :: output_unit
      LOGICAL                                            :: conv_int_diff, conv_max_diff

      !
      nspins = 1
      If (opt_dfet_mat%open_shell_embed) nspins = 2
      !
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
         opt_dfet_mat%int_dens_diff(ispin) = pw_integrate_function(fun=rho_r%pw, oprt='ABS')
         opt_dfet_mat%max_dens_diff(ispin) = max_dens_diff(rho_r)
      ENDDO

      ! Release stuff
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)

      ! Compute target convergence values
      If (opt_dfet_mat%open_shell_embed) THEN ! Open shell
         diff_int = opt_dfet_mat%int_dens_diff(1) + opt_dfet_mat%int_dens_diff(2)
         diff_spin_int = opt_dfet_mat%int_dens_diff(1) - opt_dfet_mat%int_dens_diff(2)
         diff_max = opt_dfet_mat%max_dens_diff(1) + opt_dfet_mat%int_dens_diff(2)
         diff_spin_max = abs(opt_dfet_mat%max_dens_diff(1) - opt_dfet_mat%max_dens_diff(2))
      ELSE ! Closed shell
         diff_int = opt_dfet_mat%int_dens_diff(1)
         diff_max = opt_dfet_mat%max_dens_diff(1)
      ENDIF

      !
      ! Converged?
      !

      ! Maximum value criterium
      ! Open shell
      IF (opt_dfet_mat%open_shell_embed) THEN
         IF ((diff_max .LE. opt_dfet_mat%conv_max) .AND. (diff_spin_max .LE. opt_dfet_mat%conv_max_spin)) THEN
            conv_max_diff = .TRUE.
         ELSE
            conv_max_diff = .FALSE.
         ENDIF
      ELSE
         ! Closed shell
         IF (diff_max .LE. opt_dfet_mat%conv_max) THEN
            conv_max_diff = .TRUE.
         ELSE
            conv_max_diff = .FALSE.
         ENDIF
      ENDIF

      ! Integrated value criterium
      ! Open shell
      IF (opt_dfet_mat%open_shell_embed) THEN
         IF ((diff_int .LE. opt_dfet_mat%conv_int) .AND. (diff_spin_int .LE. opt_dfet_mat%conv_int_spin)) THEN
            conv_int_diff = .TRUE.
         ELSE
            conv_int_diff = .FALSE.
         ENDIF
      ELSE
         ! Closed shell
         IF (diff_int .LE. opt_dfet_mat%conv_int) THEN
            conv_int_diff = .TRUE.
         ELSE
            conv_int_diff = .FALSE.
         ENDIF
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
            " Maximum density difference                = ", diff_max
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Convergence limit for max. density diff.  = ", opt_dfet_mat%conv_max

         IF (opt_dfet_mat%open_shell_embed) THEN

            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Maximum spin density difference           = ", diff_spin_max
            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Convergence limit for max. spin dens.diff.= ", opt_dfet_mat%conv_max_spin

         ENDIF

         IF (conv_max_diff) THEN
            WRITE (UNIT=output_unit, FMT="(T2,2A)") &
               " Convergence in max. density diff.    =     ", &
               "             YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,2A)") &
               " Convergence in max. density diff.    =     ", &
               "              NO"
         END IF

         ! Integrated abs. value of density
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Integrated density difference             = ", diff_int
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Conv. limit for integrated density diff.  = ", opt_dfet_mat%conv_int
         IF (opt_dfet_mat%open_shell_embed) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Integrated spin density difference        = ", diff_spin_int
            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Conv. limit for integrated spin dens.diff.= ", opt_dfet_mat%conv_int_spin
         ENDIF

         IF (conv_int_diff) THEN
            WRITE (UNIT=output_unit, FMT="(T2,2A)") &
               " Convergence in integrated density diff.    =     ", &
               "             YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,2A)") &
               " Convergence in integrated density diff.    =     ", &
               "              NO"
         END IF

      ENDIF

      IF ((opt_dfet_mat%converged) .AND. (output_unit > 0)) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") REPEAT("*", 79)
         WRITE (UNIT=output_unit, FMT="(T2,A,T25,A,T78,A)") &
            "***", "EMBEDDING POTENTIAL OPTIMIZATION COMPLETED", "***"
         WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("*", 79)
      END IF

   END SUBROUTINE check_dfet_mat

END MODULE optimize_dfet_mat_potential
