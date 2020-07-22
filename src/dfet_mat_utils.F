MODULE dfet_mat_utils

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix,&
                                              write_fm_with_basis_info
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type,&
                                              cp_fm_set_submatrix
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_deallocate_matrix, dbcsr_get_block_p, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_set, dbcsr_setname, dbcsr_type, dbcsr_get_info
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush
   USE message_passing,                 ONLY: mp_sum
   USE orbital_pointers,                ONLY: nso
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE scf_control_types,               ONLY: scf_control_type
   USE particle_list_types,             ONLY: particle_list_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'population_analyses'

   PUBLIC :: mulliken_pot_partitioning

   CONTAINS

! **************************************************************************************************
!> \brief  Hacked clone of Mulliken population analysis for embedding potetnial
!>
!> \param qs_env ...
!> \param output_unit ...
!> \param print_level ...
!> \date    07.2020
!> \author  Vladimir Rybkin (hacked)
!> \author  Matthias Krack (MK)
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE mulliken_pot_partitioning(qs_env, v_mat, output_unit, print_level, do_reduce, thresh, i_iter, real_atoms)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: output_unit, print_level

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mulliken_pot_partitioning', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, handle, iatom, isgf, ispin, &
                                                            jatom, jsgf, natom, nsgf, nspin, sgfa, &
                                                            sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf_atom
      LOGICAL                                            :: print_gop
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pot_cont, v_block
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: v_mat
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: atom_cont
      REAL(KIND=dp)                                      :: thresh
      INTEGER, DIMENSION(:)                              :: real_atoms
      INTEGER                                            :: i_iter
      LOGICAL                                            :: do_reduce

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set)
      NULLIFY (qs_kind_set)
      NULLIFY (pot_cont)
      NULLIFY (particle_set)
      NULLIFY (v_block)
      NULLIFY (para_env)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      para_env=para_env)

      CPASSERT(ASSOCIATED(atomic_kind_set))
      CPASSERT(ASSOCIATED(qs_kind_set))
      CPASSERT(ASSOCIATED(particle_set))

      nspin = SIZE(v_mat, 1)

      ! Get the total number of contracted spherical Gaussian basis functions
      CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)
      ALLOCATE (first_sgf_atom(natom))
      first_sgf_atom(:) = 0

      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf_atom)

      ! Provide an array to store the potential contributions for each spin
      ALLOCATE (pot_cont(nsgf, nspin))
      pot_cont(:, :) = 0.0_dp

      ! Write headline
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
            '!-----------------------------------------------------------------------------!'
         WRITE (UNIT=output_unit, FMT="(T22,A)") "Mulliken Atomic Contributions to V(embed)"
      END IF

      ! Build Mulliken population matrix for each spin
      DO ispin = 1, nspin
         CALL dbcsr_iterator_start(iter, v_mat(1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, iatom, jatom, v_block, blk)
            IF (.NOT. (ASSOCIATED(v_block))) CYCLE

            sgfb = first_sgf_atom(jatom)
            DO jsgf = 1, SIZE(v_block, 2)
               DO isgf = 1, SIZE(v_block, 1)
                  pot_cont(sgfb, ispin) = pot_cont(sgfb, ispin) + abs(v_block(isgf, jsgf))
               END DO ! isgf
               sgfb = sgfb + 1
            END DO ! jsgf
            IF (iatom /= jatom) THEN
               sgfa = first_sgf_atom(iatom)
               DO isgf = 1, SIZE(v_block, 1)
                  DO jsgf = 1, SIZE(v_block, 2)
                     pot_cont(sgfa, ispin) = pot_cont(sgfa, ispin) + abs(v_block(isgf, jsgf))
                  END DO ! jsgf
                  sgfa = sgfa + 1
               END DO ! isgf
            END IF ! dbcsr_iterator
         END DO
         CALL dbcsr_iterator_stop(iter)

      END DO ! ispin

      CALL mp_sum(pot_cont, para_env%group)

      ! Allocate atomic contribution
      ALLOCATE(atom_cont(natom))
      ! Write atomic populations and charges
      IF (output_unit > 0) THEN
         print_gop = (print_level > 1) ! Print also orbital populations
         CALL write_pot_cont(pot_cont, atomic_kind_set, qs_kind_set, particle_set, &
                             output_unit, print_gop, atom_cont, para_env)
      END IF

      ! Reduce potential
      If (do_reduce) CALL reduce_potential(qs_env, v_mat, atom_cont, natom, thresh,&
                          real_atoms, output_unit, i_iter)
      write(*,*) 'atom_cont_outside', atom_cont

      ! Release stuff
      IF (ASSOCIATED(pot_cont)) THEN
         DEALLOCATE (pot_cont)
      END IF
      IF (ALLOCATED(first_sgf_atom)) THEN
         DEALLOCATE (first_sgf_atom)
      END IF
      !IF (ASSOCIATED(atom_cont)) THEN
         DEALLOCATE (atom_cont)
      !END IF


      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(T2,A)") &
            '!-----------------------------------------------------------------------------!'
      END IF

      CALL timestop(handle)

   END SUBROUTINE mulliken_pot_partitioning

! **************************************************************************************************
!> \brief Write atomic orbital contributions to embedding potential
!>
!> \param orbpop ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param output_unit ...
!> \param print_orbital_contributions ...
!> \date    07.07.2010
!> \author  Matthias Krack (MK)
!> \author  Adopted here by Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE write_pot_cont(pot_cont, atomic_kind_set, qs_kind_set, particle_set, output_unit, &
                           print_orbital_contributions, atom_cont, para_env)

      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pot_cont
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: output_unit
      LOGICAL, INTENT(IN)                                :: print_orbital_contributions

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_pot_cont', routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=6), DIMENSION(:), POINTER            :: sgf_symbol
      INTEGER                                            :: handle, iao, iatom, ikind, iset, isgf, &
                                                            ishell, iso, l, natom, nset, nsgf, &
                                                            nspin
      INTEGER, DIMENSION(:), POINTER                     :: nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: lshell
      REAL(KIND=dp), DIMENSION(2)                        :: sum_pot_cont, totsum_pot_cont
      REAL(KIND=dp), DIMENSION(:)                        :: atom_cont
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (lshell)
      NULLIFY (nshell)
      NULLIFY (orb_basis_set)
      NULLIFY (sgf_symbol)

      CPASSERT(ASSOCIATED(pot_cont))
      CPASSERT(ASSOCIATED(atomic_kind_set))
      CPASSERT(ASSOCIATED(particle_set))

      nspin = SIZE(pot_cont, 2)

      CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

      ! Select and write headline
      IF (print_orbital_contributions) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
            "# Orbital  AO symbol  Orbital contribution                            Net value"
      ELSE
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
            "#  Atom  Element  Kind  Atomic contribution                           Net value"
      END IF

      totsum_pot_cont(:) = 0.0_dp
      atom_cont = 0.0_dp

      iao = 1
      DO iatom = 1, natom
         sum_pot_cont(:) = 0.0_dp
         NULLIFY (orb_basis_set)
         CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                              element_symbol=element_symbol, &
                              kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                   nset=nset, &
                                   nshell=nshell, &
                                   l=lshell, &
                                   sgf_symbol=sgf_symbol)
            isgf = 1
            DO iset = 1, nset
               DO ishell = 1, nshell(iset)
                  l = lshell(ishell, iset)
                  DO iso = 1, nso(l)
                     IF (nspin == 1) THEN
                        sum_pot_cont(1) = sum_pot_cont(1) + pot_cont(iao, 1)
                        IF (print_orbital_contributions) THEN
                           IF (isgf == 1) WRITE (UNIT=output_unit, FMT="(A)") ""
                           WRITE (UNIT=output_unit, &
                                  FMT="(T2,I9,2X,A2,1X,A,T30,F12.6)") &
                              iao, element_symbol, sgf_symbol(isgf), pot_cont(iao, 1)
                        END IF
                     ELSE
                        sum_pot_cont(1:2) = sum_pot_cont(1:2) + pot_cont(iao, 1:2)
                        IF (print_orbital_contributions) THEN
                           IF (isgf == 1) WRITE (UNIT=output_unit, FMT="(A)") ""
                           WRITE (UNIT=output_unit, &
                                  FMT="(T2,I9,2X,A2,1X,A,T30,F12.6)") &
                              iao, element_symbol, sgf_symbol(isgf), &
                              pot_cont(iao, 1)+pot_cont(iao,2)
                        END IF
                     END IF
                     isgf = isgf + 1
                     iao = iao + 1
                  END DO
               END DO
            END DO
            IF (nspin == 1) THEN
               totsum_pot_cont(1) = totsum_pot_cont(1) + sum_pot_cont(1)
               ! Save atom contribution
               atom_cont(iatom) = sum_pot_cont(1)
               WRITE (UNIT=output_unit, &
                      FMT="(T2,I7,5X,A2,2X,I6,T30,F12.6,T68)") &
                  iatom, element_symbol, ikind, sum_pot_cont(1)
            ELSE
               totsum_pot_cont(1:2) = totsum_pot_cont(1:2) + sum_pot_cont(1:2)
               ! Save atom contribution
               atom_cont(iatom) = sum_pot_cont(1) + sum_pot_cont(2)
               WRITE (UNIT=output_unit, &
                      FMT="(T2,I7,5X,A2,2X,I6,T30,F12.6,T68)") &
                  iatom, element_symbol, ikind, sum_pot_cont(1)+sum_pot_cont(2)
            END IF
         END IF ! atom has an orbital basis
      END DO ! next atom iatom

      CALL mp_sum(atom_cont, para_env%group)
      write(*,*) 'atom_cont_inside', atom_cont


      ! Write total sums
      IF (print_orbital_contributions) WRITE (UNIT=output_unit, FMT="(A)") ""
      IF (nspin == 1) THEN
         WRITE (UNIT=output_unit, &
                FMT="(T2,A,T42,F12.6,T68)") &
            "# Sum of atomic contributions", totsum_pot_cont(1)
      ELSE
         WRITE (UNIT=output_unit, &
                FMT="(T2,A,T28,4(1X,F12.6),/)") &
            "# Sum of atomic contributions", totsum_pot_cont(1) + totsum_pot_cont(2)
      END IF

      IF (output_unit > 0) CALL m_flush(output_unit)

      CALL timestop(handle)

   END SUBROUTINE write_pot_cont


   SUBROUTINE reduce_potential(qs_env, v_mat, atom_cont, natom, thresh, real_atoms, output_unit, i_iter)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mulliken_pot_partitioning', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, natom, i
      INTEGER                                            :: nao, i_iter
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: v_mat
      REAL(KIND=dp)                                      :: thresh, pot_total
      INTEGER                                            :: red_num_atom, pos
      REAL(KIND=dp), DIMENSION(:)                        :: atom_cont
      INTEGER(KIND=dp), DIMENSION(:), ALLOCATABLE        :: red_atom_list
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE         :: unit_block
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: atom_mask, real_atom_mask
      INTEGER                                            :: reduced_dim
      INTEGER                                            :: unit_size, start_row, start_col
      TYPE(cp_fm_type), POINTER                          :: Pt, VPt, PVPt
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      INTEGER, ALLOCATABLE, DIMENSION(:,:)               :: full_ind, red_ind
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      INTEGER, DIMENSION(:)                              :: real_atoms

      write(*,*) 'i_iter', i_iter

      CALL timeset(routineN, handle)
      ! Find the new, reduced atom list with significant contributions:
      ! real atoms of the subsystem must be forcefully included
      ALLOCATE(atom_mask(natom))
      ALLOCATE(real_atom_mask(natom))
      real_atom_mask = 0
      DO iatom = 1, size(real_atoms)
         real_atom_mask(real_atoms(iatom)) = 1
      ENDDO
      atom_mask = 0
      red_num_atom = 0
      pot_total = sum(atom_cont)
      write(*,*) 'thresh', thresh, 'pot_total', pot_total, 'atom_cont', atom_cont
      DO iatom = 1, natom
         If (((atom_cont(iatom)/pot_total) .GE. thresh) .OR. (real_atom_mask(iatom) .EQ. 1)) then
            red_num_atom = red_num_atom + 1
            atom_mask(iatom) = 1
         Endif
      ENDDO
      ALLOCATE(red_atom_list(red_num_atom))
      ! Now, fill the shorter atom list
      pos = 0
      DO iatom = 1, natom
         If (atom_mask(iatom) .EQ. 1) THEN
             pos = pos + 1
             red_atom_list(pos) = iatom
         Endif
      ENDDO
      ! Allocate indices: full and reduced
      ALLOCATE(full_ind(red_num_atom, 2))
      ALLOCATE(red_ind(red_num_atom, 2))

      CALL find_reduced_basis_size(qs_env, atom_mask, reduced_dim, full_ind, red_ind)
      write(*,*) 'reduced dimension', reduced_dim
      write(*,*) 'red_atom_list', red_atom_list

     ! Make the projectors
     CALL get_qs_env(qs_env=qs_env, &
                      para_env=para_env)
     NULLIFY (blacs_env)
     CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)
     CALL dbcsr_get_info(v_mat(1)%matrix, nfullrows_total=nao)

     ! Projector
     NULLIFY (fm_struct)
     CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=nao, ncol_global=reduced_dim)
     NULLIFY (Pt)
     CALL cp_fm_create(Pt, fm_struct, name="projector")
     CALL cp_fm_set_all(Pt, 0.0_dp)
     ! P*V
     NULLIFY (VPt)
     CALL cp_fm_create(VPt, fm_struct, name="pot_projector")
     CALL cp_fm_set_all(VPt, 0.0_dp)
     CALL cp_fm_struct_release(fm_struct)
     ! V_reduced = P*V*P^t: reduced potential in full matrix format
     NULLIFY (fm_struct)
     CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=reduced_dim, ncol_global=reduced_dim)
     NULLIFY (PVPt)
     CALL cp_fm_create(PVPt, fm_struct, name="v_red")
     CALL cp_fm_struct_release(fm_struct)
     CALL cp_fm_set_all(PVPt, 0.0_dp)

     ! Fill in the projector (transpose)
     DO iatom = 1, red_num_atom
        unit_size = red_ind(iatom,2) - red_ind(iatom,1) + 1
        ALLOCATE(unit_block(unit_size, unit_size))
        unit_block = 0.0_dp
        DO i = 1, unit_size
           unit_block(i,i) = 1.0_dp
        ENDDO
        start_row = full_ind(iatom, 1)
        start_col = red_ind(iatom, 1)
        CALL cp_fm_set_submatrix(Pt, unit_block, start_row, start_col, unit_size, unit_size)
     END DO
     write(*,*) 'projector', Pt%local_data(:,:)
     ! Multiply, first sparse by dense: V*P^t
     CALL cp_dbcsr_sm_fm_multiply(v_mat(1)%matrix, Pt, VPt, ncol=nao)
     ! Dense by dense: P*V*P^t
     CALL cp_gemm('T', 'N', reduced_dim, reduced_dim, nao, 1.0_dp, Pt, VPt, 0.0_dp, PVPt)

     ! Print coordinates for stand-alone calculation
     CALL print_coord_for_standalone(qs_env, atom_mask, real_atom_mask, output_unit)

      ! Deallocate
      DEALLOCATE(atom_mask)
      DEALLOCATE(real_atom_mask)
      DEALLOCATE(full_ind)
      DEALLOCATE(red_ind)

      CALL cp_fm_release(Pt)
      CALL cp_fm_release(VPt)
      CALL cp_fm_release(PVPt)

      CALL cp_blacs_env_release(blacs_env)

      IF (ALLOCATED(red_atom_list)) THEN
         DEALLOCATE (red_atom_list)
      END IF
      CALL timestop(handle)

   END SUBROUTINE reduce_potential




   SUBROUTINE find_reduced_basis_size(qs_env, atom_mask, reduced_dim, full_ind, red_ind)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: reduced_dim

      INTEGER                                            :: iatom, ikind, natom, nsgf, i_red_at, &
                                                            full_dim, start_red, start_full
      INTEGER, ALLOCATABLE, DIMENSION(:)               :: kind_of
      INTEGER, DIMENSION(:)                              :: atom_mask
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, ALLOCATABLE, DIMENSION(:,:)                 :: full_ind, red_ind

      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set)

      natom = SIZE(particle_set)
      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      reduced_dim = 0
      write(*,*) 'natom = ', natom
      start_red = 0
      start_full = 0
      i_red_at = 0
      full_dim = 0
 
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), nsgf=nsgf, basis_type="ORB")
         full_dim = full_dim + nsgf
         If (atom_mask(iatom) .EQ. 1) then
            reduced_dim = reduced_dim + nsgf
            i_red_at = i_red_at + 1
            full_ind(i_red_at,1) = start_full + 1
            full_ind(i_red_at,2) = full_dim
            red_ind(i_red_at,1) = start_red + 1
            red_ind(i_red_at,2) = reduced_dim
            start_red = reduced_dim
         Endif
         start_full = full_dim
      END DO

      DEALLOCATE (kind_of)

   END SUBROUTINE find_reduced_basis_size

   SUBROUTINE print_coord_for_standalone(qs_env, atom_mask, real_atom_mask, output_unit)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(qs_subsys_type), POINTER                      :: subsys
      INTEGER, DIMENSION(:)                              :: atom_mask, real_atom_mask
      CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:)        :: particles_el
      INTEGER                                            :: iat, n, output_unit

      CALL get_qs_env(qs_env=qs_env, subsys=subsys)
      CALL qs_subsys_get(subsys=subsys, particles=particles)

      ! Prepare the file
      IF (output_unit > 0) THEN
          WRITE (UNIT=output_unit, FMT="(T22,A)") "Coordinates for stand-alone embedded calculation"

         n = particles%n_els
         ALLOCATE (particles_el(n))
         DO iat = 1, n
            CALL get_atomic_kind(particles%els(iat)%atomic_kind, element_symbol=particles_el(iat))
         END DO

         ! Print coordinates only of the included atoms to file
         DO iat = 1, SIZE(particles_el)
            If (atom_mask(iat) .EQ. 1) then ! Is atom included?
               If (real_atom_mask(iat) .EQ. 1) then ! Is it real?
                  WRITE (output_unit, '(a4,4f12.6)') particles_el(iat), particles%els(iat)%r(:)
               Else ! Ghost
                  WRITE (output_unit, '(a2,a8,4f12.6)') particles_el(iat),'_ghost', particles%els(iat)%r(:)
               Endif
            ENDIF
         ENDDO

         DEALLOCATE (particles_el)

      ENDIF ! Should output

   END SUBROUTINE print_coord_for_standalone

END MODULE dfet_mat_utils
