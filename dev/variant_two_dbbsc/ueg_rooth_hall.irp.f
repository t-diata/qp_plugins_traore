subroutine ueg_Roothaan_Hall_SCF

BEGIN_DOC
! copy from scf_utils/roothaan_hall_scf.irp.f
! Roothaan-Hall algorithm for SCF Hartree-Fock calculation
END_DOC

  implicit none

  double precision               :: energy_SCF,energy_SCF_previous,Delta_energy_SCF
  double precision               :: max_error_DIIS,max_error_DIIS_alpha,max_error_DIIS_beta
  double precision, allocatable  :: Fock_matrix_DIIS(:,:,:),error_matrix_DIIS(:,:,:)

  integer                        :: iteration_SCF,dim_DIIS,index_dim_DIIS

  logical                        :: converged
  integer                        :: i,j,m
  logical, external              :: qp_stop
  double precision, allocatable :: mo_coef_save(:,:), S(:,:)

  PROVIDE ao_md5 mo_occ level_shift

  allocate(mo_coef_save(ao_num,mo_num),                          &
      Fock_matrix_DIIS (ao_num,ao_num,max_dim_DIIS),                 &
      error_matrix_DIIS(ao_num,ao_num,max_dim_DIIS)                  &
      )

  Fock_matrix_DIIS = 0.d0
  error_matrix_DIIS = 0.d0
  mo_coef_save = 0.d0

  call write_time(6)

  print*,'Energy of the guess = ', ueg_SCF_energy
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '  N ', 'Energy  ', 'Energy diff  ',  'DIIS error  ', 'Level shift   '
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'

! Initialize energies and density matrices
  energy_SCF_previous = ueg_SCF_energy
  Delta_energy_SCF    = 1.d0
  iteration_SCF       = 0
  dim_DIIS            = 0
  max_error_DIIS      = 1.d0


!
! Start of main SCF loop
!
  PROVIDE FPS_SPF_matrix_AO ueg_Fock_matrix_AO 

  converged = .False.
  do while ( .not.converged .and. (iteration_SCF < n_it_SCF_max) )

! Increment cycle number

    iteration_SCF += 1
    if(frozen_orb_scf)then
     call initialize_mo_coef_begin_iteration
    endif

! Current size of the DIIS space

    dim_DIIS = min(dim_DIIS+1,max_dim_DIIS)

    if ( (scf_algorithm == 'DIIS').and.(dabs(Delta_energy_SCF) > 1.d-6) )  then

      ! Store Fock and error matrices at each iteration
      do j=1,ao_num
        do i=1,ao_num
          index_dim_DIIS = mod(dim_DIIS-1,max_dim_DIIS)+1
          Fock_matrix_DIIS (i,j,index_dim_DIIS) = ueg_Fock_matrix_AO(i,j)
          error_matrix_DIIS(i,j,index_dim_DIIS) = FPS_SPF_matrix_AO(i,j)
        enddo
      enddo

      ! Compute the extrapolated Fock matrix

      call extrapolate_Fock_matrix(                                    &
          error_matrix_DIIS,Fock_matrix_DIIS,                          &
          ueg_Fock_matrix_AO,size(ueg_Fock_matrix_AO,1),                       &
          iteration_SCF,dim_DIIS                                       &
          )

      ueg_Fock_matrix_AO_alpha = ueg_Fock_matrix_AO*0.5d0
      ueg_Fock_matrix_AO_beta  = ueg_Fock_matrix_AO*0.5d0
      TOUCH ueg_Fock_matrix_AO_alpha ueg_Fock_matrix_AO_beta

    endif

    MO_coef = eigenvectors_Fock_matrix_MO
    if(frozen_orb_scf)then
     call reorder_core_orb
     call initialize_mo_coef_begin_iteration
    endif

    TOUCH MO_coef

!   Calculate error vectors

    max_error_DIIS = maxval(Abs(FPS_SPF_Matrix_MO))

!   SCF energy

    energy_SCF = ueg_SCF_energy
    Delta_Energy_SCF = energy_SCF - energy_SCF_previous
    if ( (SCF_algorithm == 'DIIS').and.(Delta_Energy_SCF > 0.d0) ) then
      ueg_Fock_matrix_AO(1:ao_num,1:ao_num) = Fock_matrix_DIIS (1:ao_num,1:ao_num,index_dim_DIIS)
      ueg_Fock_matrix_AO_alpha = ueg_Fock_matrix_AO*0.5d0
      ueg_Fock_matrix_AO_beta  = ueg_Fock_matrix_AO*0.5d0
      TOUCH ueg_Fock_matrix_AO_alpha ueg_Fock_matrix_AO_beta
    endif

    double precision :: level_shift_save
    level_shift_save = level_shift
    mo_coef_save(1:ao_num,1:mo_num) = mo_coef(1:ao_num,1:mo_num)
    do while (Delta_energy_SCF > 0.d0)
      mo_coef(1:ao_num,1:mo_num) = mo_coef_save
      if (level_shift <= .1d0) then
        level_shift = 1.d0
      else
        level_shift = level_shift * 3.0d0
      endif
      TOUCH mo_coef level_shift
      mo_coef(1:ao_num,1:mo_num) = eigenvectors_Fock_matrix_MO(1:ao_num,1:mo_num)
      if(frozen_orb_scf)then
        call reorder_core_orb
        call initialize_mo_coef_begin_iteration
      endif
      TOUCH mo_coef
      Delta_Energy_SCF = ueg_SCF_energy - energy_SCF_previous
      energy_SCF = ueg_SCF_energy
      if (level_shift-level_shift_save > 40.d0) then
        level_shift = level_shift_save * 4.d0
        SOFT_TOUCH level_shift
        exit
      endif
      dim_DIIS=0
    enddo
    level_shift = level_shift * 0.5d0
    SOFT_TOUCH level_shift
    energy_SCF_previous = energy_SCF

    converged = ( (max_error_DIIS <= threshold_DIIS_nonzero) .and. &
                  (dabs(Delta_energy_SCF) <= thresh_SCF) )

!   Print results at the end of each iteration

    write(6,'(I4, 1X, F16.10, 1X, F16.10, 1X, F16.10, 1X, F16.10, 1X, I3)')  &
      iteration_SCF, energy_SCF, Delta_energy_SCF, max_error_DIIS, level_shift, dim_DIIS

!   Write data in JSON file

    call lock_io
    if (iteration_SCF == 1) then
      write(json_unit, json_dict_uopen_fmt)
    else
      write(json_unit, json_dict_close_uopen_fmt)
    endif
    write(json_unit, json_int_fmt) 'iteration', iteration_SCF
    write(json_unit, json_real_fmt) 'energy', energy_SCF
    write(json_unit, json_real_fmt) 'delta_energy_SCF', Delta_energy_SCF
    write(json_unit, json_real_fmt) 'max_error_DIIS', max_error_DIIS
    write(json_unit, json_real_fmt) 'level_shift', level_shift
    write(json_unit, json_int_fmt) 'dim_DIIS', dim_DIIS
    call unlock_io

    if (Delta_energy_SCF < 0.d0) then
      call save_mos
      write(json_unit, json_true_fmt) 'saved'
    else
      write(json_unit, json_false_fmt) 'saved'
    endif

    call lock_io
    if (converged) then
      write(json_unit, json_true_fmtx) 'converged'
    else
      write(json_unit, json_false_fmtx) 'converged'
    endif
    call unlock_io

    if (qp_stop()) exit

  enddo
  write(json_unit, json_dict_close_fmtx)

 if (iteration_SCF < n_it_SCF_max) then
   mo_label = 'Canonical'
 endif
!
! End of Main SCF loop
!

  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,*)
  if (converged) then
    write(6,*) 'SCF converged'
    call sleep(1)   ! When too fast, the MOs aren't saved.
  endif


  if(.not.frozen_orb_scf)then
   call mo_as_eigvectors_of_mo_matrix(ueg_Fock_matrix_mo,size(ueg_Fock_matrix_mo,1), &
      size(ueg_Fock_matrix_mo,2),mo_label,1,.true.)
   call restore_symmetry(ao_num, mo_num, mo_coef, size(mo_coef,1), 1.d-10)
   call orthonormalize_mos
  endif


  ! Identify degenerate MOs and force them on the axes
  allocate(S(ao_num,ao_num))
  i=1
  do while (i<mo_num)
    j=i
    m=1
    do while ( (j<mo_num).and.(ueg_fock_matrix_diag_mo(j+1)-ueg_fock_matrix_diag_mo(i) < 1.d-8) )
      j += 1
      m += 1
    enddo
    if (m>1) then
      call dgemm('N','T',ao_num,ao_num,m,1.d0,mo_coef(1,i),size(mo_coef,1),mo_coef(1,i),size(mo_coef,1),0.d0,S,size(S,1))
      call pivoted_cholesky( S, m, -1.d0, ao_num, mo_coef(1,i))
    endif
    i = j+1
  enddo


  call save_mos

  call write_double(6, Energy_SCF, 'SCF energy')

  call write_time(6)

end

