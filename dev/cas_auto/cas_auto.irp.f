      program cas_auto
      implicit none
      BEGIN_DOC
      ! TODO : Put the documentation of the program here
      END_DOC
       read_wf = .true.
       touch read_wf
       ! total one-e integrals
       io_mo_one_e_integrals = "None"
       touch io_mo_one_e_integrals
       ! Vne integrals on the MO basis
       io_mo_integrals_n_e = "None"
       touch io_mo_integrals_n_e
       ! kinetic integrals on the MO basis
       io_mo_integrals_kinetic = "None"
       touch io_mo_integrals_kinetic
       ! Vne integrals on the AO basis
       io_ao_integrals_n_e = "None"
       touch io_ao_integrals_n_e
       ! kinetic integrals on the AO basis
       io_ao_integrals_kinetic = "None"
       touch io_ao_integrals_kinetic
 
       print*, 'Entanglement measurement'
       print*, s_entanglement

       end


 BEGIN_PROVIDER [double precision, s_entanglement, (n_act_orb)]
      implicit none
      BEGIN_DOC
      ! Entropy based entanglement measure according to "Automated
      ! Selection of Active Spaces", Stein & Reiher, JCTC 2016, 12
      ! 176--1771, Eq. 3
      ! The CI version is provided by "Entanglement measure of
      ! Electron-Electron Correlation in Quantum Chemistry
      ! Calculations", Huang & Kais, CPL 2005, 413, 1-5, Eq. 11
      END_DOC
      integer :: i,j

      ! Compute single-orbital entropy
      ! s_i = \omega_alpha ln \omega_alpha
      do i=1, n_act_orb
        if (eig_for_measure(i) < 1e-10) then
                s_entanglement(i) = 0.d0
        else
                s_entanglement(i) = eig_for_measure(i)*dlog(eig_for_measure(i))/dlog(2.d0)
        end if
      enddo


 END_PROVIDER

 BEGIN_PROVIDER [double precision, s_two_orb_entanglement, (n_act_orb)]
      implicit none
      BEGIN_DOC
      ! Entropy based entanglement measure according to "Automated
      ! Selection of Active Spaces", Stein & Reiher, JCTC 2016, 12
      ! 176--1771, Eq. 4
      END_DOC
      integer :: i,j

      ! Compute two-orbital entropy
      ! s_ij = \omega_alpha ln \omega_alpha
      do i=1, n_act_orb
        if (eig_for_measure(i) < 1e-10) then
                s_entanglement(i) = 0.d0
        else
                s_entanglement(i) = eig_for_measure(i)*dlog(eig_for_measure(i))/dlog(2.d0)
        end if
      enddo


 END_PROVIDER


 BEGIN_PROVIDER [double precision, eig_for_measure, (n_act_orb)]
      implicit none
      BEGIN_DOC
      ! eigenvalues of the density matrix for the entanglement measure
      END_DOC

      integer :: i, j, info, lwork
      double precision, allocatable :: d_for_diag(:,:), w_tab(:)
      double precision, allocatable :: eig_desc(:)
      allocate(eig_desc(n_act_orb))

      lwork = max(1, 3*n_act_orb-1)
      allocate(d_for_diag(n_act_orb,n_act_orb), w_tab(lwork))
      
      d_for_diag = D0tu

      ! diagonalize the density matrix
      call dsyev('N', 'U', n_act_orb, d_for_diag, n_act_orb, eig_for_measure, w_tab, lwork, info)

      do i=0, n_act_orb-1
       j=n_act_orb-i
       eig_desc(i+1) = eig_for_measure(j)
      enddo

      eig_for_measure = eig_desc
      deallocate(d_for_diag, w_tab)

      deallocate(eig_desc)
 END_PROVIDER

 BEGIN_PROVIDER [double precision, twordm_for_measure, (n_act_orb)]
      implicit none
      BEGIN_DOC
      ! eigenvalues of the density matrix for the entanglement measure
      END_DOC

      integer :: i, j, info, lwork
      double precision, allocatable :: d_for_diag(:,:), w_tab(:)
      double precision, allocatable :: eig_desc(:)
      allocate(eig_desc(n_act_orb))

      lwork = max(1, 3*n_act_orb-1)
      allocate(d_for_diag(n_act_orb,n_act_orb), w_tab(lwork))
      
      d_for_diag = D0tu

      ! diagonalize the density matrix
      call dsyev('N', 'U', n_act_orb, d_for_diag, n_act_orb, eig_for_measure, w_tab, lwork, info)

      do i=0, n_act_orb-1
       j=n_act_orb-i
       eig_desc(i+1) = eig_for_measure(j)
      enddo

      eig_for_measure = eig_desc
      deallocate(d_for_diag, w_tab)

      deallocate(eig_desc)
 END_PROVIDER
