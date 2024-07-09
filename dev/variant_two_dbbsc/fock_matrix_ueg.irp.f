 BEGIN_PROVIDER [ double precision, ueg_Fock_matrix_ao_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ueg_Fock_matrix_ao_beta,  (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! copy of src/hartree-fock/fock_matrix_hf.irp.f
 ! Alpha Fock matrix in AO basis set
 END_DOC

 integer                        :: i,j
 do j=1,ao_num
   do i=1,ao_num
     ueg_Fock_matrix_ao_alpha(i,j) = ueg_ao_one_e_integrals(i,j) + ao_two_e_integral_alpha(i,j)
     ueg_Fock_matrix_ao_beta (i,j) = ueg_ao_one_e_integrals(i,j) + ao_two_e_integral_beta (i,j)

   enddo
 enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision, ueg_ao_one_e_integrals,(ao_num,ao_num)]
&BEGIN_PROVIDER [ double precision, ueg_ao_one_e_integrals_diag,(ao_num)]
  implicit none
  integer :: i,j,n,l
  BEGIN_DOC
 ! copy from src/ao_one_e_ints/ao_one_e_int.irp.f
 ! One-electron Hamiltonian in the |AO| basis.
  END_DOC

  IF (read_ao_one_e_integrals) THEN
     call ezfio_get_ao_one_e_ints_ao_one_e_integrals(ueg_ao_one_e_integrals)
  ELSE
        !ao_one_e_integrals = ao_integrals_n_e + ao_kinetic_integrals
        ueg_ao_one_e_integrals = ao_kinetic_integrals

  ENDIF

  DO j = 1, ao_num
    ao_one_e_integrals_diag(j) = ueg_ao_one_e_integrals(j,j)
    ueg_ao_one_e_integrals_diag(j) = ueg_ao_one_e_integrals(j,j)
  ENDDO

  IF (write_ao_one_e_integrals) THEN
       call ezfio_set_ao_one_e_ints_ao_one_e_integrals(ueg_ao_one_e_integrals)
       print *,  'ueg AO one-e integrals written to disk'
  ENDIF

END_PROVIDER

 BEGIN_PROVIDER [ double precision, ueg_Fock_matrix_mo, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, ueg_Fock_matrix_diag_mo, (mo_num)]
   implicit none
   BEGIN_DOC
   ! copy from src/hartree_fock/scf_utils
   ! Fock matrix on the MO basis.
   ! For open shells, the ROHF Fock Matrix is ::
   !
   !       |  Rcc  |  F^b  |  Fcv  |
   !       |-----------------------|
   !       |  F^b  |  Roo  |  F^a  |
   !       |-----------------------|
   !       |  Fcv  |  F^a  |  Rvv  |
   !
   ! C: Core, O: Open, V: Virtual
   !
   ! Rcc = Acc Fcc^a + Bcc Fcc^b
   ! Roo = Aoo Foo^a + Boo Foo^b
   ! Rvv = Avv Fvv^a + Bvv Fvv^b
   ! Fcv = (F^a + F^b)/2
   !
   ! F^a: Fock matrix alpha (MO), F^b: Fock matrix beta (MO)
   ! A,B: Coupling parameters
   !
   ! J. Chem. Phys. 133, 141102 (2010), https://doi.org/10.1063/1.3503173
   ! Coupling parameters from J. Chem. Phys. 125, 204110 (2006); https://doi.org/10.1063/1.2393223.
   !       cc   oo   vv
   !  A  -0.5  0.5  1.5
   !  B   1.5  0.5 -0.5
   !
   END_DOC
   integer                        :: i,j,n
   if (all_shells_closed) then
     ueg_Fock_matrix_mo = ueg_Fock_matrix_mo_alpha
   else
     ! Core
     do j = 1, elec_beta_num
       ! Core
       do i = 1, elec_beta_num
         ueg_fock_matrix_mo(i,j) = - 0.5d0 * ueg_fock_matrix_mo_alpha(i,j) &
                               + 1.5d0 * ueg_fock_matrix_mo_beta(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         ueg_fock_matrix_mo(i,j) = ueg_fock_matrix_mo_beta(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         ueg_fock_matrix_mo(i,j) =   0.5d0 * ueg_fock_matrix_mo_alpha(i,j) &
                               + 0.5d0 * ueg_fock_matrix_mo_beta(i,j)
       enddo
     enddo
     ! Open
     do j = elec_beta_num+1, elec_alpha_num
       ! Core
       do i = 1, elec_beta_num
         ueg_fock_matrix_mo(i,j) = ueg_fock_matrix_mo_beta(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         ueg_fock_matrix_mo(i,j) =   0.5d0 * ueg_fock_matrix_mo_alpha(i,j) &
                               + 0.5d0 * ueg_fock_matrix_mo_beta(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         ueg_fock_matrix_mo(i,j) = ueg_fock_matrix_mo_alpha(i,j)
       enddo
     enddo
     ! Virtual
     do j = elec_alpha_num+1, mo_num
       ! Core
       do i = 1, elec_beta_num
         ueg_fock_matrix_mo(i,j) =   0.5d0 * ueg_fock_matrix_mo_alpha(i,j) &
                               + 0.5d0 * ueg_fock_matrix_mo_beta(i,j)
       enddo
       ! Open
       do i = elec_beta_num+1, elec_alpha_num
         ueg_fock_matrix_mo(i,j) = ueg_fock_matrix_mo_alpha(i,j)
       enddo
       ! Virtual
       do i = elec_alpha_num+1, mo_num
         ueg_fock_matrix_mo(i,j) =   1.5d0 * ueg_fock_matrix_mo_alpha(i,j) &
                               - 0.5d0 * ueg_fock_matrix_mo_beta(i,j)
       enddo
     enddo
   endif

   do i = 1, mo_num
     Fock_matrix_diag_mo(i) = Fock_matrix_mo(i,i)
   enddo


   if(frozen_orb_scf)then
     integer                        :: iorb,jorb
     !       active|core|active
     !active |     | 0  |      
     !core   |  0  |    |   0
     !active |     | 0  |       
     do i = 1, n_core_orb
      iorb = list_core(i)
      do j = 1, n_act_orb
       jorb = list_act(j)
       ueg_Fock_matrix_mo(iorb,jorb) = 0.d0
       ueg_Fock_matrix_mo(jorb,iorb) = 0.d0
      enddo
     enddo
   endif

   if(no_oa_or_av_opt)then
     do i = 1, n_act_orb
       iorb = list_act(i)
       do j = 1, n_inact_orb
         jorb = list_inact(j)
         ueg_Fock_matrix_mo(iorb,jorb) = 0.d0
         ueg_Fock_matrix_mo(jorb,iorb) = 0.d0
       enddo
       do j = 1, n_virt_orb
         jorb = list_virt(j)
         ueg_Fock_matrix_mo(iorb,jorb) = 0.d0
         ueg_Fock_matrix_mo(jorb,iorb) = 0.d0
       enddo
       do j = 1, n_core_orb
         jorb = list_core(j)
         ueg_Fock_matrix_mo(iorb,jorb) = 0.d0
         ueg_Fock_matrix_mo(jorb,iorb) = 0.d0
       enddo
     enddo
   endif

END_PROVIDER

BEGIN_PROVIDER [ double precision, ueg_Fock_matrix_mo_alpha, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! copy from src/hartree_fock/scf_utils
   ! Fock matrix on the MO basis
   END_DOC
   call ao_to_mo(ueg_Fock_matrix_ao_alpha,size(ueg_Fock_matrix_ao_alpha,1), &
                 ueg_Fock_matrix_mo_alpha,size(ueg_Fock_matrix_mo_alpha,1))
END_PROVIDER


BEGIN_PROVIDER [ double precision, ueg_Fock_matrix_mo_beta, (mo_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! copy from src/hartree_fock/scf_utils
   ! Fock matrix on the MO basis
   END_DOC
   call ao_to_mo(ueg_Fock_matrix_ao_beta,size(ueg_Fock_matrix_ao_beta,1), &
                 ueg_Fock_matrix_mo_beta,size(ueg_Fock_matrix_mo_beta,1))
END_PROVIDER


BEGIN_PROVIDER [ double precision, ueg_Fock_matrix_ao, (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! copy from src/hartree_fock/scf_utils
 ! Fock matrix in AO basis set
 END_DOC

 if(frozen_orb_scf)then
   call mo_to_ao(ueg_Fock_matrix_mo,size(ueg_Fock_matrix_mo,1),              &
       ueg_Fock_matrix_ao,size(ueg_Fock_matrix_ao,1))
 else
   if (all_shells_closed.and. (level_shift == 0.)) then
     integer                        :: i,j
     do j=1,ao_num
       do i=1,ao_num
         ueg_Fock_matrix_ao(i,j) = ueg_Fock_matrix_ao_alpha(i,j)
       enddo
     enddo
   else
     call mo_to_ao(ueg_Fock_matrix_mo,size(ueg_Fock_matrix_mo,1),            &
         ueg_Fock_matrix_ao,size(ueg_Fock_matrix_ao,1))
   endif
 endif
END_PROVIDER


BEGIN_PROVIDER [ double precision, ueg_SCF_energy ]
 implicit none
 BEGIN_DOC
 ! ueg Hartree-Fock energy
 END_DOC
 ueg_SCF_energy = nuclear_repulsion

 integer                        :: i,j
 do j=1,ao_num
   do i=1,ao_num
     ueg_SCF_energy += 0.5d0 * (                                          &
         (ueg_ao_one_e_integrals(i,j) + ueg_Fock_matrix_ao_alpha(i,j) ) *  ueg_SCF_density_matrix_ao_alpha(i,j) +&
         (ueg_ao_one_e_integrals(i,j) + ueg_Fock_matrix_ao_beta (i,j) ) *  ueg_SCF_density_matrix_ao_beta (i,j) )
   enddo
 enddo
 ueg_SCF_energy += extra_e_contrib_density

END_PROVIDER

