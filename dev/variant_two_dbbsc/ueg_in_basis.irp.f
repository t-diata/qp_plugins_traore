BEGIN_PROVIDER [double precision, ueg_HF_energy]
 implicit none
 BEGIN_DOC
 ! Copy of src/hf_energy.irp.f
 ! Hartree-Fock energy containing the nuclear repulsion, and its one-body components.
 END_DOC

 integer :: i,j
 ueg_HF_energy = nuclear_repulsion
 do j=1,ao_num
   do i=1,ao_num
    ueg_HF_energy += ao_one_e_integrals(i,j) * (SCF_density_matrix_ao_alpha(i,j) + SCF_density_matrix_ao_beta (i,j) )
   enddo
 enddo

END_PROVIDER

