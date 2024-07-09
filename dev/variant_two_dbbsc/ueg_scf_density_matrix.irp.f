BEGIN_PROVIDER [double precision, ueg_SCF_density_matrix_ao_alpha, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! $C.C^t$ over $\alpha$ MOs
   END_DOC

   call dgemm('N','T',ao_num,ao_num,elec_alpha_num,1.d0, &
        mo_coef, size(mo_coef,1), &
        mo_coef, size(mo_coef,1), 0.d0, &
        ueg_SCF_density_matrix_ao_alpha, size(ueg_SCF_density_matrix_ao_alpha,1))

END_PROVIDER

BEGIN_PROVIDER [ double precision, ueg_SCF_density_matrix_ao_beta,  (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! $C.C^t$ over $\beta$ MOs
   END_DOC

   call dgemm('N','T',ao_num,ao_num,elec_beta_num,1.d0, &
        mo_coef, size(mo_coef,1), &
        mo_coef, size(mo_coef,1), 0.d0, &
        ueg_SCF_density_matrix_ao_beta, size(ueg_SCF_density_matrix_ao_beta,1))

END_PROVIDER

BEGIN_PROVIDER [ double precision, ueg_SCF_density_matrix_ao, (ao_num,ao_num) ]
   implicit none
   BEGIN_DOC
   ! Sum of $\alpha$ and $\beta$ density matrices
   END_DOC
   ASSERT (size(ueg_SCF_density_matrix_ao,1) == size(ueg_SCF_density_matrix_ao_alpha,1))
   if (elec_alpha_num== elec_beta_num) then
     ueg_SCF_density_matrix_ao = ueg_SCF_density_matrix_ao_alpha + ueg_SCF_density_matrix_ao_alpha
   else
     ASSERT (size(ueg_SCF_density_matrix_ao,1) == size(ueg_SCF_density_matrix_ao_beta ,1))
     ueg_SCF_density_matrix_ao = ueg_SCF_density_matrix_ao_alpha + ueg_SCF_density_matrix_ao_beta
   endif

END_PROVIDER
