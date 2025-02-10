program plot_stuff_3d_prog
!  BEGIN_DOC
! TODO : Put the documentation of the program here
!  END_DOC

  implicit none
  read_wf = .True.
  touch read_wf
  call routine

end

subroutine routine
  implicit none
  double precision :: r(3), dx, dm_a, dm_b, e_PBE !,eps_c_md_PBE
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  double precision :: g0_UEG_mu_inf, denom, rho_tot, grad_rho_tot
  double precision :: f_ii_val_ab_custom,two_bod_dens_ii,f_ia_val_ab_custom,two_bod_dens_ia,f_aa_val_ab_custom,two_bod_dens_aa,f_psi_cas_ab_custom(N_states)


  
  integer :: i, j, k, l, m, nx, istate    
  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
  ! File to write the density
  output=trim(ezfio_filename)//'.density'
  i_unit_output = getUnitAndOpen(output,'w')

  r(:) = nucl_coord_transp(:,1) ! starting point

  nx = 100 ! Number of points
  dx = 3. / nx
  l = 0

  write(i_unit_output,*)'x   y   z  hf  density density_gradient fpsi ec_pbe'
  do i = 1, nx
	  r(2) = nucl_coord_transp(2,1)
	  do j = 1, nx
		r(3) = nucl_coord_transp(3,1)
		do k = 1, nx
			!call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
			call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r, rho_a, rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
	        	grad_rho_a_2 = 0.d0
                	grad_rho_b_2 = 0.d0
                	grad_rho_a_b = 0.d0
                	do istate = 1, N_states
                 		do m = 1, 3
                  			grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
                  			grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
                  			grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
                 		enddo
                	enddo
               		do istate = 1, N_states
                		! convertion from (alpha,beta) formalism to (closed, open) formalism
                		call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
                		call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
                		call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE)

				! inactive-inactive part of f_psi(r1,r2)
				call give_f_ii_val_ab(r,r,f_ii_val_ab_custom,two_bod_dens_ii)
   				! inactive-active part of f_psi(r1,r2)
   				call give_f_ia_val_ab(r,r,f_ia_val_ab_custom,two_bod_dens_ia,istate)
   				! active-active part of f_psi(r1,r2)
   				call give_f_aa_val_ab(r,r,f_aa_val_ab_custom,two_bod_dens_aa,istate)
   				f_psi_cas_ab_custom(istate)    = f_ii_val_ab_custom     + f_ia_val_ab_custom     + f_aa_val_ab_custom
	       		enddo

	       		istate = 1
	       		!call ecmd_pbe_ueg_at_r(0.5, r, eps_c_md_PBE)
	       		grad_rho_tot = sqrt((grad_rho_a(1,1)+grad_rho_b(1,1))**2 + (grad_rho_a(2,1)+grad_rho_b(2,1))**2 + (grad_rho_a(3,1)+grad_rho_b(3,1))**2 ) 
	       		write(i_unit_output,'(100(F16.10,X))') r(1), r(2), r(3), HF_energy, rho_a(1)+rho_b(1), grad_rho_tot, f_psi_cas_ab_custom(istate), e_PBE  
	       		r(3) += dx
	       		l+=1
		enddo
		r(2)+=dx
		l+=1
  	enddo
	r(1)+=dx
	l+=1
         	print *, "Done: ", 100*float(l)/float(nx**3), "%"
	 enddo

end
