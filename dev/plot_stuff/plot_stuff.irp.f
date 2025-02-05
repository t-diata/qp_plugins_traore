program plot_stuff_prog
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
  double precision :: r(3), dx, dm_a, dm_b
  
  integer :: i, j, k, l, nx, istate    
  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
  ! File to write the density
  output=trim(ezfio_filename)//'.density'
  i_unit_output = getUnitAndOpen(output,'w')

  r(:) = nucl_coord_transp(:,1) ! starting point

  nx = 100 ! Number of points
  dx = 3. / nx
  l = 0
  write(i_unit_output,*)' #  x   y   z      density'
  do i = 1, nx
      r(2) = nucl_coord_transp(2,1)
	do j = 1, nx
		r(3) = nucl_coord_transp(3,1)
		do k = 1, nx
 			call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
			write(i_unit_output,'(100(F16.10,X))') r(1), r(2), r(3), dm_a+dm_b
			!write(i_unit_output,*) r(1), r(2), r(3), dm_a+dm_b
			r(3) += dx
			l+=1
		enddo
		r(2)+=dx
		l+=1
	enddo
	r(1) +=dx
	l+=1
	print *, "Done: ", 100*float(l)/float(nx**3), "%"
  enddo

end
