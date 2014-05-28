program pressure_table
	use phys_constants
	use iso_fortran_env, only : error_unit
	use utils_lib
	use neutron_eos
	use skyrme
	real :: kn_min, kn_max
	real :: kn !fm^-1
	real :: dke, dkn
	real :: nn_cm, nn_fm
	real :: T, f_n, u_n, p_n, s_n, cv_n, dpr_n, dpT_n
	real :: eden_skyrme, eden_mb77, p_skyrme
	real :: dudn, n_prev, u_prev
	integer :: N_kn, j_kn, Nnuclides
	integer :: table_id, inlist_id, master_id
	character(len=*), parameter :: inputfile = 'pressure_table.in'
	character(len=80) :: table_name, master, file_prefix
	integer :: ierr

	namelist /init/ kn_min, kn_max, N_kn, file_prefix, master

	ierr = 0
	inlist_id = alloc_iounit(ierr)
	if (ierr /= 0) then
		write(error_unit,'(a)') 'unable to allocate unit'
		stop
	end if
	
	open(unit=inlist_id,file=inputfile,status='old',action='read',iostat=ierr)
	if (ierr /= 0) then
		write(error_unit,'(a)') 'unable to open inlist'
		stop
	end if
	
	read(inlist_id,nml=init,iostat = ierr)
	if (ierr /= 0) then
		write (error_unit,'(a)') 'unable to read inlist'
		stop
	end if
	close(inlist_id)

	!create master file and allocate master_id
	write(master,'(a)') trim(file_prefix)//'mb77_pressure_master.data'
	master_id = alloc_iounit(ierr)
	if (ierr /= 0) then
		write (error_unit,'(a)') 'unable to allocate master file unit no.'
		stop
	end if
	
	!open master file
	open(master_id,file=trim(master),status='unknown',action='write',iostat=ierr)
	if (ierr /= 0) then
		write(error_unit,'(a)') 'unable to open master file'
		stop
	end if
	
	write(master_id,'(a5)') 'N(kn)'
	write(master_id,'(i5)') N_kn
	write(master_id, '((a5,tr1),8(tr1,a13))') 'j(kn)', 'kn [fm**-1]', &
	             & 'nn [cm**-3]', 'nn [fm**-3]', 'u_s [MeV/cc]', &
                     & 'u_m [MeV/cc]', 'dudn [fm**-1]', 'ps[dyn/cm**2]'&
		     & ,'p [dyn/cm**2]'

	dkn = (kn_max-kn_min)/real(N_kn-1)

	u_prev = 0.0
	n_prev = 0.0 
	dudn = 0.0 
	
	do j_kn = 1,N_kn

			kn = kn_min+dkn*real(j_kn-1)

			!find skyrme energy density
			nn_fm = 2.0*kn**3/threepisquare
			nn_cm = nn_fm*fm_to_cm**(-3.0) 
			call skyrme_eos(nn_fm,0.0,eden_skyrme) 
			
			!find skyrme pressure
			if(nn_fm .ne. 0.0) then
                        dudn = (eden_skyrme-u_prev)/(nn_fm-n_prev)  !fm**-1
	                u_prev = eden_skyrme
	                n_prev = nn_fm
			endif
			p_skyrme = nn_fm*dudn-eden_skyrme           !fm**-4
			p_skyrme = p_skyrme*fm_to_cm**(-3.0)/(5.06773d-3)
			p_skyrme = p_skyrme*mev_to_ergs

			!find mb77 pressure and energy density		
			call MB77(nn_cm,T,f_n,u_n,p_n,s_n,cv_n,dpr_n,dpt_n)

			!put in units of MeV/cc
			eden_mb77 = u_n*nn_cm*ergs_to_mev
			eden_skyrme = eden_skyrme*fm_to_cm**(-3.0)/(5.06773d-3)

			if (ierr == 0) then				   
			write(master_id,'((i5,tr1),8(tr1,es13.7))')  &
				& j_kn, kn, nn_cm, nn_fm, eden_skyrme, &
				& eden_mb77, dudn, p_skyrme, p_n
			else	! if an error occurs just reset ierr and move on
				write (error_unit,'(a)') 'error in writing file...moving on'
				ierr = 0
			end if


 	enddo

	close(master_id)
	call free_iounit(master_id)

end program pressure_table
