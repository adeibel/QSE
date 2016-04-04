program follow_chain
		use iso_fortran_env, only : error_unit, output_unit
		use phys_constants
		use utils_lib
		use mb77, only : neutron_chemical_potential, electron_pressure, neutron_pressure
		use mass_table
		use stable_table
		use alert_lib
		use rootfind	
	
		character(len=*), parameter :: default_infile = 'nucchem_chain.inlist'
		character(len=*), parameter :: default_mass_table = 'moe95_converted.data'
		character(len=*), parameter :: default_stable_table = 'stable_nuclei_moe95.data'
		character(len=*), parameter :: final_file = 'final_array.data'
		real :: kn, ke, mu_n
		real :: mu_e, ne, nn 
		real :: x1, x2, xacc
		real :: n, pres_n, pressure
		real :: pressure_start, pressure_stop, pressure_increment
		real :: mu_e_start, mu_e_stop, mu_e_increment
		real :: Sn, S2n, Sp, S2p, B, ecthresh, bthresh, VN
		real :: Snr, S2nr, Spr, S2pr, Br, ecthreshr, bthreshr, VNr
		real :: alpha(2), beta(2), gamma(2), delta(2), epsilon(2)
		real :: mu_n_range
		real :: A_sum
		real :: m_term, m_nuc
		real :: mu_n_new, n_n_new
		real, dimension(:), allocatable :: xmass
		integer, dimension(:), allocatable :: Z_int, A_int
		integer, dimension(:), allocatable :: Z_fin, A_fin 
		integer :: i, j, i_enter
		integer :: Z, A, Zr, Ar, inlist_id
		integer :: ierr, id, ios, iter, iZ, iZb, iZe, iEq(1)
		integer :: k, l, final_id
		integer, parameter :: ineg = 1, ipos = 2
		integer, parameter :: fid = output_unit, max_iterations = 100
		real, parameter :: del_m = mn_n-mp_n-me_n
		real, parameter :: mun_max = 6.5
		character(len=256) :: arg, mass_table_used, pfile, mass_table_name
		character(len=256) :: stable_table_used
		character(len=6) :: rxn
		real,dimension(:),pointer :: hist
		type(mass_table_type), pointer :: mt
		type(stable_table_type), pointer :: st
		logical, dimension(max_iterations) :: neutron_capture, dineutron_capture
		logical, dimension(max_iterations) :: neutron_emission, dineutron_emission
		logical, dimension(max_iterations) :: en_rxn, enn_rxn, ne_rxn, nne_rxn
		logical, save :: stable_table_loaded = .FALSE.
		logical, save :: mass_tables_are_loaded = .FALSE.
		logical, save :: pressure_set = .FALSE.
	
		namelist /io/ mass_table_name
		namelist /range/ Z, A, mu_e_start, mu_e_stop, &
					pressure_start, pressure_stop
				
	   ! set defaults in MeV
	   pfile = default_infile
	   mass_table_used = default_mass_table
	   stable_table_used = default_stable_table
	   Z = 26
	   A = 56
	   mu_e_start = 0.0
	   mu_e_stop = 30.0
	   mu_n = 0.0
	
	
	   ! read in the inputs
	   ierr = 0
	   inlist_id = alloc_iounit(ierr)
	   if (io_failure(ierr,'allocating unit for namelist')) stop
   
	   open(unit=inlist_id, file=pfile, iostat=ios, status="old", action="read")
	   if (io_failure(ios,'opening inlist file'//trim(pfile))) stop
   
	   read(inlist_id,nml=io,iostat=ios)
	   if (io_failure(ios,'unable to read namelist "io"')) stop
	   read(inlist_id,nml=range,iostat=ios)
	   if (io_failure(ios,'unable to read namelist "range"')) stop
	   close(inlist_id)
	   call free_iounit(inlist_id)
   
		! load the mass table
		if (mass_tables_are_loaded .eqv. .FALSE.) then
		call load_mass_table('../../../data',trim(mass_table_used),ierr)
		mass_tables_are_loaded = .TRUE. 
		if (ierr /= 0) then
			write (error_unit,'(a)') 'mass table loading error'
			stop
		end if
		end if
		mt => winvn_mass_table
	
	   ! check that we are on the table
	   if (Z < mt% Zmin .or. Z > mt% Zmax) then
		  write(error_unit,'(a,"[",2i4,"]")') 'Z must be in table range ',mt% Zmin,mt% Zmax
		  stop
	   end if
   
	   ! load stable nuclei table
	   if (stable_table_loaded .eqv. .FALSE.) then
	   !call load_stable_table('../../../data',trim(stable_table_used),ierr)
	   call load_stable_table('../../../data',trim(stable_table_used),ierr)
	   stable_table_loaded = .TRUE.
	   if (ierr /= 0) then
			write (error_unit,'(a)') 'stable table loading error'
			stop
		end if
		end if
		st => winvn_stable_table
	
	   final_id = alloc_iounit(ierr)
	   if (io_failure(ierr,'allocating unit for final array file')) stop
	   open(unit=final_id, file=final_file, iostat=ios, status='unknown') 
	   if (io_failure(ios,'opening final array file')) stop
	   write(final_id,'(3(a10,2x),4(A10,2x))') 'Pressure', 'mu_e', 'mu_n', 'Zint', 'Aint', &
				& 'Zfin', 'Afin'
   
  
	   !main loop over pressure (pushes nucleus to higher pressures)   
   
		! allocate size of arrays equal to the distribution of nuclei compressed 
   		allocate(Z_int(size(Z)), A_int(size(A)), Z_fin(size(Z)), A_fin(size(A)))
   
	   do i = 1,1000
		pressure = pressure_start*real(i) ! MeV fm**-3
	
		! end of last cascade is moved up in pressure 
		if (i >= 1 .and. Z_fin(1) > 0) then
		Z_int = Z_fin
		A_int = A_fin
		end if

		! scan through all mu_e and mu_n for given pressure
		!do j = 359,1196 !2699, 7771 !1, st% Ntable
		do j = 1,1	
	
		 ! if (st% P(j) /= pressure) cycle
	
		  if (st% P(j) == pressure) then
		  mu_e = st% mu_e(j)
		  mu_n = st% mu_n(j)
		  end if

		 mu_e = 30.383329
		 mu_n = 2.4670114897d-8
	
		 !initial nuclei accreted to the given pressure 
		 do k = 1, size(Z_int)
		  Z = Z_int(k)
		  A = A_int(k)
	 
		  neutron_capture = .FALSE.
		  neutron_emission = .FALSE.
		  dineutron_capture = .FALSE.
		  dineutron_emission = .FALSE.
		  en_rxn = .FALSE.
		  enn_rxn = .FALSE. 
		  ne_rxn = .FALSE.
		  nne_rxn = .FALSE.  
	 
		  ! loop over reactions
		  do iter = 1, max_iterations
		  !get properties of nucleus that is being pushed deeper
		  call get_nucleus_properties(Z,A,id,B,Sn,S2n,Sp,S2p,ecthresh,bthresh,VN,ierr)
	 
		  ! check for weak reactions                 
		  ! electron capture
		  Ar = A; Zr = Z-1
		  call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr, &
				&	ecthreshr,bthreshr,VNr,ierr)
		  if (ierr /= 0) then
			 write(error_unit,'(a)') 'unable to find nucleus'
			 stop
			 exit
		  end if
		  alpha(ineg) = del_m - mu_e + (B-Br)
		  if (alpha(ineg) < 0) then
			 rxn = '(e,)'
			 call print_reaction
			 A = Ar; Z = Zr
			 cycle
		  end if
	  
		  ! electron emission
		  Ar = A; Zr = Z+1
		  call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
		  if (ierr /= 0) then
			 write(error_unit,'(a)') 'unable to find nucleus'
			 stop
			 exit
		  end if
		  alpha(ipos) = mu_e - del_m + (B - Br)
		  if (alpha(ipos) < 0) then
			 rxn = '(,e)'
			 call print_reaction
			 A = Ar; Z = Zr
			 cycle
		  end if
	  
		  ! electron capture followed by neutron emission 
		  Ar = A-1; Zr = Z-1
		  call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
		  if (ierr /= 0) then
			 write(error_unit,'(a)') 'unable to find nucleus'
			 stop
			 exit
		  end if
		  beta(ineg) = mu_n - mu_e + del_m + (B - Br)
		  if (beta(ineg) < 0) then
			 en_rxn(iter) = .TRUE.
			 rxn = '(e,n)'
			 call print_reaction
			 A = Ar; Z = Zr
			 cycle
		  end if
			
		  ! electron capture followed by dineutron emission
		  Ar = A-2; Zr = Z-1
		  call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
		  if (ierr /= 0) then
			 write(error_unit,'(a)') 'unable to find nucleus'
			 stop
			 exit
		  end if
		  beta(ipos) = 2*mu_n - mu_e + del_m + (B - Br)
		  if (beta(ipos) < 0) then
			 enn_rxn(iter) = .TRUE.
			 rxn = '(e,2n)'
			 call print_reaction
			 A = Ar; Z = Zr
			 cycle
		  end if
	  
		  ! check for strong reactions     
		  ! neutron capture
		  Ar = A+1; Zr = Z
		  call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
		  if (ierr /= 0) then
			 write(error_unit,'(a)') 'unable to find nucleus'
			 stop
			 exit
		  end if
		  gamma(ineg) = -mu_n + (B-Br)
		  if (gamma(ineg) < 0) then
			 neutron_capture(iter) = .TRUE.
			 rxn = '(n,)'
			 call print_reaction
			 A = Ar; Z = Zr
			 cycle
		  end if
	  
		  ! neutron emission
		  Ar = A-1; Zr = Z
		  call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
		  if (ierr /= 0) then
			 write(error_unit,'(a)') 'unable to find nucleus'
			 stop
			 exit
		  end if
		  gamma(ipos) = mu_n + (B-Br)
		  if (gamma(ipos) < 0) then
			 neutron_emission(iter) = .TRUE.
			 rxn = '(,n)'
			 call print_reaction
			 A = Ar; Z = Zr
			 cycle
		  end if
	  
		  ! dineutron capture
		  Ar = A+2; Zr = Z
		  call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
		  if (ierr /= 0) then
			 write(error_unit,'(a)') 'unable to find nucleus'
			 stop
			 exit
		  end if
		  delta(ineg) = -2.0*mu_n + (B-Br)
		  if (delta(ineg) < 0) then
			 dineutron_capture(iter) = .TRUE.
			 rxn = '(2n,)'
			 call print_reaction
			 A = Ar; Z = Zr
			 cycle
		  end if
	  
		  ! dineutron emission
		  Ar = A-2; Zr = Z
		  call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
		  if (ierr /= 0) then
			 write(error_unit,'(a)') 'unable to find nucleus'
			 stop
			 exit
		  end if
		  delta(ipos) = 2.0*mu_n + (B-Br)
		  if (delta(ipos) < 0) then
			 dineutron_emission(iter) = .TRUE.
			 rxn = '(,2n)'
			 call print_reaction
			 A = Ar; Z = Zr
			 cycle
		  end if      
		 
		  ! neutron capture followed by electron emission
		  Ar = A+1; Zr = Z+1
		  call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
		  if (ierr /= 0) then
			 write(error_unit,'(a)') 'unable to find nucleus'
			 stop
			 exit
		  end if
		  epsilon(ineg) = mu_e - mu_n - del_m + (B-Br)
		  if (epsilon(ineg) < 0) then
			 ne_rxn(iter) = .TRUE.
			 rxn = '(n,e)'
			 call print_reaction
			 A = Ar; Z = Zr
			 cycle
		  end if
	  
		  !dineutron capture followed by electron emission
		  Ar = A+2; Zr = Z+1
		  call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,Vnr,ierr)
		  if (ierr /= 0) then
			 write(error_unit, '(a)') 'unable to find nucleus'
			 stop
			 exit
		  end if
		  epsilon(ipos) = mu_e - 2.0*mu_n - del_m + (B-Br)
		  if (epsilon(ipos) < 0) then
			 nne_rxn(iter) = .TRUE.
			 rxn = '(2n,e)'
			 call print_reaction
			 A = Ar; Z = Zr
			 cycle
		  end if 

		  ! end of nucleus stability check 
		  Z_fin(k) = Z
		  A_fin(k) = A

		  end do !end iterations over reactions

		  A_sum = 0.
		  ! size of xmass array should be equal to the number
		  !  of nuclei in the distribution being compressed
          allocate(xmass(size(Z_fin)))

		  ! begin check for mass density increase
	  	  if (Z_int(k) .ne. Z_fin(k)) then
	  	   !recalculate mu_n if strong reaction goes
		   do l = 1, size(Z_int)
            !number density of isotopes
		    m_nuc = real(mt% A(i))*amu_n       
     	    m_term = g*(m_nuc*kT/(twopi*hbarc_n**2))**(1.5)
		    xmass(i) = real(mt% A(i))*m_term*exp((mu_i(i)+mt%BE(i))/kT)/n_b		   	
	  	    A_sum = A_sum + xmass(i)
	  	   end do
	  	  n_n_new = n_b - A_sum  
	  	  x1=0.0
          x2=10.0
          xacc=1.d-15
	  	  mu_n_new = root_bisection(neutron_chemical_potential,x1,x2,xacc,ierr,hist) 
          if (io_failure(ierr,'error in bisection for mu_n')) then
          stop
          end if
	      endif ! end check for mass density increase 

		  ! begin check for baryon conservation	  	  	  
		  if (count(neutron_capture) /= count(neutron_emission)) then
		  write(*,*) count(neutron_capture), count(neutron_emission), Z, A
		  if (count(neutron_capture) > count(neutron_emission)) then

		  !scan for neutron emissions from initial nuclei array
			do l = 1, size(Z_int)
			 A = A_int(l) ; Z = Z_int(l)
			 Ar = A-1; Zr = Z 
			 call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
			 if (ierr /= 0) then
			  write(error_unit,'(a)') 'unable to find nucleus'
			  stop
			  exit
			 cycle
			 end if
			 gamma(ipos) = mu_n + (B-Br)
			 if (gamma(ipos) < 0) then
			  neutron_emission(iter) = .TRUE.
			  rxn = '(,n)'
			  call print_reaction_check
			  A = Ar; Z = Zr
			  !stop
			 end if
		   end do
		  else
		 !scan for neutron captures from initial nuclei array
		   do l = 1, size(Z_int)
			A = A_int(l) ; Z = Z_int(l)
			Ar = A+1 ; Zr = Z
			call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
			if (ierr /= 0) then
			 write(error_unit,'(a)') 'unable to find nucleus'
			 stop
			 exit
			end if
			gamma(ineg) = -mu_n + (B-Br)
			if (gamma(ineg) < 0) then
			 neutron_capture(iter) = .TRUE.
			 rxn = '(n,)'
			 call print_reaction_check
			 A = Ar; Z = Zr
			 !stop
			end if
		   enddo           
		  end if
		  end if
	  
		  if (count(dineutron_capture) /= count(dineutron_emission)) then	  
		  if (count(dineutron_capture) > count(dineutron_emission)) then
		  !scan for dineutron emissions from initial nuclei array
		  else
		  !scan for dineutron captures from initial nuclei array
		  end if
		  end if

		  if (count(en_rxn) /= count(ne_rxn)) then
		  if (count(en_rxn) > count(ne_rxn)) then
		  !scan back for ne_rxn from initial nuclei array 
		  else
		  !scan for en_rxn from initial nuclei array  	  	  
		  end if
		  endif
	  
		  if (count(enn_rxn) /= count(nne_rxn)) then
		  if (count(enn_rxn) > count(nne_rxn)) then
		  !scan back for 2ne_rxn from initial nuclei array
		  else
		  !scan back for e2n_rxn from initial nuclei array 
		  endif
		  endif
		  ! end of check for baryon conservation 

! Z and A after reactions stored in Z_fin, A_fin arrays
!		  Z_fin(k) = Z
!		  A_fin(k) = A		 
!		end do  ! end of nuclei loop 
  
  		  ! output chain to file 	
		  do l= 1, size(Z_fin)
		  write(final_id,'(3(e10.5,2x),4(I10,2x))') pressure, mu_e, mu_n, Z_int(l), A_int(l),&
				& Z_fin(l), A_fin(l)
		  enddo 
		  close(final_id)
		  stop
	
	   end do  ! end of mu loop
	  end do	! end of pressure loop 
   
		contains	
	
	  subroutine get_nucleus_properties(Z,A,id,BE,Sn,S2n,Sp,S2p,ecthresh,bthresh,VN,ierr)
		integer, intent(in) :: Z,A
		integer, intent(out) :: id
		real, intent(out) :: BE,Sn,S2n,Sp,S2p,ecthresh,bthresh,VN
		integer, intent(out) :: ierr
		integer :: N
		
		ierr = 0
		N = A-Z
		id = mass_table_index(Z,N,ierr)
		if (ierr /= 0) return
	
		! set the properties
		BE = mt% BE(id)
		Sn = mt% Sn(id)
		S2n = mt% S2n(id)
		Sp = mt% Sp(id)
		S2p = mt% S2p(id)
		ecthresh = mt% Ec(id)
		bthresh = mt% beta(id)
		VN = mt% VN(id)	
	  end subroutine get_nucleus_properties
	
	  subroutine print_reaction()
		character(len=*), parameter :: form = '(2(e12.5),4i4)'
		write (fid,form) mu_n,mu_e,Z,A,Zr,Ar
	  end subroutine print_reaction
	
	  subroutine print_reaction_check()
		character(len=*), parameter :: form = '((a6),2(i6))'
		write(fid,form) rxn, Z, A
	  end subroutine print_reaction_check
	
	  function io_failure(ierr,message)
      	integer, intent(in) :: ierr
      	character(len=*), intent(in) :: message
      	logical :: io_failure
      	if (ierr == 0) then
         io_failure = .FALSE.
         return
      	end if
      	write (error_unit,'(a,i0)') 'ERROR: '//trim(message)//'; ierr = ',ierr
      	io_failure = .TRUE.
      end function io_failure
   
      function neutron_k(x)
      	real, intent(in) :: x
      	real :: neutron_k
      	neutron_k = neutron_pressure(x) - pres_n
      end function neutron_k
   
      function neutron_k_negative(x)
      	real, intent(in) :: x
      	real :: neutron_k_negative 
      	neutron_k_negative = abs(neutron_pressure(x)) - abs(pres_n)
      end function neutron_k_negative

end program follow_chain
