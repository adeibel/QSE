program follow_chain
	use iso_fortran_env, only : error_unit, output_unit
	use phys_constants
	use utils_lib
	use mb77, only : neutron_chemical_potential, electron_pressure, neutron_pressure
	use mass_table
	use pressure_table
	use alert_lib
	use rootfind	
	
	character(len=*), parameter :: default_infile = 'moe_chain.inlist'
	character(len=*), parameter :: default_dist_file = 'nuclei_distribution.data'
	character(len=*), parameter :: default_mass_table = 'moe95_converted.data'
	character(len=*), parameter :: default_stable_table = 'stable_nuclei_moe95.data'
	character(len=*), parameter :: default_pressure_table = 'pressure_moe95.data'
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
	integer :: i, j, i_enter
	integer :: Z, A, Zr, Ar, inlist_id
	integer :: dist_id
	integer :: ierr, id, ios, iter, iZ, iZb, iZe, iEq(1)
	integer :: Z_int(1), A_int(1)
	integer :: Z_fin(1), A_fin(1)
	!integer :: Z_int(5549), A_int(5549)
	!integer :: Z_fin(5549), A_fin(5549)
	integer :: k, l, final_id
	integer, parameter :: ineg = 1, ipos = 2
	integer, parameter :: fid = output_unit, max_iterations = 100
	real, parameter :: del_m = mn_n-mp_n-me_n
	real, parameter :: mun_max = 6.5
	character(len=256) :: arg, mass_table_used, pfile, mass_table_name
	character(len=256) :: stable_table_used, pressure_table_used, dist_file
	character(len=6) :: rxn
	real,dimension(:),pointer :: hist
	type(mass_table_type), pointer :: mt
	type(pressure_table_type), pointer :: pt
	logical, dimension(max_iterations) :: neutron_capture, dineutron_capture
	logical, dimension(max_iterations) :: neutron_emission, dineutron_emission
	logical, dimension(max_iterations) :: en_rxn, enn_rxn, ne_rxn, nne_rxn
	logical, save :: stable_table_loaded = .FALSE.
	logical, save :: mass_tables_are_loaded = .FALSE.
	logical, save :: pressure_set = .FALSE.
	logical, save :: pressure_table_loaded = .FALSE.
	logical, save :: dist_file_loaded = .FALSE.
	logical :: outer_crust
	
	namelist /io/ pressure_table_used
	namelist /range/ Z, A, mu_e_start, mu_e_stop, &
				pressure_start, pressure_stop, outer_crust
				
   ! set defaults
   pfile = default_infile
   dist_file = default_dist_file
   mass_table_used = default_mass_table
   stable_table_used = default_stable_table
   pressure_table_used = default_pressure_table
   Z = 26
   A = 56
   mu_e_start = 0.0 !MeV
   mu_e_stop = 30.0 !MeV
   mu_n = 0.0 !MeV
   
   !read in array of initial nuclei distribution 
   ierr = 0
   dist_id = alloc_iounit(ierr)
   if (io_failure(ierr,'allocating unit for distribution file')) stop
   open(unit=dist_id, file=dist_file, iostat=ios, status="old", action="read")
   ! loop over file 
   read(dist_id,*) Z_int(i), A_int(i), abun(i)
   !
   dist_file_loaded = .TRUE.
   close(dist_id)
   call free_iounit(dist_id)
	
   ! read in the inputs
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

	!load pressure table
	if (pressure_table_loaded .eqv. .false.) then
	call load_pressure_table('../../../data',trim(pressure_table_used), ierr)
	pressure_table_loaded = .true.
	if (ierr /= 0) then
		write(error_unit,'(a)') 'pressure table loading error'
		stop
	end if
	end if
	pt => winvn_pressure_table
	
   final_id = alloc_iounit(ierr)
   if (io_failure(ierr,'allocating unit for final array file')) stop
   open(unit=final_id, file=final_file, iostat=ios, status='unknown') 
   if (io_failure(ios,'opening final array file')) stop
   write(final_id,'(3(a10,2x),4(A10,2x))') 'Pressure', 'mu_e', 'mu_n', 'Zint', 'Aint', &
   			& 'Zfin', 'Afin'
   
   !main loop over pressure (pushes nucleus to higher pressures)   
   
    ! make table of nuclei the first initial nuclei array 
    do i = 1, size(Z_int)
   	Z_int(i) =  Z !mt% Z
	A_int(i) =  A !mt% A
    end do
   
   do i = 1, pt% Ntable
	pressure = pt% P(i)	! MeV fm**-3
	mu_n = pt% mu_n(i)
	mu_e = pt% mu_e(i)
	write(*,*) pressure, mu_e, mu_n
	write(*,*) pt% P(i), pt% mu_e(i), pt% mu_n(i)
	
	! end of last cascade is moved up in pressure 
	if (i >= 1 .and. Z_fin(1) > 0) then
	Z_int = Z_fin
	A_int = A_fin
	end if
 		
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
	 
      ! loop over rxns 
      do iter = 1, max_iterations
   	  !get properties of nucleus that is being pushed deeper
	  call get_nucleus_properties(Z,A,id,B,Sn,S2n,Sp,S2p,ecthresh,bthresh,VN,ierr)
	 
      ! check for weak reactions                 
      ! electron capture
      Ar = A; Zr = Z-1
      if (Ar-Zr < mt% Nmin .or. Ar-Zr > mt% Nmax) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit
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
      if (Ar-Zr < mt% Nmin .or. Ar-Zr > mt% Nmax) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit
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
      if (Ar-Zr < mt% Nmin .or. Ar-Zr > mt% Nmax) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit      
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
      if (Ar-Zr < mt% Nmin .or. Ar-Zr > mt% Nmax) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit      
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
      if (ierr /= 0) then
         write(error_unit,'(a)') 'unable to find nucleus'
         stop
         exit
      end if
      beta(ipos) = 2.0*mu_n - mu_e + del_m + (B - Br)
      if (beta(ipos) < 0) then
      	 enn_rxn(iter) = .TRUE.
         rxn = '(e,2n)'
         call print_reaction
         A = Ar; Z = Zr
         cycle
      end if
      
      if (outer_crust .eqv. .true.) cycle

      ! check for strong reactions     
      ! neutron capture
      Ar = A+1; Zr = Z
      if (Ar-Zr < mt% Nmin .or. Ar-Zr > mt% Nmax) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit      
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
      if (Ar-Zr < mt% Nmin .or. Ar-Zr > mt% Nmax) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit      
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
      if (Ar-Zr < mt% Nmin .or. Ar-Zr > mt% Nmax) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit      
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
      if (Ar-Zr < mt% Nmin .or. Ar-Zr > mt% Nmax) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit      
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
      if (Ar-Zr < mt% Nmin .or. Ar-Zr > mt% Nmax) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit      
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
      if (Ar-Zr < mt% Nmin .or. Ar-Zr > mt% Nmax) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit      
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

	  Z_fin(k) = Z
	  A_fin(k) = A
	                     
	goto 1

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
         ! write(error_unit,'(a)') 'unable to find nucleus'
         ! stop
         ! exit
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

1 continue 

      end do  ! end of iteration loop
        
  	  Z_fin(k) = Z
	  A_fin(k) = A
         
    end do  ! end of nuclei loop 
  
  	  do l= 1, size(Z_fin)
	  write(final_id,'(3(e10.5,2x),4(I10,2x))') pressure, mu_e, mu_n, Z_int(l), A_int(l),&
	  		& Z_fin(l), A_fin(l)
	  enddo 
	  close(final_id)
	    
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
