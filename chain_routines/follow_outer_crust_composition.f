module outer_crust

  contains 

  subroutine follow_outer_crust_composition
	  use iso_fortran_env, only : error_unit, output_unit
	  use phys_constants
	  use utils_lib
	  use mb77, only : neutron_chemical_potential, electron_pressure, neutron_pressure
	  use mass_table
	  use pressure_table
	  use alert_lib
	  use rootfind	
	  use ash_table

	  character(len=*), parameter :: datadir = '../../../data'	
	  character(len=*), parameter :: default_dir = 'crust_eos'
	  character(len=*), parameter :: default_infile = 'moe_chain.inlist'
	  character(len=*), parameter :: default_dist_file = 'nuclei_distribution.data'
	  character(len=*), parameter :: default_mass_table = 'moe95_converted.data'
	  character(len=*), parameter :: default_stable_table = 'stable_nuclei_moe95.data'
	  character(len=*), parameter :: default_pressure_table = 'pressure_moe95.data'
	  character(len=*), parameter :: ash_table_name = 'ash_initial.data'
	  character(len=*), parameter :: final_file = 'ash_final.data'
	  real :: kn, ke, mu_n
	  real :: mu_e, ne, nn 
	  real :: x1, x2, xacc
	  real :: n, pres_n, pressure
	  real :: pressure_start, pressure_stop, pressure_increment
	  real :: mu_e_start, mu_e_stop, mu_e_increment
	  real :: Sn, S2n, Sp, S2p, B, ecthresh, bthresh, VN
	  real :: Snr, S2nr, Spr, S2pr, Br, ecthreshr, bthreshr, VNr
	  real :: alpha(2), beta(2), gamma(2), delta(2), epsilon(2)
	  real, dimension(:), pointer :: abun_initial, abun_final
	  integer, dimension(:), pointer :: Z_initial, A_initial
	  integer, dimension(:), pointer :: Z_final, A_final, N_final
	  integer :: i, j, i_enter
	  integer :: Z, A, Zr, Ar, Nr, inlist_id
	  integer :: dist_id, index
	  integer :: ierr, id, ios, iter, iZ, iZb, iZe, iEq(1)
	  integer :: k, l, final_id
	  integer, parameter :: ineg = 1, ipos = 2
	  integer, parameter :: fid = output_unit 
	  integer, parameter :: max_iterations = 100
	  real, parameter :: del_m = mn_n-mp_n-me_n
	  character(len=256) :: arg, mass_table_used, pfile, mass_table_name, filename
	  character(len=256) :: stable_table_used, pressure_table_used, dist_file
	  character(len=6) :: rxn
	  real,dimension(:),pointer :: hist
	  type(mass_table_type), pointer :: mt
	  type(pressure_table_type), pointer :: pt
	  type(ash_table_type), pointer :: at
	  logical, dimension(max_iterations) :: neutron_capture, dineutron_capture
	  logical, dimension(max_iterations) :: neutron_emission, dineutron_emission
	  logical, dimension(max_iterations) :: en_rxn, enn_rxn, ne_rxn, nne_rxn
	  logical, save :: stable_table_loaded = .FALSE.
	  logical, save :: mass_table_loaded = .FALSE.
	  logical, save :: pressure_set = .FALSE.
	  logical, save :: pressure_table_loaded = .FALSE.
	  logical, save :: ash_table_loaded = .FALSE.
	
	  namelist /io/ pressure_table_used
	  namelist /range/ mu_e_start, mu_e_stop, &
				pressure_start, pressure_stop
				
      ! set defaults
      pfile = default_infile
      dist_file = default_dist_file
      mass_table_used = default_mass_table
   	  stable_table_used = default_stable_table
   	  pressure_table_used = default_pressure_table
      mu_e_start = 0.0 !MeV
      mu_e_stop = 30.0 !MeV
      mu_n = 0.0 !MeV
   	
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
   	  if (mass_table_loaded .eqv. .FALSE.) then
	  call load_mass_table('../../../data',trim(mass_table_used),ierr)
	  mass_table_loaded = .TRUE. 
	  if (ierr /= 0) then
		write (error_unit,'(a)') 'mass table loading error'
		stop
	  end if
	  end if
	  mt => winvn_mass_table
	  write(*,*) 'Mass table loaded...'

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
	  write(*,*) 'Pressure table loaded...'
	
	  !load ash table
	  if (ash_table_loaded .eqv. .FALSE.) then
	  call load_ash_table('../../../data', trim(ash_table_name), ierr)
	  ash_table_loaded = .TRUE.
	  if (ierr /= 0) then
	  write(error_unit, '(a)') trim(alert_message)
	  stop
	  end if
	  end if
	  at => winvn_ash_table
	  write(*,*) 'Ash table loaded...'	

	  allocate(Z_initial(at% Ntable), A_initial(at% Ntable), abun_initial(at% Ntable), &
	  &		  Z_final(at% Ntable), N_final(at% Ntable), A_final(at% Ntable), abun_final(at% Ntable))
	
	  ! open output file for composition at end of outer crust
      final_id = alloc_iounit(ierr)
      if (io_failure(ierr,'allocating unit for final array file')) stop
      filename = trim(datadir)//'/'//default_dir//'/'//trim(final_file)
      open(unit=final_id, file=filename, iostat=ios, status='unknown') 
      if (io_failure(ios,'opening final array file')) stop
      write(final_id,'(a10)') 'Pressure'
      write(final_id, '(5(a10,2x))') 'Z', 'N', 'A', 'BE', 'Y' 
   
      !main loop over pressure (pushes nucleus to higher pressures)    
	  
	  ! store values from ash table into arrays
   	  Z_initial = at% Z
   	  A_initial = at% A
   
      ! scan through pressure values from pressure table
      do i = 1, pt% Ntable
	  pressure = pt% P(i)	! MeV fm**-3
	  if (pressure > 3.615d-4) exit
	  mu_n = pt% mu_n(i)
 	  mu_e = pt% mu_e(i)
 	  write(*,*) '--------'
	
	  ! end of last cascade is moved up in pressure 
	  if (i > 1 .and. Z_final(1) > 0) then
	  Z_initial = Z_final
	  A_initial = A_final
	  end if
 		!initial nuclei accreted to the given pressure 
	 	do k = 1, at% Ntable
		Z = Z_initial(k)
	    A = A_initial(k)
	 
	  ! set defaults
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
      Nr = Ar-Zr
!      index= Zr - mt%Zmin + 1
!      if (Nr < mt% Nmin(index) .or. Nr > mt% Nmax(index)) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr, &
      		&	ecthreshr,bthreshr,VNr,ierr)
      if (ierr /= 0) then
         write(error_unit,'((a,2x),2(i4,2x))') 'unable to find nucleus', Zr, Ar
!         stop
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
      Nr = Ar-Zr
!      index= Zr - mt%Zmin + 1
!      if (Nr < mt% Nmin(index) .or. Nr > mt% Nmax(index)) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
      if (ierr /= 0) then
         write(error_unit,'((a,2x),2(i4,2x))') 'unable to find nucleus', Zr, Ar
!         stop
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
      Nr = Ar-Zr
!      index= Zr - mt%Zmin + 1
!      if (Nr < mt% Nmin(index) .or. Nr > mt% Nmax(index)) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit    
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
      if (ierr /= 0) then
         write(error_unit,'((a,2x),2(i4,2x))') 'unable to find nucleus', Zr, Ar
!         stop
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
      Nr = Ar-Zr
!      index= Zr - mt%Zmin + 1
!      if (Nr < mt% Nmin(index) .or. Nr > mt% Nmax(index)) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
      if (ierr /= 0) then
         write(error_unit,'((a,2x),2(i4,2x))') 'unable to find nucleus', Zr, Ar
!         stop
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
      
      cycle
      
      ! check for strong reactions     
      ! neutron capture
      Ar = A+1; Zr = Z
      Nr = Ar-Zr
!      index= Zr - mt%Zmin + 1
!      if (Nr < mt% Nmin(index) .or. Nr > mt% Nmax(index)) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit  
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
      if (ierr /= 0) then
         write(error_unit,'((a,2x),2(i4,2x))') 'unable to find nucleus', Zr, Ar
!         stop
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
      Nr = Ar-Zr
!      index= Zr - mt%Zmin + 1
!      if (Nr < mt% Nmin(index) .or. Nr > mt% Nmax(index)) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit     
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
      if (ierr /= 0) then
         write(error_unit,'((a,2x),2(i4,2x))') 'unable to find nucleus', Zr, Ar
!         stop
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
      Nr = Ar-Zr
!      index= Zr - mt%Zmin + 1
!      if (Nr < mt% Nmin(index) .or. Nr > mt% Nmax(index)) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit      
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
      if (ierr /= 0) then
         write(error_unit,'((a,2x),2(i4,2x))') 'unable to find nucleus', Zr, Ar
!         stop
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
      Nr = Ar-Zr
!      index= Zr - mt%Zmin + 1
!      if (Nr < mt% Nmin(index) .or. Nr > mt% Nmax(index)) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit     
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
      if (ierr /= 0) then
         write(error_unit,'((a,2x),2(i4,2x))') 'unable to find nucleus', Zr, Ar
!         stop
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
      Nr = Ar-Zr
!      index= Zr - mt%Zmin + 1
!      if (Nr < mt% Nmin(index) .or. Nr > mt% Nmax(index)) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit 
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
      if (ierr /= 0) then
         write(error_unit,'((a,2x),2(i4,2x))') 'unable to find nucleus', Zr, Ar
 !        stop
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
      Nr = Ar-Zr
!      index= Zr - mt%Zmin + 1
!      if (Nr < mt% Nmin(index) .or. Nr > mt% Nmax(index)) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit   
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,Vnr,ierr)
      if (ierr /= 0) then
      	 write(error_unit, '((a,2x),2(i4,2x))') 'unable to find nucleus', Zr, Ar
!      	 stop
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
      
      end do  ! end of iteration loop      
  	 
  	 ! add Z and A to final array after chain 
  	  Z_final(k) = Z
	  A_final(k) = A      
	  N_final(k) = A-Z   
      end do  ! end of nuclei loop 
      
      ! print final ash file
      write(final_id,'(e10.5)') pressure !add final pressure reached by outer crust solver
  	  do l= 1, at% Ntable
	  call get_nucleus_properties(Z_final(l),A_final(l),id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,Vnr,ierr)
	  write(final_id,'(3(I10,2x),2(e10.5,2x))') Z_final(l), N_final(l), A_final(l), Br, at% Y(l)
	  enddo 
	  close(final_id)  	  
      end do	! end of pressure loop 
      
      call ash_table_shutdown
      write(*,*) 'Finished with outer crust...'
   
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
	 end subroutine get_nucleus_properties
	
	 subroutine print_reaction()
		character(len=*), parameter :: form = '(3(e12.5,2x),2(i4,2x),(a6,2x),2(i4,2x))'
		write (fid,form) pressure,mu_n,mu_e,Z,A,rxn,Zr,Ar
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

	end subroutine follow_outer_crust_composition

end module outer_crust