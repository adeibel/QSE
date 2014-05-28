! finds the stable nuclei at a given pressure for as many mu_e, mu_n
! combinations as desired 

program pressure_scan
   use iso_fortran_env, only : error_unit
   use utils_lib
   use alert_lib
   use phys_constants
   use check_stability
   use mb77
   use rootfind

   character(len=*), parameter :: default_pfile = 'pressure_scan.inlist'
   character(len=*), parameter :: equil_file = 'equil_nuclei_moe95.data'
   character(len=*), parameter :: stable_file = 'stable_nuclei_moe95.data'
   character(len=*), parameter :: pressure_file = 'pressure_moe95.data'
   character(len=*), parameter :: mu_file = 'mu_max.data'
   character(len=80) :: pfile, mass_table_name
   real :: mu_e_start, pressure_start
   real :: mu_e_stop, pressure_stop
   real :: mu_e_increment, pressure_increment
   real :: mu_e, mu_n, pressure, pres_n
   real :: ke, kn, n
   real :: kT
   real :: x1, x2, xacc
   real :: ke_limit, kn_limit, ne_limit, nn_limit
   real,dimension(:),pointer :: hist
   integer :: i, j, ios, ierr
   integer :: inlist_id, equil_id, stable_id, mu_id, pressure_id
   integer :: Z_best(1), A_best(1)
   logical :: inner_crust, outer_crust, no_stable
   logical, save :: ke_solved = .FALSE.
   logical, save :: kn_solved = .FALSE.

   namelist /io/ mass_table_name
   namelist /range/ mu_e_start, pressure_start, &
				mu_e_stop, pressure_stop, kT, &
				inner_crust, outer_crust
     
   ! set the name of the inlist file name  
   if (command_argument_count() == 1) then
      call get_command_argument(1, pfile)
   else
      pfile = default_pfile
   end if
     
   ! set defaults in MeV
   mass_table_name = 'moe95_converted.data'
   !mass_table_name = 'mass_mb77_1_1.data'
   mu_e_start = 30.0
   pressure_start = 1.5d-5*hbarc_n
   mu_e_stop = 100.0
   pressure_stop = 1.2d-1*hbarc_n 
   kT = 1.0d-2  
     
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
    
   if (inner_crust .eqv. outer_crust) then
   write(*,*) 'Must choose between inner and outer crust in inlist'
   stop
   end if
   
   mu_id = alloc_iounit(ierr)
   if (io_failure(ierr,'allocating unit for mu max file')) stop
   open(unit=mu_id, file=mu_file, iostat=ios, status='unknown') 
   if (io_failure(ios,'opening mu max file')) stop
   write(mu_id,'(3(a10,2x))') 'Pressure', 'mu_e_max', 'mu_n_max'
    
   equil_id = alloc_iounit(ierr)
   if (io_failure(ierr,'allocating unit for equilibrium nuclei file')) stop
   open(unit=equil_id, file=equil_file, iostat=ios, status="unknown")
   if (io_failure(ios,'opening equilibrium nuclei file')) stop   
   write(equil_id,'(5(a10,2x))') 'Pressure', 'mu_e', 'mu_n', 'Z', 'A' 
   
   stable_id = alloc_iounit(ierr)
   if (io_failure(ierr,'allocating unit for stable nuclei file')) stop
   open(unit=stable_id, file=stable_file, iostat=ios, status="unknown")
   if (io_failure(ios,'opening stable nuclei file')) stop 
   write(stable_id,'(3A13,3A5,4A12,2x,A12,2x,A12)') 'P[MeV fm**-3]', 'mu_e [MeV]', 'mu_n [MeV]',&
   		'Z', 'A', 'N', 'BE [MeV]', 'BE/A [MeV]', 'Sn [MeV]', 'Sp [MeV]', 'eden[MeV**4]', &
   		'rho [MeV**4]' 

	pressure_id = alloc_iounit(ierr)
	if (io_failure(ierr,'allocating unit for pressure file')) stop
	open(unit=pressure_id, file=pressure_file, iostat=ios, status="unknown")
	if (io_failure(ios,'opening pressure file')) stop
	write(pressure_id, '(3(A13,2x))') 'P[MeV fm**-3]', 'mu_e [MeV]', 'mu_n [MeV]'

	pressure = pressure_start

do j= 1, 500
	pressure = pressure*2. !1.d-1     !*hbarc_n*real(j)

	no_stable = .FALSE.

	 do i=0,100
	   
	   ! find P(ke) that can support the total pressure
	   if (ke_solved .eqv. .FALSE.) then 	 
	   x1=0.0
	   x2=10.0
	   xacc=1.d-20
	   ke_limit=root_bisection(p_e,x1,x2,xacc,ierr,hist) !returns in fm**-1	 
	   if (ierr /= 0) then
       write(*,*) 'Error in bisection for ke_limit wave vector'
	   stop
       endif  
       ke_solved = .TRUE.
       
       ne_limit=ke_limit**3/threepisquare 
       
       ! find P(kn) that can support the total pressure
       x1=0.0
       x2=10.0
       xacc=1.d-20
       kn_limit=root_bisection(p_n,x1,x2,xacc,ierr,hist) !returns in fm**-1
       if (ierr /= 0) then
       write(*,*) 'Error in bisection for kn_limit wave vector'
       stop
       endif
       
       nn_limit=2.0*kn_limit**3/threepisquare 
       
       write(mu_id,'(3(es12.5,2x))') pressure, electron_chemical_potential(ne_limit), &
       			& neutron_chemical_potential(nn_limit)
       			
 
       endif
	   
	   ke = ke_limit-ke_limit/100.*real(i)		  ! fm**-1
	   mu_e = sqrt((ke*hbarc_n)**2 + me_n**2)     ! MeV

	   !mu_e = mu_e_start+mu_e_increment*real(i)  !MeV
	   !ke = sqrt(mu_e**2-me_n**2)/hbarc_n        !put in fm**-1
	   	  
       pres_n = pressure - electron_pressure(ke)
       
      if (pres_n == 0.) then
      mu_n = 0.0
      endif
       
      if (pres_n < 0.) cycle
       
      if (pres_n > 0.) then
      x1=0.0
      x2=10.0
      xacc=1.d-20
      !need neutron pressure for this next rootfind
      kn = root_bisection(neutron_k, x1, x2, xacc, ierr, hist) !returns in fm**-1
       if (ierr /= 0) then
       write(*,*) 'Error in bisection for kn wave vector'
	   stop
       endif
  	  n=2.0*kn**3/threepisquare 
  	  mu_n = neutron_chemical_potential(n) !returns in MeV
  	  end if 
  	  
  	  if (outer_crust .eqv. .TRUE.) then
  	  mu_n = 0.0
  	  end if
  	  
   !  call stable_nuclei(mu_e, mu_n, kT, mass_table_name, Z_best(1), &
   !  		A_best(1), pressure, stable_id, mu_id )
	
    !  	 if (no_stable .eqv. .TRUE.) exit	
      				
    !  write(equil_id,'(3(es12.5,2x),2(I6,2x))') pressure, mu_e, mu_n, Z_best(1), A_best(1)
    ! enddo 
    
    if (outer_crust .eqv. .true.) then
    pressure = electron_pressure(ke_limit)
    mu_e = electron_chemical_potential(ne_limit)
    mu_n = 0.
    end if
    
    write(pressure_id,'(3(es13.5,2x))') pressure, mu_e, mu_n
    
    if(outer_crust .eqv. .true.) then
    exit
    end if
    
    enddo
    
     ke_solved = .FALSE.
    
    enddo 
    
       close(equil_id)
       call free_iounit(equil_id)
       close(stable_id)
       call free_iounit(stable_id)
       close(mu_id)
       call free_iounit(mu_id)
       close(pressure_id)
       call free_iounit(pressure_id)
    
 contains
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
   
   function p_e(x)
      real, intent(in) :: x
      real :: p_e     
      p_e = electron_pressure(x) - pressure     
    end function p_e  
    
   function p_n(x)
      real, intent(in) :: x
      real :: p_n    
      p_n = neutron_pressure(x) - pressure     
    end function p_n
    
   
end program pressure_scan
