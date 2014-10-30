      module inner_crust
 
      use num_def
      use num_lib
      use iso_fortran_env, only: error_unit, oid=>output_unit
      use alert_lib
      use utils_lib
      use mass_table
      use phys_constants
      use eos_table
      use dist_table
      use rootfind      
      use crust_eos_lib 
      use eos_lib     

      implicit none
      
      logical, parameter :: dbg = .false.
      
      ! dimensions
      integer, parameter :: nz = 1 ! number of zones
      integer :: nvar ! number of variables per zone
	  integer :: neq  ! number of equations (nvar*nz)

      ! information about the bandwidth of the jacobian matrix
      ! we have a square matrix, so set number zones sub and super both to 0
      integer, parameter :: stencil_zones_subdiagonal = 0
      integer, parameter :: stencil_zones_superdiagonal = 0
      integer :: m1 ! number of subdiagonals
      integer :: m2 ! number of superdiagonals 

      ! equation residuals
      real*8, pointer, dimension(:,:) :: equ
      ! equ(i) has the residual for equation i, i.e., the difference between
      ! the left and right hand sides of the equation.  So the goal
      ! of the solver is to find an approximate solution that makes the
      ! magnitude of the residuals small.

      ! the primary variables
      real*8, pointer, dimension(:,:) :: x ! new vector of primaries
      real*8, pointer, dimension(:,:) :: xold ! old vector of primaries
      real*8, pointer, dimension(:,:) :: dx ! increment vector -- on input is initial guess for x - xold
      real*8, pointer, dimension(:,:) :: xscale ! typical values

      ! the secondary variables
      integer, parameter :: nsec = 0 ! number of secondaries per zone
      real*8, pointer :: t1, t2
      integer, parameter :: ldy = nz+2 ! leading dimension of y, >= nz
      real*8, pointer, dimension(:,:) :: y ! the values

      logical :: first_step, nonconv, numerical_jacobian, doing_jacobian
      real*8 :: tol_correction_norm, tol_max_correction, tol_residual_norm, epsder, &
         dt1_dx0, dt1_dx1, dt1_dx2, dt2_dx0, dt2_dx1, dt2_dx2
      integer :: matrix_type, matrix_type_current

      integer, parameter :: max_lid = 50000, max_lrd = 100000
      integer, target :: ipar_decsol(max_lid)
      real*8, target :: rpar_decsol(max_lrd)

	  type qse_table_type
		integer :: Ntable	! number of table entries
		integer, dimension(:), pointer :: Z, N, A 
		real, dimension(:), pointer :: BE !binding energy 
		real, dimension(:), pointer :: Y !abundance fraction 
	  end type qse_table_type      
	  type(qse_table_type), target :: winvn_qse_table  
      
      type(mass_table_type), pointer :: mt
      type(eos_table_type), pointer :: et
      type(dist_table_type), pointer :: dt
      type(qse_table_type), pointer :: qt
      real, pointer, dimension(:) :: mu_i, ni
      real,dimension(:),pointer :: hist      
      real*8 :: mu_e, mu_n, Y_e, Y_n 
   	  real*8 :: n_b, n_e, n_n
   	  real*8 :: p_ext
      real*8 :: ke, kn
      real :: lambda_baryon, lambda_charge, lambda_pressure
      real :: Z_bar, A_bar
      real :: x1, x2, xacc
 	  real, save :: zsum_save = 0.
 	  real :: asum, zsum
      real :: kT
      real :: pres_n
      real :: Pn
      real :: n_b_prev
      integer :: i, j  
                          
      contains
      
      subroutine follow_inner_crust_composition
         use mtx_lib
         use mtx_def
     
      integer :: ierr, liwork, lwork, lid, lrd, which_decsol
      integer :: which_decsol_in, decsol     
      integer, dimension(:), pointer :: iwork
      integer, parameter :: lrpar = 0, lipar = 0
      integer, target :: ipar(lipar)
      integer :: i, j, ios
      integer :: inlist_id, output_id 
      integer :: mu_table_input_id, mu_table_output_id
      integer :: y_output_id, ash_id, tov_id      
      real*8, dimension(:), pointer :: work 
      real*8, target :: rpar(lrpar)   
      real :: p_ext_start, n_b_start
      real :: Y_sum
      real :: eps_const, eps_sol
      real :: rho
      real :: mterm, m_nuc, m_star 
      real :: Z_average, A_average, BE_average
      real, dimension(:), pointer :: fac1, fac2, mass_frac
      real, parameter :: g = 1.d0      
      character(len=64) :: decsol_option_name 	  
 	  character(len=*), parameter :: dist_table_name = 'ash_final.data'
	  character(len=*), parameter :: mass_table_name = 'moe95_converted.data'
 	  character(len=*), parameter :: eos_table_name = 'ashes_acc.txt'
	  character(len=*), parameter :: y_output_file = 'y_output.data'
	  character(len=*), parameter :: tov_output_file = 'tov_acc.data'
      character(len=*), parameter :: output_file = 'qse_output.data'
   	  character(len=*), parameter :: default_infile = 'qse.inlist'
   	  character(len=64), parameter :: data_dir = '../../../data/crust_eos'
	  character(len=64), parameter :: mu_table_name = 'mu_table_old.data'
      character(len=80) :: infile	
      logical :: do_numerical_jacobian
      logical, save :: eos_table_is_loaded = .FALSE.


      namelist /range/ P_ext_start, n_b_start, kT, &
      	&	do_numerical_jacobian, which_decsol_in
      	
      qt => winvn_qse_table
     
      ! set the name of the inlist file name  
      if (command_argument_count() == 1) then
       call get_command_argument(1, infile)
      else
       infile = default_infile
      end if
      
      ! set defaults 
      P_ext_start = 5.d-10!fm^-4
      n_b_start = 5.0d-8 !fm^-3
      kT = 1.0d-2  !MeV
      do_numerical_jacobian = .true.
      decsol = lapack
    
      ! read in the inputs
 	  ierr = 0
  	  inlist_id = alloc_iounit(ierr)
 	  if (io_failure(ierr,'allocating unit for namelist')) stop
   
  	  open(unit=inlist_id, file=infile, iostat=ios, status="old", action="read")
  	  if (io_failure(ios,'opening inlist file'//trim(infile))) stop   
 	  read(inlist_id,nml=range,iostat=ios)
	  if (io_failure(ios,'unable to read namelist "range"')) stop
	  close(inlist_id)
	  call free_iounit(inlist_id)
 	  decsol = which_decsol_in 
 
 	  ! load mass table 
      if (mass_table_is_loaded .eqv. .FALSE.) then
      call load_mass_table('../../../data',trim(mass_table_name),ierr)
      mass_table_is_loaded = .TRUE.
      if (ierr /= 0) then
      write(error_unit,'(a)') trim(alert_message)
      stop
      end if
      end if
      mt => winvn_mass_table    
 	  write(*,*) 'Mass table loaded...'
 	  
 	  ! load eos table
 	  if (eos_table_is_loaded .eqv. .FALSE.) then
 	  call load_eos_table('../../../data',trim(eos_table_name), ierr)
 	  eos_table_is_loaded = .TRUE.
 	  if (ierr /= 0) then
 	  write(error_unit, '(a)') trim(alert_message)
 	  stop
 	  end if
 	  end if
 	  et => winvn_eos_table
 	  write(*,*) 'EOS table loaded...'

	  !load dist table
	  if (dist_table_is_loaded .eqv. .FALSE.) then
	  call load_dist_table('../../../data', trim(dist_table_name), ierr)
	  dist_table_is_loaded = .TRUE.
	  if (ierr /= 0) then
	  write(error_unit, '(a)') trim(alert_message)
	  stop
	  end if
	  end if
	  dt => winvn_dist_table
	  write(*,*) 'Dist table loaded...'
 	  
 	  ! allocate units and open output files   
 	  output_id = alloc_iounit(ierr)
  	  if (io_failure(ierr,'allocating unit for output nuclei file')) stop
  	  open(unit=output_id, file=output_file, iostat=ios, status="unknown")
   	  if (io_failure(ios,'opening output nuclei file')) stop 
	  write(output_id,'(7A13)') 'n_b [fm^-3]', 'k_e [fm^-1]', &
   			& 'k_n [fm^-1]', 'mu_e [MeV]', 'mu_n [MeV]', 'Z_bar', 'A_bar'

	  y_output_id = alloc_iounit(ierr)
	  if (io_failure(ierr, 'allocating unit for y fractions file')) stop
	  open(unit=y_output_id, file = y_output_file, iostat=ios, status="unknown")
	  if (io_failure(ios,'opening y_output file')) stop
	  write(y_output_id,'(8(A12),2x)') 'Pressure', 'n_b', 'Ye', 'Yn', 'Zbar', 'Abar', 'mu_e', 'mu_n'
	  
	  tov_id = alloc_iounit(ierr)
	  if (io_failure(ierr, 'allocating unit for tov file')) stop
	  open(unit=tov_id, file=tov_output_file, iostat=ios, status="unknown")
	  if (io_failure(ios,'opening tov file')) stop	  
	  write(tov_id,'(5(A15),2x)') 'P[MeV*fm^-3]', 'rho[MeV*fm^-3]', 'eps[MeV^-4]', 'mue[MeV]', 'mun[MeV]' 
      
  	  ! solve for qse distribution at each pressure	  
  	  do i=1, et% Ntable
  	  	 mu_n = (et% mun(i))*hbarc_n - mn_n
  	     mu_e = (et% mue(i))*hbarc_n
  	     p_ext = (et% pr(i))*hbarc_n !p_ext_start*hbarc_n*real(i)  
		        
         which_decsol = which_decsol_in
         call decsol_option_str(which_decsol, decsol_option_name, ierr)
         if (ierr /= 0) return

         if (which_decsol == mkl_pardiso) then
            if (.not. okay_to_use_mkl_pardiso()) which_decsol = lapack
         end if
         
         numerical_jacobian = do_numerical_jacobian
		 		 
		! electron wave vector fm^-1
		if (mu_e > 0.) then
        x1=0.0
        x2=10.
        xacc=1.d-15
        ke=root_bisection(ke_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
        if (io_failure(ierr,'Error in bisection for ke wave vector')) stop
        n_e = ke**3/threepisquare               
        Y_e = n_e/n_b   
        else
        mu_e = abs(mu_e)
        x1=0.0
        x2=10.
        xacc=1.d-15
        ke=root_bisection(ke_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
        if (io_failure(ierr,'Error in bisection for ke wave vector')) stop
        n_e = ke**3/threepisquare               
        Y_e = n_e/n_b   
        mu_e = -mu_e
        end if		 
		 
		n_b = n_e/(et% Ye(i))
		         
         if (p_ext < dt% P) cycle
        ! if (i<1500) cycle 
                
         write(*,*) 'Pressure =', p_ext
         write(*,*) 'n_b =', n_b
         write(*,*) 'mu_e=', mu_e         
		 write(*,*) 'mu_n=', mu_n
!
!		 if (electron_pressure(ke) > p_ext) then
!		 p_ext = electron_pressure(ke)
!		 pn = 0.
!		 mu_n = 1.d-20
!		 else
!		 pn = p_ext - electron_pressure(ke)
!		
!         x1=0.0
!         x2=10.
!         xacc=1.d-15
!         kn=root_bisection(kn_solve_pressure,x1,x2,xacc,ierr,hist) !returns in fm**-1     	 
!     	 write(*,*) 'neutron pressure =', neutron_pressure(kn), pn
!		 write(*,*) 'electron pressure =', electron_pressure(ke)
!		 
!		 write(*,*) p_ext, electron_pressure(ke)+neutron_pressure(kn)
!		 write(*,*) p_ext-electron_pressure(ke)-neutron_pressure(kn)	 
!		 !stop
!		 end if		 

 		 ! check stability of distribution against current chemical potentials
		 call check_all_rxns(mu_e, mu_n)
		 	 
	     ! set some dimensions
	     ! make nvar equal to number of entries in mass table + 3 for 3
	     ! lagrange multipliers and constraint equations		 
	     nvar = qt% Ntable +3 
	     neq = nz*nvar
         m1 = (stencil_zones_subdiagonal+1)*nvar-1 ! number of subdiagonals
         m2 = (stencil_zones_superdiagonal+1)*nvar-1  ! number of superdiagonals
         allocate(fac1(qt% Ntable), fac2(qt% Ntable), mu_i(qt% Ntable), ni(qt% Ntable), &
         		& mass_frac(qt% Ntable))
         allocate(equ(nvar,nz), x(nvar,nz), xold(nvar,nz), dx(nvar,nz), xscale(nvar,nz), y(ldy, nsec), stat=ierr)
         if (ierr /= 0) stop 1

      	 ! get mean properties of distribution      
	   	 A_average = 0.
	  	 Z_average = 0.
	  	 Y_sum = 0. 
      	 do j=1,qt% Ntable
      	 Y_sum = (qt% Y(j))+Y_sum
         Z_average = (qt% Y(j))*real(qt% Z(j)) + Z_average
      	 A_average = (qt% Y(j))*real(qt% A(j)) + A_average
      	 enddo
     	 Z_average = Z_average/Y_sum
     	 A_average = A_average/Y_sum      	 
                  
         ! xold is the initial guess for nuclei chemical potentials         
		 ! sets mass fractions to 1d-20 for unpopulated nuclei
		 do j=1,qt% Ntable
!		 mass_frac(j) = real(qt% A(j))*(qt% Y(j))/A_average
		 m_star = mn_n-mp_n-me_n
		 m_nuc = real(qt% A(j))*amu_n
         mterm = g*(m_nuc*kT/(twopi*hbarc_n**2))**(1.5)
!        fac1(j) = real(qt% A(j))/n_b/A_average 
!		 fac1(j) = real(qt% A(j))/n_b
 		 fac1(j) = 1./n_b
         fac2(j) = mterm	
         ! set mass fractions from abundance fractions	 
		 xold(j,1) = log((qt% Y(j))/fac1(j)/fac2(j))*kT-qt%BE(j)
		 !doesnt seem to matter
		 ! try dividing up y equally between all isotopes
		 !xold(j,1) = log((1./(qt% Ntable))/fac1(j)/fac2(j))*kT-qt%BE(j)
		 end do          
         
  		 ! initial values of additional variables 
		! xold(qt% Ntable + 1,1) = n_b
		 xold(qt% Ntable + 1,1) = mu_n 
		
         dx = 0 ! a not very good starting "guess" for the solution
         x = xold

         ! controls
         lid = max_lid
         lrd = max_lrd
         first_step = .true.
         tol_correction_norm=1d-9 ! upper limit on magnitude of average scaled correction
         tol_max_correction = 1d99
         tol_residual_norm = 1d99
         epsder = 1d-6 ! relative variation to compute derivatives
         doing_jacobian = .false.
         matrix_type = square_matrix_type
         
         call newton_work_sizes(m1, m2, nvar, nz, nsec, matrix_type, lwork, liwork, ierr)
         if (ierr /= 0) stop 1
         
         allocate(work(lwork), iwork(liwork), stat=ierr)
         if (ierr /= 0) stop 1
         
         work = 0
         iwork = 0
         
         iwork(r_tol_residual_norm) = 1
         iwork(i_try_really_hard) = 1 ! try really hard for first model
         iwork(i_model_number) = 1
             
         xold = x-dx
         dx = x-xold
         
         
         if (which_decsol == lapack) then
            call do_newt(lapack_decsol, null_decsolblk, null_decsols)
         else if (which_decsol == block_thomas_lapack) then
            call do_newt(null_decsol, block_thomas_lapack_decsolblk, null_decsols)
         else if (which_decsol == mkl_pardiso) then
            call do_newt(null_decsol, null_decsolblk, mkl_pardiso_decsols)
         end if

		 ! if qse solver fails to converge
         if (nonconv) then
         write(*, *) 'failed to converge'
         write(*,*) p_ext, n_b, y_e, y_n, Z_bar, A_bar, mu_e, mu_n
         n_b_prev = n_b
	     end if

         if (iwork(i_debug) /= 0) then
            write(*, *) 'num_jacobians', iwork(i_num_jacobians)
            write(*, *) 'num_solves', iwork(i_num_solves)
         end if
         
         ! if qse solver converges     
		 if (nonconv .eqv. .FALSE.) then
		 write(*,*) 'converged'
  		 A_average = 0.
	  	 Z_average = 0.
	  	 BE_average = 0.
	  	 Y_sum = 0. 
      	 ! get average mass number of distribution 
      	 do j=1,qt% Ntable
      	 Y_sum = (qt% Y(j))+Y_sum
         Z_average = (qt% Y(j))*real(qt% Z(j)) + Z_average
      	 A_average = (qt% Y(j))*real(qt% A(j)) + A_average
      	 BE_average = (qt% Y(j))*(qt% BE(j)) + BE_average
      	 enddo
     	 Z_average = Z_average/Y_sum
     	 A_average = A_average/Y_sum  
     	 BE_average = BE_average/Y_sum
     	 
     	 write(*,*) Y_sum, Z_average, A_average
     	 
     	 !pressure check
         x1=0.0
         x2=10.
         xacc=1.d-15
         kn=root_bisection(kn_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1     	 
     	 write(*,*) 'neutron pressure =', neutron_pressure(kn)
         x1=0.0
         x2=10.
         xacc=1.d-15
         ke=root_bisection(ke_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
		 write(*,*) 'electron pressure =', electron_pressure(ke)
		 
		 write(*,*) p_ext, electron_pressure(ke)+neutron_pressure(kn)
		 write(*,*) p_ext-electron_pressure(ke)-neutron_pressure(kn)	 
		! stop		 
		 		 
		 !call chem_equil_check(mu_n,p_ext,rho,eps_const)
		 !n_b = x(dt% Ntable+1, 1)
         !call get_energy_density(Z_average,A_average,mu_e,mu_n,BE_average,kT,eps_sol)
         eps_sol = 0.
         write(tov_id,'(5(es15.5,2x))') p_ext, n_b*amu_n, eps_sol, mu_e, mu_n
         write(y_output_id,'(8(es12.5,2x))') p_ext, n_b, y_e, y_n, Z_average, A_average, mu_e, mu_n
         write(*,'(8(es12.5,2x))') p_ext, n_b, y_e, y_n, z_average, a_average, mu_e, mu_n
		 n_b_prev = n_b
		 write(*,*) 'eps_check =', eps_const
		 end if
         
         deallocate(iwork, work)
         deallocate(equ, x, xold, dx, xscale, y)
         	
         enddo 
         
        close(y_output_id)  
        
!        do j=1,qt%Ntable
!        write(*,*) qt% Z(j), qt% A(j),qt% Y(j)
!        enddo
         
         contains         
                  
         subroutine do_newt(decsol, decsolblk, decsols)
            interface
               include "mtx_decsol.dek"
               include "mtx_decsolblk.dek"
               include "mtx_decsols.dek"
            end interface
            real*8, pointer :: AF(:,:)
            AF => null()
            call newton( &
               nz, nvar, x, xold, &
               matrix_type, m1, m2, &
               decsol, decsolblk, decsols, &
               lrd, rpar_decsol, lid, ipar_decsol, which_decsol, &
               tol_correction_norm, &
               set_primaries, set_secondaries, default_set_xscale, &
               default_Bdomain, xdomain, eval_equations, &
               size_equ, default_sizeB, default_inspectB, &
               enter_setmatrix, exit_setmatrix, failed_in_setmatrix, default_force_another_iter, &
               xscale, equ, ldy, nsec, y, work, lwork, iwork, liwork, AF, &
               lrpar, rpar, lipar, ipar, &
               nonconv, ierr)
            if (ierr /= 0) stop 1
            deallocate(AF)
         end subroutine do_newt
         
      end subroutine follow_inner_crust_composition

         
      subroutine set_primaries(nvar, nz, x, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: nvar, nz
         real*8, pointer :: x(:,:) ! (nvar, nz)
         integer, intent(in) :: lrpar, lipar
         real*8, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr
         ierr = 0
		 do j = 1, qt% Ntable
		 mu_i(j) = x(j,1)
		 enddo	 
		 lambda_baryon = x(qt% Ntable+1, 1)
		 lambda_charge = x(qt% Ntable+2, 1)
		 lambda_pressure = x(qt% Ntable+3, 1)
      end subroutine set_primaries
      

      subroutine set_secondaries(ivar, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: ivar
         integer, intent(in) :: lrpar, lipar
         real*8, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr
         logical, parameter :: skip_partials = .true.
         call set_sec(ivar, skip_partials, lrpar, rpar, lipar, ipar, ierr)
      end subroutine set_secondaries
      

      subroutine set_sec(ivar, skip_partials, lrpar, rpar, lipar, ipar, ierr)
         use const_def, only:pi
         integer, intent(in) :: ivar
         logical, intent(in) :: skip_partials
         integer, intent(in) :: lrpar, lipar
         real*8, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine set_sec
      

      subroutine eval_equations(iter, nvar, nz, x, xscale, equ, lrpar, rpar, lipar, ipar, ierr)
         use const_def, only:pi
         integer, intent(in) :: iter, nvar, nz
         real*8, pointer, dimension(:,:) :: x, xscale, equ ! (nvar, nz)
         integer, intent(in) :: lrpar, lipar
         real*8, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr
      	 real*8, dimension(nvar*nz, nvar*nz) :: A ! square matrix for jacobian
		 logical, parameter :: skip_partials = .true.			
		 call eval_equ(nvar, nz, equ, skip_partials, A, lrpar, rpar, lipar, ipar, ierr)         
      end subroutine eval_equations
      

      subroutine eval_equ(nvar, nz, equ, skip_partials, A, lrpar, rpar, lipar, ipar, ierr)
         use nucchem_def, only: composition_info_type
         use phys_constants, only: ln10,avogadro,boltzmann
         use electron_eos

         integer, intent(in) :: nvar, nz
         real*8, pointer, dimension(:,:) :: equ
		 logical, intent(in) :: skip_partials
      	 real*8, dimension(nvar*nz, nvar*nz) :: A 
         integer, intent(in) :: lrpar, lipar
         real*8, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr      
         real :: chi, rho      
		 real :: m_star 
		 real :: m_nuc, m_nuc1
		 real :: m_term, m_term1		 
		 real :: ni_Zsum, ni_Asum
		 real :: der_Zsum, der_Asum
		 real :: Asum, Zsum
		 real :: As(nvar), Zs(nvar)
		 real :: Zi, Ai
		 real :: pressure
		 real :: Zbar, Abar
		 real :: n_nin
         real, parameter :: g=1.0d0		 
		 real, parameter :: n_0 = 1.740151d-1
		 real, parameter :: n_1 = -1.577306d-2
         real :: nl, iso    
         real :: phi, phi_sum    
         real :: R_n, R_ws
                
         ierr = 0
		 ! solve for electron wave vector 
		 if (mu_e > 0.) then
         x1=0.0
         x2=10.
         xacc=1.d-15
         ke=root_bisection(ke_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
         if (io_failure(ierr,'Error in bisection for ke wave vector')) stop
         n_e = ke**3/threepisquare               
         Y_e = n_e/n_b   
         else
         mu_e = abs(mu_e)
         x1=0.0
         x2=10.
         xacc=1.d-15
         ke=root_bisection(ke_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
         if (io_failure(ierr,'Error in bisection for ke wave vector')) stop
         n_e = ke**3/threepisquare               
         Y_e = n_e/n_b   
         mu_e = -mu_e
         end if
	     ! solve for neutron wave vector
		 if (mu_n > 0.) then
         x1=0.0
         x2=10.
         xacc=1.d-15
         kn=root_bisection(kn_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
         if (io_failure(ierr,'Error in bisection for kn wave vector')) stop
         n_n = 2.0*kn**3/threepisquare              
		 else
		 mu_n = abs(mu_n)
         x1=0.0
         x2=10.
         xacc=1.d-15
         kn=root_bisection(kn_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
         if (io_failure(ierr,'Error in bisection for kn wave vector')) stop
         n_n = 2.0*kn**3/threepisquare              
		 mu_n = -mu_n
		 end if		

         rho = (n_b*amu_n)*(mev_to_ergs/clight2)/(1.d-39) ! cgs                      

		 ! zero out secondary variables
		 Asum = 0. ; Zsum = 0. 
		 Ai = 0. ; Zi = 0.
		 ni_Asum = 0. ; ni_Zsum = 0.
		 phi_sum=0.
		
		 do i = 1, qt% Ntable
		  iso = 1.0-(2.*real(qt% Z(i))/real(qt% A(i)))
		  n_nin = (n_0+n_1*iso**2)*(1.0+0.9208*iso)/2.	
		  R_n = ((real(qt% A(i))-real(qt% Z(i)))/n_nin/pi*0.75)**onethird 
		  R_ws = (real(qt% Z(i))/n_e/pi*0.75)**onethird		 
          !number density of isotopes
		  m_star = mn_n-mp_n-me_n !does not contain m_e because mu_e has rest mass in definition 
		  m_nuc = real(qt%A(i))*amu_n 
     	  m_term = g*(twopi*hbarc_n**2/(m_nuc*kT))**(-3.0/2.0)
     	  ni(i) = m_term*exp((mu_i(i)+qt%BE(i))/kT)
     	  phi = 1.25*pi*R_n**3*ni(i)
     	  phi_sum = phi+phi_sum
		  !for baryon conservation
		  as(i) = real(qt% A(i))*m_term*exp((mu_i(i)+qt%BE(i))/kT)	 
		  Asum = Asum + as(i) 
		  ni_Asum = ni_Asum + m_term*exp((mu_i(i)+qt%BE(i))/kT)/n_b	
		  Ai = Ai + as(i)/n_b
		  !for charge conservation
		  zs(i) = real(qt% Z(i))*m_term*exp((mu_i(i)+qt%BE(i))/kT)
		  Zsum = Zsum + zs(i)
		  ni_Zsum = ni_Zsum + m_term*exp((mu_i(i)+qt%BE(i))/kT)/n_b	
		  Zi = Zi + zs(i)/n_b
		  !detailed balance
		  equ(i,1) = real(qt% Z(i))*(mu_n-mu_e+m_star)+(real(qt% A(i))-real(qt% Z(i)))*mu_n-mu_i(i) 
		 enddo

		 Zbar = Zi/ni_Zsum
		 Abar = Ai/ni_Asum
		 Z_bar = Zbar
		 A_bar = Abar
 		 chi = phi_sum 
	     Y_n = n_n*(1.-chi)/n_b   
  		 !baryon and charge conservation 
         !equ(qt% Ntable+1,1) = Zsum - n_e
         !equ(qt% Ntable+1,1) = Asum - n_b + n_n*(1.0-chi)  
      end subroutine eval_equ
               
      subroutine eval_jacobian(ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: ldA ! leading dimension of A
         real*8 :: A(ldA, nvar*nz) ! the jacobian matrix
         integer, intent(inout) :: idiag 
         integer, intent(in) :: lrpar, lipar
         real*8, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr         
		 logical, parameter :: skip_partials = .false.			
		 real :: chi
      	 real :: m_term, m_nuc
      	 real, parameter :: g = 1.0d0
      	 real :: asum2, zsum2	
      	 real :: sumn, sume 
      	 real :: Pressure, P_ext		 
         ierr = 0        
      end subroutine eval_jacobian
      
      subroutine enter_setmatrix(iter, nvar, nz, neqs, x, xold, xscale, xder, need_solver_to_eval_jacobian, &
            ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: iter, nvar, nz, neqs
         real*8, pointer, dimension(:,:) :: x, xold, xscale, xder ! (nvar, nz)
         logical, intent(out) :: need_solver_to_eval_jacobian
         integer, intent(in) :: ldA ! leading dimension of A
         real*8, pointer, dimension(:,:) :: A ! (ldA, neqs)
         integer, intent(inout) :: idiag 
         integer, intent(in) :: lrpar, lipar
         real*8, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr
         integer :: i, j
         if (dbg) write(*, '(/, a)') 'enter_setmatrix'
         if (numerical_jacobian) then
            xder=epsder*(xscale+abs(xold))
            need_solver_to_eval_jacobian = .true.
            doing_jacobian = .true.
         else
            call eval_jacobian(ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
            need_solver_to_eval_jacobian = .false.
         end if
      end subroutine enter_setmatrix
      

      subroutine exit_setmatrix(iter, nvar, nz, neqs, dx, ldA, A, idiag, xscale, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: ldA ! leading dimension of A
         integer, intent(in) :: iter, nvar, nz, neqs ! number of equations, 2nd dimension of A
         integer, intent(inout) :: idiag ! row of A with the matrix diagonal entries
         real*8, pointer, dimension(:,:) :: dx, A
         real*8, pointer, dimension(:,:) :: xscale ! (nvar, nz)
         integer, intent(in) :: lrpar, lipar
         real*8, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr
         integer :: i, j
         if (dbg) write(*, '(a, /)') 'exit_setmatrix'
         ierr = 0
         doing_jacobian = .false.
         return
         write(*, *)
         write(*, *) 'result of numerical jacobian'
         do j=1, nvar ! variable
            do i=1, nvar ! equation
               write(*, '(e20.12, 3x)', advance = 'no') A(i, j)
            end do
            write(*, *)
         end do
         write(*, *)
         stop 'exit_setmatrix'
      end subroutine exit_setmatrix


      subroutine failed_in_setmatrix(j, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: j
         integer, intent(in) :: lrpar, lipar
         real*8, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr
         if (dbg) write(*, '(a, /)') 'failed_in_setmatrix'
         ierr = 0
         doing_jacobian = .false.
      end subroutine failed_in_setmatrix

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
  
     function electron_chemical_potential(k) result(mu)
      use phys_constants
      real, intent(in) :: k
      real :: mu
      real :: x   
      x = k*hbarc_n/me_n      
      mu = me_n*(sqrt(1.0+x**2)-1.0)     
     end function electron_chemical_potential

	 function neutron_chemical_potential_local(k) result(mu)
		use phys_constants
		real, intent(in) :: k	! (fm**-1)
		real :: mu	! MeV
		real :: W
	    ! eq. 3.7
	    real, dimension(0:3), parameter :: cw0 = [ 1.2974, 15.0298, -15.2343, 7.4663 ]		
		W = k*(cw0(0) + k*(cw0(1) + k*(cw0(2) + k*cw0(3))))
		mu = W + onethird*  &
				& k*(cw0(0) + k*(2.0*cw0(1) + k*(3.0*cw0(2) + 4.0*k*cw0(3))))
	 end function neutron_chemical_potential_local         
 
     function kn_solve(x)
      real, intent(in) :: x
      real :: kn_solve    
      kn_solve = neutron_chemical_potential_local(x) - mu_n  
     end function kn_solve
     
     function kn_solve_pressure(x)
     	real, intent(in) :: x
     	real :: kn_solve_pressure
     	kn_solve_pressure = neutron_pressure(x) - Pn
     end function kn_solve_pressure	

     function ke_solve(x)
      real, intent(in) :: x
      real :: ke_solve    
      ke_solve = electron_chemical_potential(x) - mu_e  
     end function ke_solve  
     
     function neutron_pressure(k) result(P)
		use phys_constants
		real, intent(in) :: k   ! fm**-1
		real :: P       ! MeV fm**-3
		real :: dWdk
        real, dimension(0:3), parameter :: cw0 = [ 1.2974, 15.0298, -15.2343, 7.4663 ]			
		!dWdk = k*(cw0(0) + k*(2.0*cw0(1) + k*(3.0*cw0(2) + k*4.0*cw0(3)))) 
		dWdk =  (cw0(0) + k*(2.0*cw0(1) + k*(3.0*cw0(2) + k*4.0*cw0(3))))
		P = 2.0*onethird/threepisquare* k**4 * dWdk
     end function neutron_pressure    

     function electron_pressure(k) result(P)
		use phys_constants
		real, intent(in) :: k   ! fm**-1
		real :: P       ! MeV fm**-3
		real :: n			
		n = k**3/threepisquare
		P = 0.25*n*k*hbarc_n
     end function electron_pressure

     function lattice_pressure(Z_average,A_average,n,ne) result(P)
     	use phys_constants
     	real, intent(in) :: Z_average,A_average,n,ne
     	real :: P, p_f
     	real :: aion
     	real, parameter :: C_l = 3.40665d-3
     	p_f = (threepisquare*n)**onethird
     	aion = (Z_average/ne/pi*0.75)**onethird
       !P = -(n/3.0)*C_l*(Z_average**2/A_average**(4.0/3.0))*p_f*hbarc_n
     	P = -(9.0/10.0)*Z_average**2.0/aion*electroncharge**2.0 *hbarc_n*n
     end function lattice_pressure 
     
          
     function neutron_k(x)
      real, intent(in) :: x
      real :: neutron_k      
      neutron_k = neutron_pressure(x) - pres_n     
     end function neutron_k          
     
     subroutine size_equ(iter, nvar, nz, equ, residual_norm, residual_max, lrpar, rpar, lipar, ipar, ierr)
     	integer, intent(in) :: iter, nvar, nz
     	double precision, pointer :: equ(:,:) 
     	double precision, intent(out) :: residual_norm, residual_max
     	integer, intent(in) :: lrpar, lipar
     	double precision, intent(inout) :: rpar(lrpar)
     	integer, intent(inout) :: ipar(lipar)
     	integer, intent(out) :: ierr
     	ierr = 0
     	residual_norm = 1.d-20
     	residual_max = 1.d-15
    end subroutine size_equ 	     
       
      subroutine xdomain(iter, nvar, nz, x, dx, xold, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: iter, nvar, nz
         double precision, pointer, dimension(:,:) :: x, dx, xold ! (nvar, nz)
         integer, intent(in) :: lrpar, lipar
         double precision, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr
         ierr = 0
		 ! set bounds of variables
! 		 x(qt% Ntable+1, 1) = max(n_b_prev, x(qt% Ntable+1,1))
 !		 x(qt% Ntable+1, 1) = max(1.d-7, x(qt% Ntable+1,1))
! 		 x(qt% Ntable+2,1) = abs(x(qt% Ntable+2,1)) 		 
! 		 dx(qt% Ntable+1,1) = x(qt% Ntable+1,1)-xold(qt% Ntable+1,1)     
!		 dx(qt% Ntable+2,1) = x(qt% Ntable+2,1)-xold(qt% Ntable+2,1)     
 !		 x(qt% Ntable+1,1) = xold(qt%Ntable+1,1)+dx(qt%Ntable+1,1)     
!		 x(qt% Ntable+2,1) = xold(qt%Ntable+2,1)+dx(qt%Ntable+2,1)     	 

! 		 if (x(qt% Ntable+1,1) < n_b_prev) then
! 		 x(qt% Ntable+1,1) = n_b_prev
! 		 dx(qt% Ntable+1,1) = x(qt% Ntable+1,1)-xold(qt% Ntable+1,1) 
! 		 x(qt% Ntable+1,1) = xold(qt% Ntable+1,1)+dx(qt% Ntable+1,1) 
! 		 end if
      end subroutine xdomain   
      
      subroutine chem_equil_check(mu,p,rho,eps_const)
      real :: eps_const
      real :: mu, phi_e, phi_g
      real :: n, m_term, psh_n
      real :: p, rho
      real, parameter :: g = 1.d0   
      !real, parameter :: grav = 1.86d14 !cm/s/s   
      real, parameter :: grav = 3.4d-18
      !psh_n = p/rho/grav !kT/mn_n
      !depth = 
 !     phi_e = 0.
      !phi_g = mn_n*depth*grav
 !     phi_g = mn_n*p
!      m_term = g*(mn_n*kT/(twopi*hbarc_n**2))**(1.5)
!	  mu = kT*log(n/m_term)
 !     eps_const = mu + phi_e + phi_g
!      write(*,*) 'mu =', mu, 'phi_g =', phi_g
      end subroutine chem_equil_check
      
      subroutine check_all_rxns(mu_e, mu_n)
      real, intent(in) :: mu_e, mu_n
      real :: pressure
	  real :: Sn, S2n, Sp, S2p, B, ecthresh, bthresh, VN
	  real :: Snr, S2nr, Spr, S2pr, Br, ecthreshr, bthreshr, VNr 
	  real :: alpha(2), beta(2), gamma(2), delta(2), epsilon(2)
	  real :: B_temp, Y_temp, Yr
	  real, parameter :: del_m = mn_n-mp_n-me_n
      integer :: i, j, k, index, Ntable, index_change, index_temp
      integer :: Z, A
	  integer :: Zr, Ar, Nr
	  integer :: Z_temp, A_temp
      integer, dimension(:), allocatable :: Z_new, A_new
      real, dimension(:), allocatable :: B_new, Y_new
	  integer :: ierr, id, ios, iter, iZ, iZb, iZe, iEq(1) 	  
	  integer, parameter :: ineg = 1, ipos = 2
	  integer, parameter :: max_iterations = 100
	  integer, parameter :: qt_temp = 100000
	  integer, parameter :: pycno_thresh_Z = 4
	  logical, save :: dt_table_used = .false.
      logical, dimension(max_iterations) :: neutron_capture, dineutron_capture
	  logical, dimension(max_iterations) :: neutron_emission, dineutron_emission
	  logical, dimension(max_iterations) :: en_rxn, enn_rxn, ne_rxn, nne_rxn
	  character(len=6) :: rxn

      ! set defaults
	  neutron_capture = .FALSE.
	  neutron_emission = .FALSE.
	  dineutron_capture = .FALSE.
	  dineutron_emission = .FALSE.
	  en_rxn = .FALSE.
	  enn_rxn = .FALSE. 
	  ne_rxn = .FALSE.
	  nne_rxn = .FALSE.
      
      ! set size of do loop depending on input table
      if (dt_table_used .eqv. .false.) then
      Ntable = dt% Ntable
	  else
	  Ntable = qt% Ntable
	  end if

	  allocate(Z_new(qt_temp), A_new(qt_temp), B_new(qt_temp), Y_new(qt_temp))
	  index = 1

	  do j = 1, Ntable
	 
	  ierr=0
	  if (dt_table_used .eqv. .false.) then
	  Z = dt% Z(j)
	  A = dt% A(j)
	  B = dt% BE(j)
	  Yr = dt% Y(j)
	  else
	  Z = qt% Z(j)
	  A = qt% A(j)
	  B = qt% BE(j)
	  Yr = qt% Y(j)
	  end if
	  
	  ! add nucleus to new table before it's checked for rxns 
	  Z_new(index) = Z
	  A_new(index) = A
	  B_new(index) = B 
	  Y_new(index) = Yr
	  Z_temp = Z
	  A_temp = A
	  B_temp = B	 
	  Y_temp = Yr 
!      index = index+1
      	 
      ! loop over rxns 
      do iter = 1, max_iterations
   	  !get properties of nucleus that is being pushed deeper
	  call get_nucleus_properties(Z,A,id,B,Sn,S2n,Sp,S2p,ecthresh,bthresh,VN,ierr)
      
!	  if (Z/=Z_new(index) .and. A/=A_new(index) .and. ierr==0 &
!	  	& .and. Z/=Z_temp .or. A/=A_temp) then      
	  if (Z/=Z_temp .or. A/=A_temp) then
	  index = index+1
	  Z_new(index) = Z
	  A_new(index) = A
	  B_new(index) = B
	  Y_new(index) = 0. !if nucleus is new to distribution, Y=0.
	  Z_temp = Z
	  A_temp = A
	  B_temp = B
	  Y_temp = 0.
      end if

      if (index > qt_temp) then
      write(*,*) 'need to allocate more space for qse_table_type'
      stop
      end if
      
      !check for pycnonuclear reactions
      if (Z .le. pycno_thresh_Z) then
	  Zr = 2.0*Z ; Ar = 2.0*A
	  Nr = Ar-Zr
	  call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr, &
	  		& ecthreshr,bthreshr,VNr,ierr)
	  if (ierr /= 0) exit	  
	  rxn = 'pyc'
	  write(*,*) rxn, Z, A, Zr, Ar
	  call print_reaction_check
	  A = Ar ; Z = Zr
	  cycle      
      end if

      ! check for weak reactions                 
      ! electron capture
      Ar = A; Zr = Z-1
      Nr = Ar-Zr
!      index= Zr - mt%Zmin + 1
!      if (Nr < mt% Nmin(index) .or. Nr > mt% Nmax(index)) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr, &
      		&	ecthreshr,bthreshr,VNr,ierr)
      if (ierr /= 0) exit
      alpha(ineg) = del_m - mu_e + (B-Br)
      if (alpha(ineg) < 0) then
         rxn = '(e,)'
         call print_reaction_check
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
      if (ierr /= 0) exit
      alpha(ipos) = mu_e - del_m + (B - Br)
      if (alpha(ipos) < 0) then
         rxn = '(,e)'
         call print_reaction_check
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
      if (ierr /= 0) exit
      beta(ineg) = mu_n - mu_e + del_m + (B - Br)
      if (beta(ineg) < 0) then
      	 en_rxn(iter) = .TRUE.
         rxn = '(e,n)'
         call print_reaction_check
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
      if (ierr /= 0) exit
      beta(ipos) = 2.0*mu_n - mu_e + del_m + (B - Br)
      if (beta(ipos) < 0) then
      	 enn_rxn(iter) = .TRUE.
         rxn = '(e,2n)'
         call print_reaction_check
         A = Ar; Z = Zr
         cycle
      end if
      
      ! check for strong reactions     
      ! neutron capture
      Ar = A+1; Zr = Z
      Nr = Ar-Zr
!      index= Zr - mt%Zmin + 1
!      if (Nr < mt% Nmin(index) .or. Nr > mt% Nmax(index)) exit
      if (Zr < mt% Zmin .or. Zr > mt% Zmax) exit  
      call get_nucleus_properties(Zr,Ar,id,Br,Snr,S2nr,Spr,S2pr,ecthreshr,bthreshr,VNr,ierr)
      if (ierr /= 0) exit
      gamma(ineg) = -mu_n + (B-Br) 
      if (gamma(ineg) < 0) then
         neutron_capture(iter) = .TRUE.
         rxn = '(n,)'
         call print_reaction_check
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
      if (ierr /= 0) exit
      gamma(ipos) = mu_n + (B-Br) 
      if (gamma(ipos) < 0) then
         neutron_emission(iter) = .TRUE.
         rxn = '(,n)'
         call print_reaction_check
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
      if (ierr /= 0) exit
      delta(ineg) = -2.0*mu_n + (B-Br) 
      if (delta(ineg) < 0) then
      	 dineutron_capture(iter) = .TRUE.
         rxn = '(2n,)'
         call print_reaction_check
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
      if (ierr /= 0) exit
      delta(ipos) = 2.0*mu_n + (B-Br) 
      if (delta(ipos) < 0) then
         dineutron_emission(iter) = .TRUE.
         rxn = '(,2n)'
         call print_reaction_check
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
      if (ierr /= 0) exit
      epsilon(ineg) = mu_e - mu_n - del_m + (B-Br) 
      if (epsilon(ineg) < 0) then
      	 ne_rxn(iter) = .TRUE.
         rxn = '(n,e)'
         call print_reaction_check
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
      if (ierr /= 0) exit
      epsilon(ipos) = mu_e - 2.0*mu_n - del_m + (B-Br) 
      if (epsilon(ipos) < 0) then
      	 nne_rxn(iter) = .TRUE.
      	 rxn = '(2n,e)'
      	 call print_reaction_check
      	 A = Ar; Z = Zr
      	 cycle
      end if 
      
      exit
      
      end do  ! end of iteration loop      
      end do ! end of dt% table loop
      
      ! make duplicate entries zero
	  do j=1,index-1
	   do k=1,index-1
	    if (Z_new(k) == Z_new(j) .and. A_new(k) == A_new(j)  &
	    	& .and. j/=k .and. Y_new(j) == 0.) then
	    Z_new(k) = 0
	    A_new(k) = 0
	    B_new(k) = 0.
	    Y_new(k) = 0.
	    end if 
	   end do
	  end do

      ! take care of abundances (need to be nonzero for log space)
      do j=1,index-1
      if (Y_new(j) == 0.) then
      Y_new(j) = 1.d-20
      end if
      end do
	
	  ! find the number of duplicate entries
	  index_change = 0
	  do j = 1, index-1
	  if (Z_new(j) == 0) then
	  index_change = index_change+1
	  end if
	  end do

	  ! allocate memory without space for duplicates
      index_temp = 0
      if (dt_table_used .eqv. .false.) then
      qt% Ntable = index-1-index_change	  	
	  allocate(qt% Z(qt% Ntable), qt% A(qt% Ntable), qt% BE(qt% Ntable), &
	  			& qt% Y(qt% Ntable))
	  else
      call alloc_qse_table(qt% Ntable, index-1-index_change)
	  end if
	  
	  do j=1,index-1
	   if (Z_new(j) == 0) then
	   cycle
	   else
	   index_temp = index_temp + 1
	   qt% Z(index_temp) = Z_new(j)
	   qt% A(index_temp) = A_new(j)
 	   qt% BE(index_temp) = B_new(j)
 	   qt% Y(index_temp) = Y_new(j)	   
	   end if
	  end do  

	  deallocate(Z_new, A_new, B_new, Y_new)
	  
	  ! flip logical for dt% 
	  if (dt_table_used .eqv. .false.) then
	  dt_table_used = .true.
	  call dist_table_shutdown
	  end if
	  
!	  write(*,*) '-------'
!	  write(*,*) qt% Ntable
!	  do j=1,qt% Ntable
!	  write(*,*) qt% Z(j), qt% A(j)
!	  end do
!	  stop

      end subroutine check_all_rxns       

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
		!write (*,form) pressure,mu_n,mu_e,Z,A,rxn,Zr,Ar
	 end subroutine print_reaction
	 	
	 subroutine print_reaction_check()
		character(len=*), parameter :: form = '((a6),2(i6))'
		character(len=6) :: rxn
		!write(*,form) rxn, Z, A
	 end subroutine print_reaction_check
	 
	 subroutine alloc_qse_table(Ntable, Ntable_new)
	 	real, dimension(:), allocatable :: BEi_temp, Yi_temp
		integer, dimension(:), allocatable :: Zi_temp, Ai_temp
		integer :: Ntable, Ntable_new	
		allocate(Zi_temp(Ntable), Ai_temp(Ntable), &
				& BEi_temp(Ntable), Yi_temp(Ntable))
		Zi_temp(1:Ntable) = qt% Z(1:Ntable)
		Ai_temp(1:Ntable) = qt% A(1:Ntable)
		BEi_temp(1:Ntable) = qt% BE(1:Ntable)
		Yi_temp(1:Ntable) = qt% Y(1:Ntable)
		deallocate(qt% Z, qt% A, qt% BE, qt% Y)
		allocate(qt% Z(Ntable_new), qt% A(Ntable_new), &
				& qt% BE(Ntable_new), qt% Y(Ntable_new))
		qt% Z(1:Ntable) = Zi_temp(1:Ntable)
		qt% A(1:Ntable) = Ai_temp(1:Ntable)
		qt% BE(1:Ntable) = BEi_temp(1:Ntable)
		qt% Y(1:Ntable) = Yi_temp(1:Ntable)
		deallocate(Zi_temp, Ai_temp, BEi_temp, Yi_temp)
	 end subroutine
          
    end module inner_crust