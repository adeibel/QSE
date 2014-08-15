      module qse_solver
      use num_def
      use num_lib
      use iso_fortran_env, only: error_unit, oid=>output_unit
      use alert_lib
      use utils_lib
      use phys_constants
      use mass_table 
      use eos_table
      use rootfind      
      use crust_eos_lib      

      implicit none
      
      logical, parameter :: dbg = .false.

      ! dimensions
      integer, parameter :: nz = 1 ! number of zones
      integer, parameter :: nvar = 28+2 !5284+2  ! number of variables per zone
      integer, parameter :: neq = nz*nvar

      ! information about the bandwidth of the jacobian matrix
      ! we have a square matrix, so set number zones sub and super both to 0
      integer, parameter :: stencil_zones_subdiagonal = 0
      integer, parameter :: stencil_zones_superdiagonal = 0
      integer, parameter :: m1 = (stencil_zones_subdiagonal+1)*nvar-1 ! number of subdiagonals
      integer, parameter :: m2 = (stencil_zones_superdiagonal+1)*nvar-1 ! number of superdiagonals

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
      
      integer :: i, j  

	  ! for crust
      type(mass_table_type), pointer :: mt
      type(eos_table_type), pointer :: et
      real*8 :: mu_e, mu_n, mu_i(28)
      real*8 :: mu_e_prev, Y_e_prev
      real*8 :: Y_e, Y_n 
   	  real*8 :: n_b, n_e, n_n
   	  real*8 :: p_ext
      real*8 :: ke, kn
      real :: Z_bar, A_bar
      real,dimension(:),pointer :: hist
      real :: x1, x2, xacc
      integer :: eos_handle
 	  real, save :: zsum_save = 0.
 	  real :: asum, zsum
      real :: kT
      real :: mu_n_prev
      real :: ke_prev, kn_prev
      real :: pres_n
      real :: ni(28)
            real :: Pn
   

                    
      contains
      
      subroutine do_test_newton
         use mtx_lib
         use mtx_def
        
      integer :: ierr, liwork, lwork, lid, lrd, which_decsol
      integer, dimension(:), pointer :: iwork
      real*8, dimension(:), pointer :: work 
      integer, parameter :: lrpar = 0, lipar = 0
      integer, target :: ipar(lipar)
      real*8, target :: rpar(lrpar)
      character (len=64) :: decsol_option_name
      !namelist
      real :: p_ext_start
      real :: n_b_start
      logical :: have_mu_table
      logical :: do_numerical_jacobian
      integer :: which_decsol_in, decsol    
 
 	  ! for crust
 	  character(len=*), parameter :: mass_table_name = 'nucchem_andrew2.data'
 	  character(len=*), parameter :: eos_table_name = 'ashes_acc.txt'
! 	  character(len=*), parameter :: mass_table_name = 'nucchem.data'   
!	  character(len=*), parameter :: mass_table_name = 'nucchem_trunc.data'
	  character(len=*), parameter :: y_output_file = 'y_output.data'
      character(len=*), parameter :: output_file = 'qse_output.data'
   	  character(len=*), parameter :: abundance_file = 'qse_abun.data'
   	  character(len=*), parameter :: default_infile = 'qse.inlist'
   	  character(len=64), parameter :: data_dir = '../../../data/crust_eos'
!	  character(len=64), parameter :: eos_table_name = 'helm_table.dat'
	  character(len=64), parameter :: mu_table_name = 'mu_table_old.data'
      character(len=80) :: infile	
      integer :: i, j, ios
      integer :: inlist_id, output_id, abundance_id  
      integer :: mu_table_input_id, mu_table_output_id
      integer :: y_output_id
      logical, save :: mass_table_is_loaded = .FALSE.
      logical, save :: eos_table_is_loaded = .FALSE.
      logical, save :: ye_set = .FALSE.
      real :: mterm, fac1(28), fac2(28), m_nuc, m_star
      real, parameter :: g = 1.d0

      namelist /range/ P_ext_start, n_b_start, kT, have_mu_table, &
      	&	do_numerical_jacobian, which_decsol_in
     
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
      have_mu_table = .false.
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
 	  write(*,*) 'Loaded mass table'
 	  
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
 	  write(*,*) 'EOS table loaded'
 	  
 !	  write(*,*) et% nb(1), et% rho(1), et% Ye(1), et% Yn(1), et% fr(1), et% Z_bar(1), &
! 	  		& et% A_bar(1), et% Q(1), et% fr_x(1), et% a_ionic(1), et% heat(1), &
! 	  		& et% heat_full(1), et% nn(1), et% pr(1), et% gb(1), et% ne(1), &
! 	  		& et% mun(1), et% mue(1)
! 	  
! 	  stop
   
 	  output_id = alloc_iounit(ierr)
  	  if (io_failure(ierr,'allocating unit for output nuclei file')) stop
  	  open(unit=output_id, file=output_file, iostat=ios, status="unknown")
   	  if (io_failure(ios,'opening output nuclei file')) stop 
	  write(output_id,'(7A13)') 'n_b [fm^-3]', 'k_e [fm^-1]', &
   			& 'k_n [fm^-1]', 'mu_e [MeV]', 'mu_n [MeV]', 'Z_bar', 'A_bar'

	  abundance_id = alloc_iounit(ierr)
	  if (io_failure(ierr, 'allocating unit for abundance file')) stop
	  open(unit=abundance_id, file=abundance_file, iostat=ios, status="unknown")
	  if (io_failure(ios,'opening abundance file')) stop
	  write(abundance_id,'(A13)') 'n_i [fm^-3]'

	  y_output_id = alloc_iounit(ierr)
	  if (io_failure(ierr, 'allocating unit for y fractions file')) stop
	  open(unit=y_output_id, file = y_output_file, iostat=ios, status="unknown")
	  write(y_output_id,'(8(A12),2x)') 'Pressure', 'n_b', 'Ye', 'Yn', 'Zbar', 'Abar', 'mu_e', 'mu_n'

  	  ! solve for qse distribution at each n_b  	  
  	  do i=1, et% Ntable
  	     n_b = et% nb(i) !n_b_start !/1000.
  	     p_ext = (et% pr(i))*hbarc_n !p_ext_start*hbarc_n*real(i)  
		! mu_e = (et% mue(i))*hbarc_n

!		 write(*,*) 'P_ext=', P_ext
!  	     write(*,*) 'n_b =', n_b
!  	     write(*,*) 'numerical jacobian? =', do_numerical_jacobian
!  	     write(*,*) 'have mu table? =', have_mu_table
!         write(*,*) 'decsol option', decsol     
              
         which_decsol = which_decsol_in
         call decsol_option_str(which_decsol, decsol_option_name, ierr)
  !       write(*,*) which_decsol, trim(decsol_option_name), ierr
         if (ierr /= 0) return
!         write(*,*) 'use ' // trim(decsol_option_name)

         if (which_decsol == mkl_pardiso) then
            if (.not. okay_to_use_mkl_pardiso()) which_decsol = lapack
         end if
         
         allocate(equ(nvar,nz), x(nvar,nz), xold(nvar,nz), dx(nvar,nz), xscale(nvar,nz), y(ldy, nsec), stat=ierr)
         if (ierr /= 0) stop 1

55 continue 

         numerical_jacobian = do_numerical_jacobian
         if (have_mu_table .eqv. .true.) then
		 mu_table_input_id = alloc_iounit(ierr)
		 if (io_failure(ierr, 'allocating unit for mu table read')) stop		 
		 open(mu_table_input_id,file=mu_table_name, iostat=ios, status='unknown')
		 do j = 1, mt% Ntable + 2
		 read(mu_table_input_id,*) xold(j,1)
		 enddo		 
		 close(mu_table_input_id)
         call free_iounit(mu_table_input_id)		 
		 else 


!         xold(mt% Ntable+1,1) = (log(fac1*fac2)*kT+real(mt% Z(867))*m_star&
!         	& + mt%BE(867)+real(mt%A(867))*xold(mt% Ntable+2,1))/real(mt% Z(867))
!         xold(mt% Ntable+1,1) = (log(fac1*fac2)*kT+real(mt% Z(867))*m_star&
!        	& +real(mt%A(867))*xold(mt% Ntable+3,1))/real(mt% Z(867))

		 ! sets mass fractions to 1d-30
		 do j=1,mt% Ntable
		 m_star = mn_n-mp_n-me_n
		 m_nuc = real(mt% A(j))*amu_n
         mterm = g*(m_nuc*kT/(twopi*hbarc_n**2))**(1.5)
         fac1(j) = real(mt% A(j))/n_b
         fac2(j) = mterm		 
		 !xold(j,1) = log((1./16.)/fac1(j)/fac2(j))*kT-mt%BE(j)
		 xold(j,1) = log((1.d-20)/fac1(j)/fac2(j))*kT-mt%BE(j)
		 end do 
		 
		 ! set mass fraction to 1.0 for one nucleus
		 xold(14,1) = log((1.d0)/fac1(14)/fac2(14))*kT-mt%BE(14)
		 !xold(8,1) =  log((1.d0)/fac1(8)/fac2(8))*kT-mt%BE(8)
	!	 xold(mt% Ntable+1,1) = -(xold(8,1))/(26.0) + m_star
		 
		 !mu_e = -(xold(8,1))/26 + m_star		 
		 mu_e = -(xold(14,1))/28. + m_star
		 !mu_e = (et% mue(i))*hbarc_n
		 !xold(mt% Ntable+2,1) = n_b
		 mu_n = mu_e/1000.

		! initial values of additional variables 
		! xold(mt% Ntable+1,1) = n_b
		 xold(mt% Ntable+2,1) = mu_n 
		 xold(mt% Ntable+1,1) = mu_e
			
!		 write(*,*) mu_e, mu_n, n_b
		 
		 end if

         dx = 0 ! a not very good starting "guess" for the solution
         !dx = x/1.d-3
         x = xold
         
         lid = max_lid
         lrd = max_lrd
                  
         first_step = .true.
         
         tol_correction_norm=1d-9
         !tol_correction_norm = 1d-12
         !tol_correction_norm = 1d-14 ! upper limit on magnitude of average scaled correction
         tol_max_correction = 1d99
         tol_residual_norm = 1d99
         
         !epsder = 1d-8
         epsder = 1d-6 ! relative variation to compute derivatives
         !epsder = p_ext
         
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
             
         xold = x-dx; dx = x-xold
         
         if (which_decsol == lapack) then
            call do_newt(lapack_decsol, null_decsolblk, null_decsols)
         else if (which_decsol == block_thomas_lapack) then
            call do_newt(null_decsol, block_thomas_lapack_decsolblk, null_decsols)
         else if (which_decsol == mkl_pardiso) then
            call do_newt(null_decsol, null_decsolblk, mkl_pardiso_decsols)
         end if

         if (nonconv) then
            write(*, *) 'failed to converge'
        write(*,*) p_ext, n_b, y_e, y_n, Z_bar, A_bar, mu_e, mu_n
            
	     endif
!			mu_table_output_id = alloc_iounit(ierr)
!	  		if (io_failure(ierr, 'allocating unit for mu table file for output')) stop
!	        open(unit=mu_table_output_id, file=mu_table_name, iostat=ios, status="unknown")
!	        if (io_failure(ios,'opening mu table file for output')) stop
!	        do j = 1,mt%Ntable
!	        write(mu_table_output_id,*) mu_i(j)
!	        enddo 
!	        write(mu_table_output_id,*) Y_e
!	        write(mu_table_output_id,*) Y_n
!	        close(mu_table_output_id) 
!	        call free_iounit(mu_table_output_id) 
!	        have_mu_table = .true. 
!	        goto 55
 !           stop 2
 !        end if

         if (iwork(i_debug) /= 0) then
            write(*, *) 'num_jacobians', iwork(i_num_jacobians)
            write(*, *) 'num_solves', iwork(i_num_solves)
         end if
         
         deallocate(iwork, work)
         deallocate(equ, x, xold, dx, xscale, y)
         
      
!         write(*,*) 'finished n_b', i

!		do j = 1, mt% Ntable
!		write(*,*) mt% Z(j), mt% A(j), ni(j)
!		enddo

		if (nonconv .eqv. .FALSE.) then
        write(y_output_id,'(8(es12.5,2x))') p_ext, n_b, y_e, y_n, Z_bar, A_bar, mu_e, mu_n
		end if
	
 !        stop
         
       !  have_mu_table = .true. 
  		 have_mu_table = .false.
  
         enddo 
         
        close(y_output_id)  
         
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
         
      end subroutine do_test_newton

         
      subroutine set_primaries(nvar, nz, x, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: nvar, nz
         real*8, pointer :: x(:,:) ! (nvar, nz)
         integer, intent(in) :: lrpar, lipar
         real*8, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr
         ierr = 0
		 do i = 1, mt% Ntable
		 mu_i(i) = x(i,1)
		 enddo
		 mu_e = x(mt% Ntable+1,1)	 
		 mu_n = x(mt% Ntable+2,1)
		 !n_b = x(mt% Ntable+1,1)			 
      	 ! p_ext = x(mt% Ntable+3,1)
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
		 !for loop over nuclei abundances
         real, parameter :: g=1.0d0
		 real :: m_star 
		 real :: m_nuc, m_nuc1
		 real :: m_term, m_term1		 
		 !for equations in log space
		 real :: sum_lnZ(16), sum_lnA(16) 
		 real :: sum_lnZ_total, sum_lnZ_final 
		 real :: sum_lnA_total, sum_lnA_final
		 real :: logZ_exponent
		 real :: logA_exponent 
		 real :: ni_Zsum, ni_Asum
		 real :: der_Zsum, der_Asum
		 real :: Asum, Zsum
		 real :: As(28), Zs(28)
		 real :: Zi, Ai
		 real :: pressure
		 real :: Zbar, Abar
!		 real :: ni(16)
		 ! for chi
		 real :: n_nin
		 real, parameter :: n_0 = 1.740151d-1
		 real, parameter :: n_1 = -1.577306d-2
         real :: nl       
         real :: iso    
         real :: phi, phi_sum    
         real :: R_n, R_ws
                
         ierr = 0

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

		if (mu_n > 0.) then
        x1=0.0
        x2=10.
        xacc=1.d-15
        kn=root_bisection(kn_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
        if (io_failure(ierr,'Error in bisection for kn wave vector')) stop
        n_n = 2.0*kn**3/threepisquare              
        !Y_n = n_n*(1.-chi)/n_b   
		Y_n = n_n/n_b
		else
		mu_n = abs(mu_n)
        x1=0.0
        x2=10.
        xacc=1.d-15
        kn=root_bisection(kn_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
        if (io_failure(ierr,'Error in bisection for kn wave vector')) stop
        n_n = 2.0*kn**3/threepisquare              
        
        !Y_n = n_n*(1.-chi)/n_b   
		Y_n = n_n/n_b
		mu_n = -mu_n
		end if		

	     !chi = use_default_nuclear_size
	     !chi = 0.95
         rho = (n_b*amu_n)*(mev_to_ergs/clight2)/(1.d-39) ! cgs                      

		! enter with mu_n and mu_e from initial conditions
		
		 Asum = 0. ; Zsum = 0. 
		 Ai = 0. ; Zi = 0.
		 ni_Asum = 0. ; ni_Zsum = 0.
		 phi_sum=0.

		
		 do i = 1, mt% Ntable
		  iso = 1.0-(2.*real(mt% Z(i))/real(mt% A(i)))
		  n_nin = (n_0+n_1*iso**2)*(1.0+0.9208*iso)/2.	
		  R_n = (real(mt% N(i))/n_nin/pi*0.75)**onethird 
		  R_ws = (real(mt% Z(i))/n_e/pi*0.75)**onethird		 
          !number density of isotopes
		  m_star = mn_n-mp_n-me_n !does not contain m_e because mu_e has rest mass in definition 
		  m_nuc = real(mt%A(i))*amu_n !real(mt% Z(i))*mp_n+real(mt% N(i))*mn_n !-mt%BE(i)    
     	  m_term = g*(twopi*hbarc_n**2/(m_nuc*kT))**(-3.0/2.0)
     	  ni(i) = m_term*exp((mu_i(i)+mt%BE(i))/kT)
     	  phi = 1.25*pi*R_n**3*ni(i)
     	  !phi = (R_n/R_ws)**3
!     	  write(*,*) 'phi=', phi
     	  phi_sum = phi+phi_sum
		  !for baryon conservation
		  as(i) = real(mt% A(i))*m_term*exp((mu_i(i)+mt%BE(i))/kT)	 
		  Asum = Asum + as(i) 
		  ni_Asum = ni_Asum + m_term*exp((mu_i(i)+mt%BE(i))/kT)/n_b	
		  Ai = Ai + as(i)/n_b
		  !for charge conservation
		  zs(i) = real(mt% Z(i))*m_term*exp((mu_i(i)+mt%BE(i))/kT)
		  Zsum = Zsum + zs(i)
		  ni_Zsum = ni_Zsum + m_term*exp((mu_i(i)+mt%BE(i))/kT)/n_b	
		  Zi = Zi + zs(i)/n_b
		  !detailed balance
		  equ(i,1) = real(mt% Z(i))*(mu_n-mu_e+m_star)+real(mt% N(i))*mu_n-mu_i(i) !+mt%BE(i)
!		  write(*,*) equ(i,1)
		 enddo
		 
		! chi = phi_sum

		 Zbar = Zi/ni_Zsum
		 Abar = Ai/ni_Asum
		 Z_bar = Zbar
		 A_bar = Abar

		 !chi = use_default_nuclear_size
	     !chi = 3.d-3
		 chi = phi_sum 
		 ! chi = use_default_nuclear_size

	!	 if (Zsum > 0.) then
!		 n_e = Zsum
!		 ke = (n_e*threepisquare)**onethird
!		 mu_e = electron_chemical_potential(ke)
	!	 end if

		 !Pn = P_ext - electron_pressure(ke) - lattice_pressure(Zbar,Abar,n_b,n_e)
!		 n_n = (n_b-Asum)/(1.-chi)
!		 kn = (threepisquare*n_n/2.0)**onethird
!		 mu_n = neutron_chemical_potential(kn)

!		if (pn > 0.) then
!        x1=0.0
!        x2=10.
!        xacc=1.d-15
!        kn=root_bisection(kn_solve_pressure,x1,x2,xacc,ierr,hist) !returns in fm**-1
!        if (io_failure(ierr,'Error in bisection for kn wave vector')) then
!        write(*,*) 'Pn=', Pn
!        end if
!
!		mu_n = neutron_chemical_potential(kn)
!		n_n = 2.0*kn**3/threepisquare
!		else
!		mu_n = 0.
!		n_n = 0.
!		end if


		 !if (Asum < n_b) then		 
!		 n_n = (n_b - Asum)/(1.0-chi)
!		 kn = (n_n*threepisquare/2.0)**onethird
!		 mu_n = neutron_chemical_potential(kn)
		 !end if

		  !equ(mt% Ntable+1, 1) = chi*Asum - n_b + n_n*(1.0-chi)
		!  equ(mt% Ntable+1, 1) = Asum - n_b + n_n !*(1.0-chi)
!		 equ(mt% Ntable+1,1) = neutron_pressure(kn)+electron_pressure(ke) &
		 	!& - P_ext
!		 	& +lattice_pressure(Zbar,Abar,n_b,n_e) - P_ext
		 
	!	 y_e = n_e/n_b

		  
  		 !baryon and charge conservation 
        equ(mt% Ntable+1,1) = Zsum - n_e

         equ(mt% Ntable+2,1) = Asum - n_b + n_n !*(1.0-chi)  
 !	     equ(mt% Ntable+1, 1) = Zsum/n_b - n_e/n_b
 	    ! equ(mt% Ntable+3, 1) = chi*ni_Asum - 1.0 + y_n 
!        equ(mt% Ntable+3, 1) = Ai - 1.0 + y_n
 !        equ(mt% Ntable+1, 1) = electron_pressure(ke)+neutron_pressure(kn) &
            !& - P_ext
!	 		& + lattice_pressure(Zbar,Abar,n_b,n_e) - P_ext

!		write(*,*) 'mu_e =', mu_e, 'mu_n=', mu_n
!		write(*,*) 'Asum=', Asum, 'n_b =', n_b, 'n_n*(1-chi)=', n_n*(1.0-chi)
!		write(*,*) 'chi=', chi
!		write(*,*) 'Zsum =', Zsum
!		write(*,*) 'Zsum-ne = ' , Zsum - n_e
!		write(*,*) 'Abar=', Abar, 'Zbar=', Zbar
!		write(*,*) 'n_b=', n_b
! 		write(*,*) 'Y_e=', Y_e, 'mu_e=', mu_e, 'n_e=', n_e, 'ke=', ke
!        write(*,*) 'Y_n=', Y_n, 'mu_n=', mu_n, 'n_n=', n_n, 'kn=', kn
!	    write(*,*) 'mu_i', mu_i(1), mu_i(16), equ(1,1), equ(16,1)
!!	    write(*,*) 'sumZ=', Zsum, 'n_e=', n_e, 'equN_1=', equ(mt% Ntable+1,1)
!!	    write(*,*) 'sumA=', Asum, 'n_b=', n_b, 'n_n=', n_n,  'equN_2=', equ(mt% Ntable+3,1)
!		write(*,*) 'pressure=', electron_pressure(ke) + neutron_pressure(kn), &
!			& 'P_ext=', P_ext, 'equN_3=', equ(mt% Ntable+1,1)
!		write(*,*) 'Pressure [fm-4] = ', p_ext/hbarc_n
!		write(*,*) 'lattice_pressure=', lattice_pressure(Zi/ni_Zsum,Ai/ni_Asum,n_b*ni_Asum/Ai,n_e)
!		write(*,*) 'lattice_pressure(n_b)=', lattice_pressure(Zi/ni_Zsum,Ai/ni_Asum,n_b,n_e)
!	    write(*,*) 'Zi=', Zi/ni_Zsum, 'Ai=', Ai/ni_Asum
!     	write(*,*) '------------------------------'                   
 	
      end subroutine eval_equ
      
      
      subroutine eval_jacobian(ldA, A, idiag, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: ldA ! leading dimension of A
         real*8 :: A(ldA, nvar*nz) ! the jacobian matrix
         ! A(idiag+q-v, v) = partial of equation(q) wrt variable(v)
         integer, intent(inout) :: idiag 
         integer, intent(in) :: lrpar, lipar
         real*8, intent(inout) :: rpar(lrpar)
         integer, intent(inout) :: ipar(lipar)
         integer, intent(out) :: ierr         
		 logical, parameter :: skip_partials = .false.			
		 real :: chi
      	!real :: mu_e_prev, mu_n_prev
      	!real :: ke_prev, kn_prev
      	real :: m_term, m_nuc
      	real, parameter :: g = 1.0d0
      	real :: asum2, zsum2	
      	real :: sumn, sume 
      	real :: Pressure, P_ext
		 
         ierr = 0        

	     chi = use_default_nuclear_size
         !rho = (n_b*amu_n)*(mev_to_ergs/clight2)/(1.d-39) ! cgs
            
		if (mu_e < 0. ) then
		mu_e = abs(mu_e)
		! electron wave vector fm^-1
        x1=0.0
        x2=10.
        xacc=1.d-15
        ke=root_bisection(ke_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
        if (io_failure(ierr,'Error in bisection for ke wave vector')) then
        write(*,*) 'mu_e < 0', mu_e
        ke = ke_prev
        mu_e = mu_e_prev
        end if
        ke=ke/hbarc_n
        n_e = ke**3/threepisquare               
        Y_e = n_e/n_b   
        mu_e = -abs(mu_e)
		end if
		
		if (mu_e > 0.) then                
        ! electron wave vector fm^-1
        x1=0.0
        x2=10.
        xacc=1.d-15
        ke=root_bisection(ke_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
        if (io_failure(ierr,'Error in bisection for ke wave vector')) then
        write(*,*) 'mu_e > 0', mu_e, ke
        ke= ke_prev
        mu_e = mu_e_prev
        end if
        ke = ke/hbarc_n
        n_e = ke**3/threepisquare               
        Y_e = n_e/n_b   
        end if
        
        if (mu_e == 0.) then
        ke = 0. 
        n_e = 0. 
        Y_e = 0.
        end if

		if (mu_n < 0.) then
		mu_n = abs(mu_n)
        x1=0.0
        x2=10.
        xacc=1.d-15
        kn=root_bisection(kn_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
        if (io_failure(ierr,'Error in bisection for kn wave vector')) then
        write(*,*) 'mu_n<0', 'mu_n=', mu_n, kn
        kn = kn_prev
        mu_n = mu_n_prev
        end if   
        kn = kn/hbarc_n
        n_n = 2.0*kn**3/threepisquare !-2.0*kn**3/threepisquare               
        Y_n = n_n*(1.-chi)/n_b   
		mu_n = -abs(mu_n) 
		end if
		
		if (mu_n > 0.) then
        x1=0.0
        x2=10.
        xacc=1.d-15
        kn=root_bisection(kn_solve,x1,x2,xacc,ierr,hist) !returns in fm**-1
        if (io_failure(ierr,'Error in bisection for kn wave vector')) then
        write(*,*) 'mu_n>0', 'mu_n=', mu_n, kn
        kn = kn_prev
        mu_n = mu_n_prev
        end if
        kn = kn/hbarc_n
        n_n = 2.0*kn**3/threepisquare              
        Y_n = n_n*(1.-chi)/n_b   
		!Y_n = n_n/n_b
		end if

		if (mu_n == 0.) then
		kn=0.
		n_n = 0.
		Y_n = 0.
		end if 
		 
		 do i=1, mt% Ntable
		 	do j=1, mt% Ntable
		 	 ! diagonal jacobian => d(equ)/dmu_i	 	
		 	 if(i==j) then
		 	 A(i,j)= -1.0
		 	 else
		 	 A(i,j)=0.0
		 	 endif
		    enddo
		 enddo
		 
		 asum2 = 0. ; zsum2 = 0.
		 sumn = 0. ; sume = 0.
		 
		 do i=1, mt% Ntable    
		 
		  m_nuc = real(mt% Z(i))*mp_n + real(mt% N(i))*mn_n         
     	  m_term = g*(twopi*hbarc_n**2/(m_nuc*kT))**(-3.0/2.0)
		    
  		    !last two rows of jacobian, derivatives wrt the conservation equations 
     		! n=0 term
     		A(mt% Ntable+1, i) = real(mt% Z(i))*m_term*exp((mu_i(i)+mt%BE(i))/kT)/kT/n_b	
     		A(mt% Ntable+2, i) = real(mt% A(i))*m_term*exp((mu_i(i)+mt%BE(i))/kT)/kT/n_b	

	asum2 = asum2 + real(mt% A(i))**2*m_term*exp((mu_i(i)+mt%BE(i))/kT)/kT/n_b
	zsum2 = zsum2 + real(mt% Z(i))**2*m_term*exp((mu_i(i)+mt%BE(i))/kT)/kT/n_b
	sume = sume + real(mt% Z(i))*real(mt%A(i))*m_term*exp((mu_i(i)+mt%BE(i))/kT)/kT/n_b

		    !last two columns 
			A(i, mt% Ntable+1) = 0. !-real(mt% Z(i)) !-2.*real(mt% Z(i))		 ! MeV		    
			A(i, mt% Ntable+2) = 0. !-real(mt% A(i))  !real(mt% A(i)) 		 ! MeV

		end do

			A(mt% Ntable+1, mt% Ntable+1) = -zsum2			
			A(mt% Ntable+1, mt% Ntable+2) = sume
			A(mt% Ntable+2, mt% Ntable+2) = asum2 
			A(mt% Ntable+2, mt% Ntable+1) = -sume
			
			A(:, mt% Ntable+3) = 0. 

!			A(mt% Ntable+1, mt% Ntable+1) = 0.		
!			A(mt% Ntable+1, mt% Ntable+2) = 0.
!			A(mt% Ntable+2, mt% Ntable+2) = 0. 
!			A(mt% Ntable+2, mt% Ntable+1) = 0.
			
		do i = 1, mt% Ntable+3
		A(:, i) = A(:, i)*xscale(i,1)
		end do	
			

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

	 function neutron_chemical_potential(k) result(mu)
		use phys_constants
		real, intent(in) :: k	! (fm**-1)
		real :: mu	! MeV
		real :: W
	    ! eq. 3.7
	    real, dimension(0:3), parameter :: cw0 = [ 1.2974, 15.0298, -15.2343, 7.4663 ]		
		W = k*(cw0(0) + k*(cw0(1) + k*(cw0(2) + k*cw0(3))))
		mu = W + onethird*  &
				& k*(cw0(0) + k*(2.0*cw0(1) + k*(3.0*cw0(2) + 4.0*k*cw0(3))))
	 end function neutron_chemical_potential            
 
     function kn_solve(x)
      real, intent(in) :: x
      real :: kn_solve    
      kn_solve = neutron_chemical_potential(x) - mu_n  
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
         ! if the variables must be non-negative, you can do the following:
         !x = max(0d0, x)
         !dx = x-xold
         !x = xold+dx
 		 ! set mu_e, mu_n, and n_b >0
 		 x(mt% Ntable+1,1) = abs(x(mt% Ntable+1,1))
 		 x(mt% Ntable+2,1) = abs(x(mt% Ntable+2,1))
  !   	 x(mt% Ntable+3,1) = abs(x(mt% Ntable+3,1))
 		 dx(mt% Ntable+1,1) = x(mt% Ntable+1,1)-xold(mt% Ntable+1,1)     
		 dx(mt% Ntable+2,1) = x(mt% Ntable+2,1)-xold(mt% Ntable+2,1)     
 !		 dx(mt% Ntable+3,1) = x(mt% Ntable+3,1)-xold(mt% Ntable+3,1)
 		 x(mt% Ntable+1,1) = xold(mt%Ntable+1,1)+dx(mt%Ntable+1,1)     
		 x(mt% Ntable+2,1) = xold(mt%Ntable+2,1)+dx(mt%Ntable+2,1)     
!		 x(mt% Ntable+3,1) = xold(mt%Ntable+3,1)+dx(mt%Ntable+3,1)     
 		 
 		 ! set mu_n<0 in the outer crust
 !		 x(mt% Ntable+2,1) = -abs(x(mt%Ntable+2,1))
 !		 dx(mt% Ntable+2,1) = x(mt% Ntable+2,1)-xold(mt% Ntable+2,1) 
 !		 x(mt% Ntable+2,1) = xold(mt%Ntable+2,1)+dx(mt%Ntable+2,1)  
 		 
 		 ! set mu_e > 20.
! 		 x(mt% Ntable+1, 1) = max(20.,x(mt% Ntable+1,1)) 
! 		 dx(mt% Ntable+1, 1) =  x(mt% Ntable+1,1)-xold(mt% Ntable+1,1) 
! 		 x(mt% Ntable+1,1) = xold(mt%Ntable+1,1)+dx(mt%Ntable+1,1)  
      end subroutine xdomain          
          
    end module qse_solver
