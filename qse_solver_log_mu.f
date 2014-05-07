      module qse_solver
      use num_def
      use num_lib
      use iso_fortran_env, only: error_unit, oid=>output_unit
      use alert_lib
      use utils_lib
      use phys_constants
      use mass_table 
      use rootfind      
      use crust_eos_lib       

      implicit none
      
      logical, parameter :: dbg = .false.

      ! dimensions
      integer, parameter :: nz = 1 ! number of zones
      integer, parameter :: nvar = 5549+2  ! number of variables per zone
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
      real*8 :: mu_e, mu_n, mu_i(5549), n_i(5549)
      real*8 :: Y_e, Y_n 
   	  real*8 :: n_b, n_e, n_n
      real*8 :: ke, kn
      real,dimension(:),pointer :: hist
      real :: x1, x2, xacc
      integer :: eos_handle
 	  real, save :: zsum_save = 0.
 	  real :: asum, zsum
      real :: kT
		 real :: ke_prev, kn_prev
		 real :: mu_e_prev, mu_n_prev
            
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
      real :: n_b_start
      logical :: have_mu_table
      logical :: do_numerical_jacobian
      integer :: which_decsol_in, decsol    
 
 	  ! for crust
 	  character(len=*), parameter :: mass_table_name = 'nucchem.data'   
      character(len=*), parameter :: output_file = 'qse_output.data'
   	  character(len=*), parameter :: abundance_file = 'qse_abun.data'
   	  character(len=*), parameter :: default_infile = 'qse.inlist'
   	  character(len=64), parameter :: data_dir = '../../../data/crust_eos'
	  character(len=64), parameter :: eos_table_name = 'helm_table.dat'
	  character(len=64), parameter :: mu_table_name = 'mu_table_old.data'
      character(len=80) :: infile	
      integer :: i, j, ios
      integer :: inlist_id, output_id, abundance_id  
      integer :: mu_table_input_id, mu_table_output_id
      logical, save :: mass_table_is_loaded = .FALSE.
      real :: mterm, fac1, fac2, m_nuc, m_star
      real, parameter :: g = 1.d0


      namelist /range/ n_b_start, kT, have_mu_table, &
      	&	do_numerical_jacobian, which_decsol_in
     
      ! set the name of the inlist file name  
      if (command_argument_count() == 1) then
       call get_command_argument(1, infile)
      else
       infile = default_infile
      end if
      
      ! set defaults 
      n_b_start = 5.0d-10 !fm^-3
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

  	  ! solve for qse distribution at each n_b  	  
  	  do i=1,1
  	     n_b = n_b_start*real(i)  

  	     write(*,*) 'n_b =', n_b
  	     write(*,*) 'numerical jacobian? =', do_numerical_jacobian
  	     write(*,*) 'have mu table? =', have_mu_table
         write(*,*) 'decsol option', decsol     
              
         which_decsol = which_decsol_in
         call decsol_option_str(which_decsol, decsol_option_name, ierr)
         write(*,*) which_decsol, trim(decsol_option_name), ierr
         if (ierr /= 0) return
         write(*,*) 'use ' // trim(decsol_option_name)

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

 		 xold(mt% Ntable+2,1) = 0.
		 m_star = mn_n-mp_n-me_n
		 m_nuc = real(mt% A(867))*amu_n  		         
         mterm = g*(m_nuc*kT/(twopi*hbarc_n**2))**(1.5)
         fac1 = real(mt% A(867))/n_b
         fac2 = mterm
!         xold(mt% Ntable+1,1) = (log(fac1*fac2)*kT+real(mt% Z(867))*m_star&
!         	& + mt%BE(867)+real(mt%A(867))*xold(mt% Ntable+2,1))/real(mt% Z(867))
         xold(mt% Ntable+1,1) = (log(fac1*fac2)*kT+real(mt% Z(867))*m_star&
        	& +real(mt%A(867))*xold(mt% Ntable+2,1))/real(mt% Z(867))
		 
		 
		 do j=1,mt% Ntable
		 xold(j,1) = -mt% BE(j) !- 1.0
		 end do		 
		 
	 
		 end if

         dx = 0 ! a not very good starting "guess" for the solution
         x = xold
         
         lid = max_lid
         lrd = max_lrd
                  
         first_step = .true.
         
         tol_correction_norm = 1d-9 ! upper limit on magnitude of average scaled correction
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
			mu_table_output_id = alloc_iounit(ierr)
	  		if (io_failure(ierr, 'allocating unit for mu table file for output')) stop
	        open(unit=mu_table_output_id, file=mu_table_name, iostat=ios, status="unknown")
	        if (io_failure(ios,'opening mu table file for output')) stop
	        do j = 1,mt%Ntable
	        write(mu_table_output_id,*) mu_i(j)
	        enddo 
	        write(mu_table_output_id,*) Y_e
	        write(mu_table_output_id,*) Y_n
	        close(mu_table_output_id) 
	        call free_iounit(mu_table_output_id) 
	        have_mu_table = .true. 
	        goto 55
            stop 2
         end if

         if (iwork(i_debug) /= 0) then
            write(*, *) 'num_jacobians', iwork(i_num_jacobians)
            write(*, *) 'num_solves', iwork(i_num_solves)
         end if
         
         deallocate(iwork, work)
         deallocate(equ, x, xold, dx, xscale, y)
         
         if (nonconv) then !stop 1
         have_mu_table = .true.
         goto 55
         end if
         
         write(*,*) 'finished n_b', i
         have_mu_table = .true. 
         
 		write(*,*) 'mu_e=', mu_e
 		write(*,*) 'n_e=', n_e 
        write(*,*) 'mu_n=', mu_n
        write(*,*) 'n_n=', n_n 
        write(*,*) 'Y_n=', Y_n
        write(*,*) 'Y_e=', Y_e
		write(*,*) mu_i(1)
	    write(*,*) equ(1,1), kn, ke, n_e, n_n
	    write(*,*) equ(mt% Ntable+1,1), equ(mt% Ntable+2,1)
                  
         stop
         
         enddo 
         
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
               default_Bdomain, default_xdomain, eval_equations, &
               default_size_equ, default_sizeB, default_inspectB, &
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
		 mu_e = x((mt% Ntable)+1,1)		 
		 mu_n = x((mt% Ntable)+2,1)
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
		 !for analytical jacobian
	     real, dimension(0:3), parameter :: cw0 = [ 1.2974, 15.0298, -15.2343, 7.4663 ]		
		 real :: dmudk_n, dmudk_e, dkdn_n, dkdn_e
		 real :: dmudne, dmudnn
		 !for equations in log space
		 real :: sum_lnZ(5549), sum_lnA(5549) 
		 real :: sum_lnZ_total, sum_lnZ_final 
		 real :: sum_lnA_total, sum_lnA_final
		 real :: logZ_exponent
		 real :: logA_exponent 
		 real :: ni_Zsum, ni_Asum
		 real :: n_i(5549)


         ierr = 0
	     chi = use_default_nuclear_size
         rho = (n_b*amu_n)*(mev_to_ergs/clight2)/(1.d-39) ! cgs

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
        n_e = ke**3/threepisquare               
        Y_e = n_e/n_b   
        mu_e = -mu_e
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
        n_n = 2.0*kn**3/threepisquare              
        Y_n = n_n*(1.-chi)/n_b   
		!Y_n = n_n/n_b
		end if
	
		if (mu_n == 0.) then
		kn=0.
		n_n = 0.
		Y_n = 0.
		end if 

		 sum_lnA = 0. ; sum_lnA_total = 0. ; sum_lnA_final = 0. 
		 sum_lnZ = 0. ; sum_lnZ_total = 0. ; sum_lnZ_final = 0. 

		 do i = 2, mt% Ntable
          !number density of isotopes
		  m_star = mn_n-mp_n-me_n !does not contain m_e because mu_e has rest mass in definition 
		  m_nuc = real(mt% Z(i))*mp_n + real(mt% N(i))*mn_n         
     	  m_term = g*(twopi*hbarc_n**2/(m_nuc*kT))**(-3.0/2.0)
     	  m_nuc1 = real(mt% Z(1))*mp_n + real(mt% N(1))*mn_n 
     	  m_term1 = g*(twopi*hbarc_n**2/(m_nuc1*kT))**(-3.0/2.0)
		  !for baryon conservation
		  sum_lnA(i) = log(real(mt%A(i))*m_term) + (mu_i(i)+abs(mt%BE(i)))/kT
		  sum_lnA(1) = log(real(mt%A(1))*m_term1) + (mu_i(1)+abs(mt%BE(1)))/kT
		  sum_lnA(i) = exp(sum_lnA(i)-sum_lnA(1))
		  sum_lnA_total = sum_lnA(i) + sum_lnA_total		  
		  !for charge conservation
		  sum_lnZ(i) = log(real(mt%Z(i))*m_term) + (mu_i(i)+abs(mt%BE(i)))/kT
 		  sum_lnZ(1) = log(real(mt%Z(1))*m_term1) + (mu_i(1)+abs(mt%BE(1)))/kT		
		  sum_lnZ(i) = exp(sum_lnZ(i)-sum_lnZ(1))
		  sum_lnZ_total = sum_lnZ(i) + sum_lnZ_total
		  !detailed balance
		  equ(i,1) = real(mt% Z(i))*(mu_n-mu_e+m_star)+real(mt% N(i))*mu_n-mu_i(i) !-abs(mt% BE(i)) 
		 enddo

          equ(1,1) = real(mt% Z(1))*(mu_n-mu_e+m_star)+real(mt% N(1))*mu_n-mu_i(1) !-abs(mt% BE(1))
         	
		  sum_lnA_final = sum_lnA(1) + log(1.0+sum_lnA_total)
    	  sum_lnZ_final = sum_lnZ(1) + log(1.0+sum_lnZ_total) 
		  		  
  		 !baryon and charge conservation 
         equ(mt% Ntable+1,1) = sum_lnZ_final - log(n_e) 
         equ(mt% Ntable+2,1) = sum_lnA_final - log(n_b-n_n) !- log(1.0 - Y_n/(1.0-chi))        

 		write(*,*) 'Y_e=', Y_e, 'mu_e=', mu_e, 'n_e=', n_e 
        write(*,*) 'Y_n=', Y_n, 'mu_n=', mu_n, 'n_n=', n_n
	    write(*,*) 'mu_i', mu_i(1), mu_i(5549)
	    write(*,*) 'sumZ=', sum_lnZ_final, 'log(n_e)=', log(n_e), 'equN_1=', equ(mt% Ntable+1,1)
	    write(*,*) 'sumA=', sum_lnA_final, 'log(n_b)=', log(n_b),  'equN_2=', equ(mt% Ntable+2,1)
     	write(*,*) '------------------------------'                   
                          
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
         real :: chi, rho      
		 !for loop over nuclei abundances
         real, parameter :: g=1.0d0
		 real :: m_star 
		 real :: m_nuc, m_nuc1
		 real :: m_term, m_term1		 
		 !for analytical jacobian
	     real, dimension(0:3), parameter :: cw0 = [ 1.2974, 15.0298, -15.2343, 7.4663 ]		
		 real :: dmudk_n, dmudk_e, dkdn_n, dkdn_e
		 real :: dmudne, dmudnn
		 !for equations in log space
		 real :: sum_lnZ(5549), sum_lnA(5549) 
		 real :: sum_lnZ_total, sum_lnZ_final 
		 real :: sum_lnA_total, sum_lnA_final
		 real :: logZ_exponent
		 real :: logA_exponent 
		 real :: ni_Zsum, ni_Asum
		 real :: n_i(5549)

		 
         ierr = 0        
	     chi = use_default_nuclear_size
         rho = (n_b*amu_n)*(mev_to_ergs/clight2)/(1.d-39) ! cgs

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
        n_e = ke**3/threepisquare               
        Y_e = n_e/n_b   
        mu_e = -mu_e
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
        n_n = 2.0*kn**3/threepisquare              
        Y_n = n_n*(1.-chi)/n_b   
		!Y_n = n_n/n_b
		end if
	
		if (mu_n == 0.) then
		kn=0.
		n_n = 0.
		Y_n = 0.
		end if 


		 sum_lnA = 0. ; sum_lnA_total = 0. ; sum_lnA_final = 0. 
		 sum_lnZ = 0. ; sum_lnZ_total = 0. ; sum_lnZ_final = 0. 


		 do i = 2, mt% Ntable
          !number density of isotopes
		  m_nuc = real(mt% Z(i))*mp_n + real(mt% N(i))*mn_n         
     	  m_term = g*(twopi*hbarc_n**2/(m_nuc*kT))**(-3.0/2.0)
		  !for baryon conservation
		  sum_lnA(i) = log(real(mt%A(i))*m_term) + (mu_i(i)+abs(mt%BE(i)))/kT
		  sum_lnA(i) = exp(sum_lnA(i)-sum_lnA(1))
		  sum_lnA_total = sum_lnA(i) + sum_lnA_total		  
		  !for charge conservation
		  sum_lnZ(i) = log(real(mt%Z(i))*m_term) + (mu_i(i)+abs(mt%BE(i)))/kT
		  sum_lnZ(i) = exp(sum_lnZ(i)-sum_lnZ(1))
		  sum_lnZ_total = sum_lnZ(i) + sum_lnZ_total
		 enddo

		n_i = 0. ; ni_Zsum = 0. ; ni_Asum = 0.
		dmudnn = 0. ; dmudne = 0.

		do i = 1, mt% Ntable
		  m_nuc = real(mt% Z(i))*mp_n + real(mt% N(i))*mn_n         
     	  m_term = g*(twopi*hbarc_n**2/(m_nuc*kT))**(-3.0/2.0)
     	  n_i(i) = m_term*exp((mu_i(i)+abs(mt%BE(i)))/kT)
     	  ni_Zsum = ni_Zsum + real(mt%Z(i))*n_i(i)
     	  ni_Asum = ni_Asum + real(mt%A(i))*n_i(i)
		enddo

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
		 
		 do i=1, mt% Ntable    
		    
            dkdn_e = (1.0/3.0)*(n_e*threepisquare)**(-2.0/3.0)*(threepisquare)			
			dmudk_e = (me_n/2.0)*(1.0+((ke*hbarc_n)/me_n)**2.0)**(-1.0/2.0) &
				& *2.0*(ke*hbarc_n/me_n)*hbarc_n/me_n
			dmudk_e = (hbarc_n**2.0*ke/me_n)*(1.0+((ke*hbarc_n)/me_n)**2.0)**(-1.0/2.0)
										
			dkdn_n = (1.0/3.0)*(n_n*threepisquare/2.0)**(-2.0/3.0)*(threepisquare/2.0)	
			!dmudk_n = (cw0(0) + 2.0*kn*(cw0(1) + kn*(3.0*cw0(2) + 4.0*kn*cw0(3)))) + onethird*  &
			!	& (cw0(0) + kn*(4.0*cw0(1) + kn*(9.0*cw0(2) + 16.0*kn*cw0(3))))				
			dmudk_n = (4.0/3.0)*cw0(0)+(10.0/3.0)*cw0(1)*kn &
				&	  + 6.0*cw0(2)*kn**2.0+(28.0/3.0)*cw0(3)*kn**3.0

			dmudne=dkdn_e*dmudk_e
			dmudnn=dkdn_n*dmudk_n

            if (n_n .eq. 0.) then
            dkdn_n = 0.
            dmudk_n = 0.
            end if
            if (n_e .eq. 0.) then
            dkdn_e = 0.
            dmudk_e = 0. 
            end if	
            	
		    m_nuc = real(mt% Z(i))*mp_n + real(mt% N(i))*mn_n         
     	    m_term = g*(twopi*hbarc_n**2/(m_nuc*kT))**(-3.0/2.0)  
 		    m_nuc1 = real(mt% Z(1))*mp_n + real(mt% N(1))*mn_n         
     	    m_term1 = g*(twopi*hbarc_n**2/(m_nuc1*kT))**(-3.0/2.0)      	              

            logZ_exponent = log(real(mt%Z(i))*m_term) + (mu_i(i)+mt%BE(i))/kT &
            			&  - log(real(mt%Z(1))*m_term1) + (mu_i(1)+mt%BE(1))/kT
            logA_exponent = log(real(mt%A(i))*m_term) + (mu_i(i)+mt%BE(i))/kT &
            			&  - log(real(mt%A(1))*m_term1) + (mu_i(1)+mt%BE(1))/kT            			
            
 			! n=1 to N terms of last two rows
      		!A(mt% Ntable+1, i) = (1.0/kT)*(1.0+sum_lnZ_total)**(-1.0) &
      		!		& *exp(logZ_exponent)
     		!with summation
       		A(mt% Ntable+1, i) = (1.0/kT)-(1.0/kT)*sum_lnZ_total &
       				& /(1.0+sum_lnZ_total)	
     		A(mt% Ntable+2, i) = (1.0/kT)*(1.0+sum_lnA_total)**(-1.0) &
     				& *exp(logA_exponent)
  
  		    !last two rows of jacobian, derivatives wrt the conservation equations 
     		! n=0 term
     		!A(mt% Ntable+1, 1) = (1.0/kT)-(1.0/kT)*sum_lnZ_total*(1.0+sum_lnZ_total)**(-1.0)
     		!A(mt% Ntable+2, 1) = (1.0/kT)-(1.0/kT)*sum_lnA_total*(1.0+sum_lnA_total)**(-1.0)

		    !last two columns 
			A(i, mt% Ntable+1) = -real(mt% Z(i)) !*dmudk_e*dkdn_e*n_b			 ! MeV		    
			A(i, mt% Ntable+2) = real(mt% A(i)) !*dmudk_n*dkdn_n*n_b/(1.0-chi) ! MeV
     
   			enddo  	
     		!last two rows
     		!A(mt% Ntable+1, i) = real(mt% Z(i))*n_i(i)/(kT*ni_Zsum)
     		!A(mt% Ntable+2, i) = real(mt% A(i))*n_i(i)/(kT*ni_Asum)
     			 	    		       	
			A(mt% Ntable+1, mt% Ntable+1) = 0. !-1.0/Y_e			
!			A(mt% Ntable+1, mt% Ntable+2) = (1.0/ni_Zsum)*real(mt%Z(i))*n_i(i) &
!				& *real(mt%A(i))/kT*dmudnn*n_b/(1.0-chi)		
				
			A(mt% Ntable+2, mt% Ntable+2) = 0. !1.0 /((1.0-chi)-Y_n)  
!			A(mt% Ntable+2, mt% Ntable+1) = (1.0/ni_Asum)*real(mt%A(i))*n_i(i) &
!				& *(-real(mt%Z(i)))/kT*dmudne*n_b	
						
			!try
			!A(i, mt% Ntable+1) = 0.
			!A(i, mt% Ntable+2) = 0.
			A(mt% Ntable+1, mt% Ntable+2) = 0.
			A(mt% Ntable+2, mt% Ntable+1) = 0. 			
						
			!try	
			!A(mt% Ntable+1, i) = 1./kT
			!A(mt% Ntable+2, i) = 1./kT
			!A(i, mt% Ntable+2) = real(mt% A(i))*dmudk_n*dkdn_n*n_b/(1.0-chi) ! MeV
			!A(mt% Ntable+2, mt% Ntable+2) = 1.0/((1.0-chi)-Y_n)

		        
         
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
		real, intent(in) :: k	! (fm**-3)
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

     function ke_solve(x)
      real, intent(in) :: x
      real :: ke_solve    
      ke_solve = electron_chemical_potential(x) - mu_e  
     end function ke_solve  
      
    end module qse_solver
