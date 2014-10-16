	module eos_lib
	   use iso_fortran_env, only : error_unit
	   use utils_lib
	   use alert_lib
	   use mass_table
	   use phys_constants
	   use rootfind
	   use MB77, only: neutron_chemical_potential
   
	 contains

	! all number densities and energy density are in nuclear units
	 subroutine get_energy_density(Z, A, mu_e, mu_n, BE, kT, energy_density)
	    real, intent(in) :: Z, A
		real, intent(in) :: mu_e, mu_n, BE, kT
		real :: n_n, n_p, n_i, n_e
		real :: k_fe, k_fn, k_fp
		real :: rest_mass
		real :: energy_density
		real :: x1, x2, xacc
		real :: r
		real,dimension(:),pointer :: hist
		real, parameter :: g=2.0
		integer :: ierr 
		!electron number density
		k_fe = sqrt(mu_e**2-me_n**2)
		n_e = k_fe**3/(3*pi**2) 
		n_i = n_e/Z
		!neutron number density rootfind
		x1=0.0
		x2=10.0
		xacc=1.d-15
		n_n = root_bisection(neutron_chem,x1,x2,xacc,ierr,hist) 
		if (ierr /= 0) then
		write(*,*) 'Error in bisection for number density of neutrons'
		stop
		endif		 		
		rest_mass = (Z*mp_n+(A-Z)*mn_n+Z*me_n)*n_i+mn_n*n_n
		energy_density = BE*n_i+abs(eden_drip(n_n)-mn_n*n_n) &
				+(abs(eden_el(n_e)-Z*n_i*me_n))+rest_mass
	 contains 
 
     ! neutron chemical potential from mb77
	 function neutron_chem(n)
	  real, intent(in) :: n
	  real :: neutron_chem
	  neutron_chem = neutron_chemical_potential(n)-abs(mu_n)   
	 end function
 
	 ! electron energy density
	 ! Roca-Maza (2008)
	 function eden_el(n_e)
		real ::  k_f, eden_el, x, y
		real ::  n_e
		k_f = (3.0*pi**2*n_e)**(1./3.)
		x = k_f/me_n
		y = sqrt(1.0+x**2)
		eden_el = ((me_n**4)/(8.0*pi**2))*(x*y*(x**2+y**2)-log(x+y)) 
	 end function eden_el

	 !energy density of quasi-free neutrons
	 function eden_drip(n_n)
		real :: eden_drip
		real :: h1, h2, h3, h4, h5, h6
		real :: t0, t1, t2, t3, eps
		real :: x0, x1, x2, x3
		real :: eta1, eta2, m_nstar, m_pstar
		real :: k_n, k_p, n 
		real :: n_n, n_p	
		! defaults
		n_p = 0.0
		! SLy4 coefficients
		t0 = -2488.913; t1 = 486.818; t2 = -546.395; t3 = 13777. 
		x0 = .834; x1 = -.3438; x2 = -1.; x3 = 1.354
		eps = 1./6. 

		n = n_n+n_p
		! effective masses
		eta1 = .25*(t1*(1.+x1/2.)+t2*(1.+x2/2.))
		eta2 = .25*(t2*(.5+x2)-t1*(.5+x1))
		m_nstar = mn_n/(1.+(5.06773d-3)*2.*(n*eta1+n_n*eta2)*mn_n)
		m_pstar = mp_n/(1.+(5.06773d-3)*2.*(n*eta1+n_p*eta2)*mp_n)
		! fermi momenta
		k_n = (3.*(pi**2)*n_n)**(1./3.) 
		k_p = (3.*(pi**2)*n_p)**(1./3.)
		! six hamiltonian terms
		h1 = (k_n**5)/(10.*m_nstar*pi**2)+n_n*mn_n
		h2 = (k_p**5)/(10.*m_pstar*pi**2)+n_p*mp_n
		h3 = (5.06773d-3)*(t0/2.)*(1.+x0/2.)*n*n
		h4 = -(5.06773d-3)*(t0/2.)*(.5+x0)*(n_n*n_n+n_p*n_p)
		h5 = (5.06773d-3)*(t3/12.)*(1.+x3/2.)*n*n*(n**eps)
		h6 = -(5.06773d-3)*(t3/12)*(.5+x3)*(n_n*n_n+n_p*n_p)*(n**eps) 
		! homogeneous bulk matter hamiltonian
		eden_drip = h1+h2+h3+h4+h5+h6
	 end function eden_drip

	 end subroutine get_energy_density

end module eos_lib