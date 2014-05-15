module check_stability
   use iso_fortran_env, only: error_unit, oid=>output_unit
   use alert_lib
   use phys_constants
   use mass_table 
   use minimum_nucleus


contains

subroutine stable_nuclei(mu_e, mu_n, kT, mass_table_name, Z_best, A_best, pressure, &
		stable_id, mu_id)
   real, intent(in) :: mu_e, mu_n, kT   !MeV
   real, intent(in) :: pressure			!MeV fm**-3
   character(len=*), intent(in) :: mass_table_name
   character(len=*), parameter :: output_file = 'stable.data'
   character(len=256) :: arg
   real :: energy_density, eden_save(10000)
   real :: mass_density, mass_save(10000)
   real :: mu_p
   real, parameter :: del_m = mn_n-mp_n-me_n
   integer :: Z, A, ierr, i, stable_id, mu_id
   integer :: Z_best(1), A_best(1)
   integer :: best_nucleus_id(1)
   type(mass_table_type), pointer :: mt
   logical, dimension(:), allocatable :: beta_stable, neutron_stable, proton_stable
   logical, dimension(:), allocatable :: totally_stable
   logical, save :: mass_table_is_loaded = .FALSE.

	mu_p = mu_n-mu_e

   if (mass_table_is_loaded .eqv. .FALSE.) then
   call load_mass_table('../../../data/',trim(mass_table_name),ierr)
   mass_table_is_loaded = .TRUE.
   if (ierr /= 0) then
      write(error_unit,'(a)') trim(alert_message)
      stop
   end if
   end if
   mt => winvn_mass_table
   
   allocate(beta_stable(mt% Ntable), neutron_stable(mt% Ntable), proton_stable(mt% Ntable))
   allocate(totally_stable(mt% Ntable))
   
   ! for equil_nuclei.data file 
   eden_save = 1.d20
   do i = 1, mt% Ntable
    beta_stable(i) = nucleus_is_beta_stable(i)
    neutron_stable(i) = nucleus_is_neutron_stable(i)
    !proton_stable(i) = nucleus_is_proton_stable(i)
    !if (beta_stable(i) .and. neutron_stable(i) .and. proton_stable(i) .eqv. .TRUE.) then
	if (beta_stable(i) .and. neutron_stable(i) .eqv. .TRUE.) then
      call energy_density_solve(mt% Z(i), mt% A(i), mu_e, mu_n, mt% BE(i), kT, energy_density, &
      			mass_density)
      eden_save(i) = energy_density
      mass_save(i) = mass_density
    end if
   end do
   

  ! find equilibrium nucleus in energy density array
  if (count(beta_stable .and. neutron_stable) > 0) then
  best_nucleus_id = minloc(eden_save)
  Z_best = mt% Z(best_nucleus_id)
  A_best = mt% A(best_nucleus_id)
  end if
  
  ! for stable_nuclei.data file
  if (mu_n .gt. 0.0) then
   do i = 1, mt% Ntable
 !  	  if (mt% Z(i) == 51 .and. mt% A(i) == 172) then
!   if (beta_stable(i) .and. neutron_stable(i) .and. proton_stable(i) .eqv. .TRUE.) then
	if (beta_stable(i) .and. neutron_stable(i) .eqv. .TRUE.) then
       write(stable_id,'(3es13.5,3i5,4f12.4,2x,es12.4,2x,es12.4)') pressure, mu_e, mu_n, mt% Z(i), mt% A(i), mt% N(i),&
         mt% BE(i), mt% BE(i)/mt% A(i),mt% Sn(i), mt% Sp(i), eden_save(i), mass_save(i)
!      end if
	totally_stable(i) = .TRUE.
      end if 
   end do    
  end if
  
  write(*,*) mu_e, mu_n, count(totally_stable)
  
return

 contains

   function nucleus_is_beta_stable(id) 
      integer, intent(in) :: id
      logical :: nucleus_is_beta_stable
      integer :: Z, A, Ar, Zr, Nr, idr, ierr
      real :: B, Br
      
      A = mt% A(id)
      Z = mt% Z(id)
      B = mt% BE(id)
      
      nucleus_is_beta_stable = .FALSE.

      ! check stablity
      ! electron capture
      Zr = Z-1
      Ar = A
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if (-mu_e + del_m + B - Br < 0) return
      
      ! electron emission
      Zr = Z+1
      Ar = A
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if (mu_e - del_m + B - Br < 0) return
           
      ! electron capture followed by neutron emission
      Zr = Z-1
      Ar = A-1
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if (mu_n - mu_e + del_m + B - Br < 0) return
      
      ! electron capture followed by dineutron emission
      Zr = Z-1
      Ar = A-2
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if (2*mu_n - mu_e + del_m + B - Br < 0) return
      
      ! if we got to here, assume the nucleus is stable
      nucleus_is_beta_stable = .TRUE.
     
   end function nucleus_is_beta_stable

   function nucleus_is_neutron_stable(id)
      integer, intent(in) :: id
      logical :: nucleus_is_neutron_stable
      integer :: Z, A, Ar, Zr, Nr, idr, ierr
      real :: B, Br
      
      A = mt% A(id)
      Z = mt% Z(id)
      B = mt% BE(id)
      
      nucleus_is_neutron_stable = .FALSE.

      ! check stablity
      ! neutron capture
      Zr = Z
      Ar = A+1
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if ( -mu_n + B - Br < 0) return

      ! neutron emission
      Zr = Z
      Ar = A-1
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if (mu_n + B - Br < 0) return

      ! dineutron capture
      Zr = Z
      Ar = A+2
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if (-2*mu_n + B - Br < 0) return

      
      ! dineutron emission
      Zr = Z
      Ar = A-2
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if (2*mu_n + B - Br < 0) return


      ! neutron capture followed by electron emission
      Zr = Z+1
      Ar = A+1
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if (mu_e - mu_n - del_m + B - Br < 0) return
      
      ! dineutron capture followed by electron emission
      Zr = Z+1
      Ar = A+2
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if (mu_e - 2.0*mu_n - del_m + B - Br < 0) return
      
      
      ! if we got to here, assume the nucleus is stable
      nucleus_is_neutron_stable = .TRUE.

   end function nucleus_is_neutron_stable
 
    function nucleus_is_proton_stable(id)
      integer, intent(in) :: id
      logical :: nucleus_is_proton_stable
      integer :: Z, A, Ar, Zr, Nr, idr, ierr
      real :: B, Br
      
      A = mt% A(id)
      Z = mt% Z(id)
      B = mt% BE(id)
      
      nucleus_is_proton_stable = .FALSE.

      ! check stablity
      ! proton capture
      Zr = Z+1
      Ar = A+1
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if ( -mu_p + B - Br < 0) return
      
       ! proton emission
      Zr = Z-1
      Ar = A-1
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if (mu_p + B - Br < 0) return

      
      ! diproton capture
      Zr = Z+2
      Ar = A+2
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if (-2*mu_p + B - Br < 0) return

      
      ! diproton emission
      Zr = Z-2
      Ar = A-2
      Nr = Ar-Zr
      idr = mass_table_index(Zr,Nr,ierr)
      if (ierr /= 0) return
      Br = mt% BE(idr)
      if (2*mu_p + B - Br < 0) return


      ! proton capture followed by electron emission?
!      Zr = Z+1
!      Ar = A+1
!      Nr = Ar-Zr
!      idr = mass_table_index(Zr,Nr,ierr)
!      if (ierr /= 0) return
!      Br = mt% BE(idr)
!      if (mu_e - mu_n - del_m + B - Br < 0) return

		nucleus_is_proton_stable = .TRUE.
      
   end function nucleus_is_proton_stable
 
    
 end subroutine stable_nuclei 

end module check_stability