module short_table
	! storage for the mass table
	type short_table_type
		integer :: Ntable	! number of table entries
		integer :: Zmin, Zmax	! minimum, maximum Z
		integer, dimension(:), pointer :: Z, A, N
		real, dimension(:), pointer :: BE	! binding energies (MeV)
		real, dimension(:), pointer :: Sn	! neutron separation energy (MeV)
		real, dimension(:), pointer :: S2n	! dineutron separation energy (MeV)
		real, dimension(:), pointer :: Sp	! proton separation energy (MeV)
		real, dimension(:), pointer :: S2p  ! diproton separation energy (MeV)
		real, dimension(:), pointer :: Ec	! electron capture threshold (MeV)
		real, dimension(:), pointer :: beta ! beta- threshold (MeV)
		real, dimension(:), pointer :: VN	! volume occupied by the nucleus (fm**3)
		! bookmarks
		integer :: Nelements	! = Zmax-Zmin+1
		!index of first nuclide of a given element
		integer, dimension(:), pointer :: Zstart	! (Nelements) 
		! minimum, maximum neutron no. for each element
		integer, dimension(:), pointer :: Nmin, Nmax	! (Nelements)
	end type short_table_type
	logical, save :: short_table_is_loaded = .FALSE.

	type(short_table_type), target :: winvn_short_table

contains

	subroutine short_table_shutdown()
		type(short_table_type), pointer :: mt
		if (.not.short_table_is_loaded) return
		mt => winvn_short_table
		deallocate(mt% Zstart)
		deallocate(mt% Nmin)
		deallocate(mt% Nmax)
		deallocate(mt% Z)
		deallocate(mt% A)
		deallocate(mt% N)
		deallocate(mt% BE)
		deallocate(mt% Sn)
		deallocate(mt% S2n)
		deallocate(mt% Sp)
		deallocate(mt% S2p)
		deallocate(mt% Ec)
		deallocate(mt% beta)
		deallocate(mt% VN)
		mt% Ntable = 0
		mt% Zmin = -1
		mt% Zmax = -1
		mt% Nelements = 0
	end subroutine short_table_shutdown

end module short_table

