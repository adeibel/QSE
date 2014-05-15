module mass_table
	! storage for the mass table
	type mass_table_type
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
	end type mass_table_type
	logical, save :: mass_table_is_loaded = .FALSE.

	type(mass_table_type), target :: winvn_mass_table

contains

	subroutine mass_table_shutdown()
		type(mass_table_type), pointer :: mt
		if (.not.mass_table_is_loaded) return
		mt => winvn_mass_table
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
	end subroutine mass_table_shutdown

	subroutine load_mass_table(datadir,datafile,ierr)
		use alert_lib
		use utils_lib
		character(len=*), intent(in) :: datadir,datafile
		integer, intent(out) :: ierr
		character(len=*), parameter :: default_dir = 'crust_eos'
		character(len=256) :: filename
		type(mass_table_type), pointer :: mt
		integer :: i, Ntab, Nele, iounit
		integer :: currentZ, id_el

		if (mass_table_is_loaded) then
			ierr = 1
			call alert(1,'mass_table_support: load_mass_table:: table is already loaded')
			return
		end if
		
		iounit = alloc_iounit(ierr)
		if (ierr /= 0) then
			call alert(ierr,'mass_table_support: load_mass_table:: unable to alloc iounit')
			return
		end if
		filename = trim(datadir)//'/'//default_dir//'/'//trim(datafile)
		open(unit=iounit, file=trim(filename), iostat=ierr, action='read')
		if (ierr /= 0) then
			call alert(ierr,'mass_table_support: load_mass_table:: unable to open '//trim(filename))
			return
		end if

		! read through the file to get number of entries
		Ntab = 0
		read(iounit,*,iostat=ierr)
		if (ierr /= 0) then
			call alert(ierr,'mass_table_support: load_mass_table:: unable to read any lines')
			return
		end if
		do
			read(iounit,*,iostat=ierr)
			if (ierr /= 0) exit
			Ntab = Ntab + 1
		end do
		
		ierr = 0
		rewind(iounit,iostat=ierr)
		if (ierr /= 0) then
			! try closing and reopening
			close(iounit)
			open(unit=iounit,file=trim(filename),iostat=ierr,action='read')
			if (ierr /= 0) then
				call alert(ierr,'unable to rewind or reopen '//trim(filename))
				return
			end if
		end if
		
		ierr = 0
		mt => winvn_mass_table
		mt% Ntable = Ntab
		! allocate the tables
		allocate(mt% Z(Ntab), mt% N(Ntab), mt% A(Ntab), mt% BE(Ntab), mt% Sn(Ntab), mt% S2n(Ntab), &
			& mt% Sp(Ntab), mt% S2p(Ntab), mt% Ec(Ntab), mt% beta(Ntab), mt% VN(Ntab))

		! now read in the table, skipping first line
		read(iounit,*,iostat=ierr)
		if (ierr /= 0) then
			call alert(ierr,'mass_table_support: load_mass_table:: unable to read any lines')
			return
		end if
		
		do i = 1, mt% Ntable
			read(iounit,*,iostat=ierr) mt% Z(i), mt% N(i), mt% A(i), mt% BE(i), mt% Sn(i), mt% S2n(i),  &
					& mt% Sp(i), mt % S2p(i), mt% EC(i), mt% beta(i), mt% VN(i)
			if (ierr /=0) then
			   call alert(ierr,'mass_table_support: load_mass_table:: unable to read lines')
			   exit
		   end if
		end do
		close(iounit)
		call free_iounit(iounit)

		! set up the bookmarks; table is assumed to be sorted!
		mt% Zmin = minval(mt% Z)
		mt% Zmax = maxval(mt% Z)
		
		! allocate the bookmarks
		Nele = mt% Zmax-mt% Zmin+1
		mt% Nelements = Nele
		allocate(mt% Zstart(Nele), mt% Nmin(Nele), mt% Nmax(Nele))
		
		currentZ = -1
		id_el = 0
		do i = 1, mt% Ntable
			if (mt% Z(i) /= currentZ) then	! found a new element
				id_el = id_el + 1
				mt% Zstart(id_el) = i
				mt% Nmin(id_el) = mt% N(i)
				if (id_el > 1) then
					mt% Nmax(id_el-1) = mt% N(i-1)
				end if
				currentZ = mt% Z(i)
			end if
		end do
		mt% Nmax(id_el) = mt% N(mt% Ntable)

		if (ierr == 0) mass_table_is_loaded = .TRUE.
	end subroutine load_mass_table

	function mass_table_index(Z, N, ierr) result(id)
		integer, intent(in) :: Z, N
		integer, intent(out) :: ierr
		integer :: id
		integer :: Zindex
		type(mass_table_type), pointer :: mt
	
		mt => winvn_mass_table
		ierr = 0
		id = -1
		if (Z < mt% Zmin .or. Z > mt% Zmax) then
			ierr = -1
			return
		end if
		Zindex = Z - mt% Zmin + 1
		if (N < mt% Nmin(Zindex) .or. N > mt% Nmax(Zindex)) then
			ierr = -2
			return
		end if
	
		id = mt% Zstart(Zindex) + N - mt% Nmin(Zindex)
	end function mass_table_index

end module mass_table