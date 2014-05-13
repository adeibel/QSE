module stable_table
	! storage for the mass table
	type stable_table_type
		integer :: Ntable	! number of table entries
		!integer :: Zmin, Zmax	! minimum, maximum Z
		integer :: Pmin, Pmax
		integer, dimension(:), pointer :: Z, A, N
		real, dimension(:), pointer :: BE	! binding energies (MeV)
		real, dimension(:), pointer :: Sn	! neutron separation energy (MeV)
		real, dimension(:), pointer :: S2n	! dineutron separation energy (MeV)
		real, dimension(:), pointer :: Sp	! proton separation energy (MeV)
		real, dimension(:), pointer :: S2p  ! diproton separation energy (MeV)
		real, dimension(:), pointer :: Ec	! electron capture threshold (MeV)
		real, dimension(:), pointer :: beta ! beta- threshold (MeV)
		real, dimension(:), pointer :: VN	! volume occupied by the nucleus (fm**3)
		real, dimension(:), pointer :: P    ! pressure coordinate
		real, dimension(:), pointer :: mu_e ! electron chemical potential
		real, dimension(:), pointer :: mu_n ! neutron chemical potential
		real, dimension(:), pointer :: eden ! energy density
		real, dimension(:), pointer :: rho  ! mass denisty
		real, dimension(:), pointer :: BEA  ! binding energy per nucleon
		! bookmarks
		integer :: Nelements	! = Zmax-Zmin+1
		!index of first nuclide of a given element
		integer, dimension(:), pointer :: Pstart	! (Nelements) 
		!real, dimension(:), pointer :: Pstart
		
		! minimum, maximum neutron no. for each element
		!integer, dimension(:), pointer :: mu_n_min, mu_n_max	! (Nelements)
		real, dimension(:), pointer :: mu_n_min, mu_n_max
		
	end type stable_table_type
	logical, save :: stable_table_is_loaded = .FALSE.

	type(stable_table_type), target :: winvn_stable_table

contains

	subroutine stable_table_shutdown()
		type(stable_table_type), pointer :: st
		if (.not.stable_table_is_loaded) return
		st => winvn_stable_table
		deallocate(st% Pstart)
		deallocate(st% mu_n_min)
		deallocate(st% mu_n_max)
		deallocate(st% Z)
		deallocate(st% A)
		deallocate(st% N)
		deallocate(st% BE)
		deallocate(st% Sn)
		deallocate(st% S2n)
		deallocate(st% Sp)
		deallocate(st% S2p)
		deallocate(st% Ec)
		deallocate(st% beta)
		deallocate(st% VN)
		deallocate(st% P)
		deallocate(st% mu_e)
		deallocate(st% mu_n)
		deallocate(st% BEA)
		deallocate(st% eden)
		deallocate(st% rho)
		st% Ntable = 0
		st% Pmin = -1
		st% Pmax = -1
		st% Nelements = 0
	end subroutine stable_table_shutdown

	subroutine load_stable_table(datadir,datafile,ierr)
		use alert_lib
		use utils_lib
		character(len=*), intent(in) :: datadir,datafile
		integer, intent(out) :: ierr
		character(len=*), parameter :: default_dir = 'crust_eos'
		character(len=256) :: filename
		type(stable_table_type), pointer :: st
		integer :: i, Ntab, Nele, iounit
		integer :: currentP, id_el

		if (stable_table_is_loaded) then
			ierr = 1
			call alert(1,'stable_table_support: load_stable_table:: table is already loaded')
			return
		end if
		
		iounit = alloc_iounit(ierr)
		if (ierr /= 0) then
			call alert(ierr,'stable_table_support: load_stable_table:: unable to alloc iounit')
			return
		end if
		filename = trim(datadir)//'/'//default_dir//'/'//trim(datafile)
		open(unit=iounit, file=trim(filename), iostat=ierr, action='read')
		if (ierr /= 0) then
			call alert(ierr,'stable_table_support: load_stable_table:: unable to open '//trim(filename))
			return
		end if

		! read through the file to get number of entries
		Ntab = 0
		read(iounit,*,iostat=ierr)
		if (ierr /= 0) then
			call alert(ierr,'stable_table_support: load_stable_table:: unable to read any lines')
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
		st => winvn_stable_table
	
		st% Ntable = Ntab
		! allocate the tables
		allocate(st% P(Ntab), st% mu_e(Ntab), st% mu_n(Ntab), st% Z(Ntab), st% A(Ntab), &
				st% N(Ntab), st% BE(Ntab), st% BEA(Ntab), st% Sn(Ntab), st% Sp(Ntab), &
				st% eden(Ntab), st% rho(Ntab))

		! now read in the table, skipping first line
		read(iounit,*,iostat=ierr)
		if (ierr /= 0) then
			call alert(ierr,'stable_table_support: load_stable_table:: unable to read any lines')
			return
		end if
		
		do i = 1, st% Ntable
			read(iounit,*,iostat=ierr) st% P(i), st% mu_e(i), st% mu_n(i), &
				& st% Z(i), st% A(i), st% N(i), st% BE(i), st% BEA(i), &
				& st% Sn(i), st% Sp(i), st% eden(i), st% rho(i)
			if (ierr /=0) then
			   call alert(ierr,'stable_table_support: load_stable_table:: unable to read lines')
			   exit
		   end if
		end do
		close(iounit)
		call free_iounit(iounit)

!fix following

		! set up the bookmarks; table is assumed to be sorted!
		st% Pmin = 1	   !minval(st% P)
		st% Pmax = st% Ntable  !maxval(st% P)
		
		! allocate the bookmarks
		Nele = st% Pmax- st% Pmin+1
		st% Nelements = Nele
		allocate(st% Pstart(Nele), st% mu_n_min(Nele), st% mu_n_max(Nele))
		
		!currentP = -1
	    currentP = 0.
		id_el = 0
		do i = 1, st% Ntable
			if (st% P(i) /= currentP) then	! found a new pressure
				id_el = id_el + 1
				st% Pstart(id_el) = i   ! Pstart can be an integer?
				st% mu_n_min(id_el) = st% mu_n(i)
				if (id_el > 1) then
					st% mu_n_max(id_el-1) = st% mu_n(i-1)
				end if
				currentP = st% P(i)
			end if
		end do
		st% mu_n_max(id_el) = st% mu_n(st% Ntable)

		if (ierr == 0) stable_table_is_loaded = .TRUE.
	end subroutine load_stable_table

	function stable_table_index(P, mu_n, ierr) result(id)
		real, intent(in) :: P, mu_n
		integer, intent(out) :: ierr
		integer :: id
		integer :: Pindex
		type(stable_table_type), pointer :: st
	
		st => winvn_stable_table
		ierr = 0
		id = -1
		if (P < st% Pmin .or.  P > st% Pmax) then
			ierr = -1
			return
		end if
		Pindex = P - st% Pmin + 1
		if (mu_n < st% mu_n_min(Pindex) .or. mu_n > st% mu_n_max(Pindex)) then
			ierr = -2
			return
		end if
	
		id = st% Pstart(Pindex) + mu_n - st% mu_n_min(Pindex)
	end function stable_table_index

end module stable_table