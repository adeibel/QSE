module pressure_table
	! storage for the mass table
	type pressure_table_type
		integer :: Ntable	! number of table entries
		!integer :: Zmin, Zmax	! minimum, maximum Z
		integer :: Pmin, Pmax
		real, dimension(:), pointer :: P    ! pressure coordinate
		real, dimension(:), pointer :: mu_e ! electron chemical potential
		real, dimension(:), pointer :: mu_n ! neutron chemical potential
		! bookmarks
		integer :: Nelements	! = Zmax-Zmin+1
		!index of first nuclide of a given element
		integer, dimension(:), pointer :: Pstart	! (Nelements) 
		!real, dimension(:), pointer :: Pstart		
		! minimum, maximum neutron no. for each element
		!integer, dimension(:), pointer :: mu_n_min, mu_n_max	! (Nelements)
		real, dimension(:), pointer :: mu_n_min, mu_n_max
		
	end type pressure_table_type
	logical, save :: pressure_table_is_loaded = .FALSE.

	type(pressure_table_type), target :: winvn_pressure_table

 contains

	subroutine pressure_table_shutdown()
		type(pressure_table_type), pointer :: pt
		if (.not.pressure_table_is_loaded) return
		pt => winvn_pressure_table
		deallocate(pt% Pstart)
		deallocate(pt% P)
		deallocate(pt% mu_n_min)
		deallocate(pt% mu_n_max)
		deallocate(pt% mu_e)
		deallocate(pt% mu_n)
		pt% Ntable = 0
		pt% Pmin = -1
		pt% Pmax = -1
		pt% Nelements = 0
	end subroutine pressure_table_shutdown

	subroutine load_pressure_table(datadir,datafile,ierr)
		use alert_lib
		use utils_lib
		character(len=*), intent(in) :: datadir,datafile
		integer, intent(out) :: ierr
		character(len=*), parameter :: default_dir = 'crust_eos'
		character(len=256) :: filename
		type(pressure_table_type), pointer :: pt
		integer :: i, Ntab, Nele, iounit
		integer :: currentP, id_el

		if (pressure_table_is_loaded) then
			ierr = 1
			call alert(1,'pressure_table_support: load_pressure_table:: table is already loaded')
			return
		end if
		
		iounit = alloc_iounit(ierr)
		if (ierr /= 0) then
			call alert(ierr,'pressure_table_support: load_pressure_table:: unable to alloc iounit')
			return
		end if
		filename = trim(datadir)//'/'//default_dir//'/'//trim(datafile)
		open(unit=iounit, file=trim(filename), iostat=ierr, action='read')
		if (ierr /= 0) then
			call alert(ierr,'pressure_table_support: load_pressure_table:: unable to open '//trim(filename))
			return
		end if

		! read through the file to get number of entries
		Ntab = 0
		read(iounit,*,iostat=ierr)
		if (ierr /= 0) then
			call alert(ierr,'pressure_table_support: load_pressure_table:: unable to read any lines')
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
		pt => winvn_pressure_table
	
		pt% Ntable = Ntab
		! allocate the tables
		allocate(pt% P(Ntab), pt% mu_e(Ntab), pt% mu_n(Ntab))

		! now read in the table, skipping first line
		read(iounit,*,iostat=ierr)
		if (ierr /= 0) then
			call alert(ierr,'pressure_table_support: load_pressure_table:: unable to read any lines')
			return
		end if
		
		do i = 1, pt% Ntable
			read(iounit,*,iostat=ierr) pt% P(i), pt% mu_e(i), pt% mu_n(i)
			if (ierr /=0) then
			   call alert(ierr,'pressure_table_support: load_pressure_table:: unable to read lines')
			   exit
		   end if
		end do
		close(iounit)
		call free_iounit(iounit)

!fix following

		! set up the bookmarks; table is assumed to be sorted!
		pt% Pmin = 1	   !minval(pt% P)
		pt% Pmax = pt% Ntable  !maxval(pt% P)
		
		! allocate the bookmarks
		Nele = pt% Pmax- pt% Pmin+1
		pt% Nelements = Nele
		allocate(pt% Pstart(Nele), pt% mu_n_min(Nele), pt% mu_n_max(Nele))
		
		!currentP = -1
	    currentP = 0.
		id_el = 0
		do i = 1, pt% Ntable
			if (pt% P(i) /= currentP) then	! found a new pressure
				id_el = id_el + 1
				pt% Pstart(id_el) = i   ! Pstart can be an integer?
				pt% mu_n_min(id_el) = pt% mu_n(i)
				if (id_el > 1) then
					pt% mu_n_max(id_el-1) = pt% mu_n(i-1)
				end if
				currentP = pt% P(i)
			end if
		end do
		pt% mu_n_max(id_el) = pt% mu_n(pt% Ntable)

		if (ierr == 0) pressure_table_is_loaded = .TRUE.
	end subroutine load_pressure_table

	function pressure_table_index(P, mu_n, ierr) result(id)
		real, intent(in) :: P, mu_n
		integer, intent(out) :: ierr
		integer :: id
		integer :: Pindex
		type(pressure_table_type), pointer :: pt
	
		pt => winvn_pressure_table
		ierr = 0
		id = -1
		if (P < pt% Pmin .or.  P > pt% Pmax) then
			ierr = -1
			return
		end if
		Pindex = P - pt% Pmin + 1
		if (mu_n < pt% mu_n_min(Pindex) .or. mu_n > pt% mu_n_max(Pindex)) then
			ierr = -2
			return
		end if
	
		id = pt% Pstart(Pindex) + mu_n - pt% mu_n_min(Pindex)
	end function pressure_table_index

end module pressure_table