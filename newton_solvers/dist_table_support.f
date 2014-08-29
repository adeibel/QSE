module dist_table
	! storage for the mass table
	type dist_table_type
		integer :: Ntable	! number of table entries
		!integer :: Zmin, Zmax	! minimum, maximum Z
		integer :: Pmin, Pmax
		real, dimension(:), pointer :: P    ! dist coordinate
		real, dimension(:), pointer :: mu_e ! electron chemical potential
		real, dimension(:), pointer :: mu_n ! neutron chemical potential
		integer, dimension(:), pointer :: Z, A
		real, dimension(:), pointer :: abun
		! bookmarks
		integer :: Nelements	! = Zmax-Zmin+1
		!index of first nuclide of a given element
		integer, dimension(:), pointer :: Pstart	! (Nelements) 
		!real, dimension(:), pointer :: Pstart		
		! minimum, maximum neutron no. for each element
		!integer, dimension(:), pointer :: mu_n_min, mu_n_max	! (Nelements)
		real, dimension(:), pointer :: mu_n_min, mu_n_max	
	end type dist_table_type
	logical, save :: dist_table_is_loaded = .FALSE.

	type(dist_table_type), target :: winvn_dist_table

 contains

	subroutine dist_table_shutdown()
		type(dist_table_type), pointer :: dt
		if (.not.dist_table_is_loaded) return
		dt => winvn_dist_table
		deallocate(dt% Pstart)
		deallocate(dt% P)
		deallocate(dt% mu_n_min)
		deallocate(dt% mu_n_max)
		deallocate(dt% mu_e)
		deallocate(dt% mu_n)
		deallocate(dt% Z)
		deallocate(dt% A)
		deallocate(dt% abun)
		dt% Ntable = 0
		dt% Pmin = -1
		dt% Pmax = -1
		dt% Nelements = 0
	end subroutine dist_table_shutdown

	subroutine load_dist_table(datadir,datafile,ierr)
		use alert_lib
		use utils_lib
		character(len=*), intent(in) :: datadir,datafile
		integer, intent(out) :: ierr
		character(len=*), parameter :: default_dir = 'crust_eos'
		character(len=256) :: filename
		type(dist_table_type), pointer :: dt
		integer :: i, Ntab, Nele, iounit
		integer :: currentP, id_el

		if (dist_table_is_loaded) then
			ierr = 1
			call alert(1,'dist_table_support: load_dist_table:: table is already loaded')
			return
		end if
		
		iounit = alloc_iounit(ierr)
		if (ierr /= 0) then
			call alert(ierr,'dist_table_support: load_dist_table:: unable to alloc iounit')
			return
		end if
		filename = trim(datadir)//'/'//default_dir//'/'//trim(datafile)
		open(unit=iounit, file=trim(filename), iostat=ierr, action='read')
		if (ierr /= 0) then
			call alert(ierr,'dist_table_support: load_dist_table:: unable to open '//trim(filename))
			return
		end if

		! read through the file to get number of entries
		Ntab = 0
		read(iounit,*,iostat=ierr)
		if (ierr /= 0) then
			call alert(ierr,'dist_table_support: load_dist_table:: unable to read any lines')
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
		dt => winvn_dist_table
	
		dt% Ntable = Ntab
		! allocate the tables
		allocate(dt% P(Ntab), dt% mu_e(Ntab), dt% mu_n(Ntab), &
			&	 dt% Z(Nele), dt% A(Nele), dt% abun(Nele))

		! now read in the table, skipping first line
		read(iounit,*,iostat=ierr)
		if (ierr /= 0) then
			call alert(ierr,'dist_table_support: load_dist_table:: unable to read any lines')
			return
		end if
		
		do i = 1, dt% Ntable
			read(iounit,*,iostat=ierr) dt% P(i), dt% mu_e(i), dt% mu_n(i)
			if (ierr /=0) then
			   call alert(ierr,'dist_table_support: load_dist_table:: unable to read lines')
			   exit
		   end if
		end do
		close(iounit)
		call free_iounit(iounit)

!fix following

		! set up the bookmarks; table is assumed to be sorted!
		dt% Pmin = 1	   !minval(dt% P)
		dt% Pmax = dt% Ntable  !maxval(dt% P)
		
		! allocate the bookmarks
		Nele = dt% Pmax- dt% Pmin+1
		dt% Nelements = Nele
		allocate(dt% Pstart(Nele), dt% mu_n_min(Nele), dt% mu_n_max(Nele))
		
		!currentP = -1
	    currentP = 0.
		id_el = 0
		do i = 1, dt% Ntable
			if (dt% P(i) /= currentP) then	! found a new dist
				id_el = id_el + 1
				dt% Pstart(id_el) = i   ! Pstart can be an integer?
				dt% mu_n_min(id_el) = dt% mu_n(i)
				if (id_el > 1) then
					dt% mu_n_max(id_el-1) = dt% mu_n(i-1)
				end if
				currentP = dt% P(i)
			end if
		end do
		dt% mu_n_max(id_el) = dt% mu_n(dt% Ntable)

		if (ierr == 0) dist_table_is_loaded = .TRUE.
	end subroutine load_dist_table

	function dist_table_index(P, mu_n, ierr) result(id)
		real, intent(in) :: P, mu_n
		integer, intent(out) :: ierr
		integer :: id
		integer :: Pindex
		type(dist_table_type), pointer :: dt
	
		dt => winvn_dist_table
		ierr = 0
		id = -1
		if (P < dt% Pmin .or.  P > dt% Pmax) then
			ierr = -1
			return
		end if
		Pindex = P - dt% Pmin + 1
		if (mu_n < dt% mu_n_min(Pindex) .or. mu_n > dt% mu_n_max(Pindex)) then
			ierr = -2
			return
		end if
	
		id = dt% Pstart(Pindex) + mu_n - dt% mu_n_min(Pindex)
	end function dist_table_index

end module dist_table