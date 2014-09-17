module dist_table
	! storage for the mass table
	type dist_table_type
		integer :: Zmin, Zmax	! minimum, maximum Z
		integer :: Nelements	! = Zmax-Zmin+1
		integer, dimension(:), pointer :: Zstart	! (Nelements) 
		integer, dimension(:), pointer :: Nmin, Nmax	! (Nelements)		
		integer :: Ntable	! number of table entries
		integer, dimension(:), pointer :: Z, N, A 
		real, dimension(:), pointer :: BE !binding energy 
		real, dimension(:), pointer :: Y !abundance fraction 
		real :: P !pressure
	end type dist_table_type
	integer :: index
	integer, dimension(:), allocatable :: Z, N, A 
	real, dimension(:), allocatable :: BE !binding energy 
	real, dimension(:), allocatable :: YY !abundance fraction 	
	logical, dimension(:), allocatable :: abundant_nucleus
	logical, save :: dist_table_is_loaded = .FALSE.

	type(dist_table_type), target :: winvn_dist_table

contains

	subroutine dist_table_shutdown()
		type(dist_table_type), pointer :: dt
		if (.not.dist_table_is_loaded) return
		dt => winvn_dist_table
		deallocate(dt% Zstart)
		deallocate(dt% Nmin)
		deallocate(dt% Nmax)
		deallocate(dt% Z)
		deallocate(dt% N)
		deallocate(dt% A)
		deallocate(dt% BE)
		deallocate(dt% Y)	
		dt% Ntable = 0
		dt% Zmin = -1
		dt% Zmax = -1		
		dt% Nelements = 0
		dist_table_is_loaded = .FALSE.
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
		integer :: currentZ, id_el

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

		! read through the file to get number of entries (skips first three lines)
		Ntab = 0
		read(iounit,*,iostat=ierr)
		read(iounit,*,iostat=ierr)
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
		allocate(Z(Ntab), N(Ntab), A(Ntab), BE(Ntab), YY(Ntab))
		
		! now read in the table, skipping first three lines
		read(iounit,*,iostat=ierr)
		read(iounit,*,iostat=ierr) 
		read(iounit,*,iostat=ierr) dt% P
		if (ierr /= 0) then
			call alert(ierr,'dist_table_support: load_dist_table:: unable to read any lines')
			return
		end if
		
		do i = 1, dt% Ntable
			read(iounit,*,iostat=ierr) 	Z(i), N(i), A(i), BE(i), YY(i) 
			if (ierr /=0) then
			   call alert(ierr,'dist_table_support: load_dist_table:: unable to read lines')
			   exit
		   end if
		end do
		close(iounit)
		call free_iounit(iounit)

		allocate(abundant_nucleus(Ntab))

		! set Y cut off, reallocate table accordingly
		do i = 1, dt% Ntable
			if (YY(i) > 1.d-20) then
			abundant_nucleus(i) = .TRUE.
			else
			abundant_nucleus(i) = .FALSE.
			end if
		end do			

		Ntab = count(abundant_nucleus)
		allocate(dt% Z(Ntab), dt% N(Ntab), dt% A(Ntab), dt% BE(Ntab), dt% Y(Ntab), &
			&		dt% Zstart(Ntab), dt% Nmin(Ntab), dt% Nmax(Ntab))
		
		index = 1		
		do i = 1, dt% Ntable
			if (abundant_nucleus(i) .eqv. .TRUE.) then
			dt% Z(index) = Z(i)
			dt% N(index) = N(i)
			dt% A(index) = A(i)
			dt% BE(index) = BE(i)
			dt% Y(index) = YY(i)
			index = index+1
			end if
	    end do		
	    dt% Ntable = count(abundant_nucleus)
				
		if (ierr == 0) dist_table_is_loaded = .TRUE.
	end subroutine load_dist_table

end module dist_table
