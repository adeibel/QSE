module ash_table
	! storage for the mass table
	type ash_table_type
		integer :: Zmin, Zmax	! minimum, maximum Z
		integer :: Nelements	! = Zmax-Zmin+1
		integer, dimension(:), pointer :: Zstart	! (Nelements) 
		integer, dimension(:), pointer :: Nmin, Nmax	! (Nelements)		
		integer :: Ntable	! number of table entries
		integer, dimension(:), pointer :: Z, N, A 
		real, dimension(:), pointer :: BE !binding energy 
		real, dimension(:), pointer :: Y !abundance fraction 
		real :: P !pressure
	end type ash_table_type
	logical, save :: ash_table_is_loaded = .FALSE.

	type(ash_table_type), target :: winvn_ash_table

contains

	subroutine ash_table_shutdown()
		type(ash_table_type), pointer :: at
		if (.not.ash_table_is_loaded) return
		at => winvn_ash_table
		deallocate(at% Zstart)
		deallocate(at% Nmin)
		deallocate(at% Nmax)
		deallocate(at% Z)
		deallocate(at% N)
		deallocate(at% A)
		deallocate(at% BE)
		deallocate(at% Y)	
		at% Ntable = 0
		at% Zmin = -1
		at% Zmax = -1		
		at% Nelements = 0
		ash_table_is_loaded = .FALSE.
	end subroutine ash_table_shutdown

	subroutine load_ash_table(datadir,datafile,ierr)
		use alert_lib
		use utils_lib
		character(len=*), intent(in) :: datadir,datafile
		integer, intent(out) :: ierr
		character(len=*), parameter :: default_dir = 'crust_eos'
		character(len=256) :: filename
		type(ash_table_type), pointer :: at
		integer :: i, Ntab, Nele, iounit
		integer :: currentZ, id_el

		if (ash_table_is_loaded) then
			ierr = 1
			call alert(1,'ash_table_support: load_ash_table:: table is already loaded')
			return
		end if
		
		iounit = alloc_iounit(ierr)
		if (ierr /= 0) then
			call alert(ierr,'ash_table_support: load_ash_table:: unable to alloc iounit')
			return
		end if
		filename = trim(datadir)//'/'//default_dir//'/'//trim(datafile)
		open(unit=iounit, file=trim(filename), iostat=ierr, action='read')
		if (ierr /= 0) then
			call alert(ierr,'ash_table_support: load_ash_table:: unable to open '//trim(filename))
			return
		end if

		! read through the file to get number of entries (skips first three lines)
		Ntab = 0
		read(iounit,*,iostat=ierr)
		read(iounit,*,iostat=ierr)
		read(iounit,*,iostat=ierr)
		if (ierr /= 0) then
			call alert(ierr,'ash_table_support: load_ash_table:: unable to read any lines')
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
		at => winvn_ash_table
		at% Ntable = Ntab

		! allocate the tables		
		allocate(at% Z(Ntab), at% N(Ntab), at% A(Ntab), at% BE(Ntab), at% Y(Ntab), &
			&		at% Zstart(Ntab), at% Nmin(Ntab), at% Nmax(Ntab))
		
		! now read in the table, skipping first three lines
		read(iounit,*,iostat=ierr)
		read(iounit,*,iostat=ierr) 
		read(iounit,*,iostat=ierr) at% P
		if (ierr /= 0) then
			call alert(ierr,'ash_table_support: load_ash_table:: unable to read any lines')
			return
		end if
		
		do i = 1, at% Ntable
			read(iounit,*,iostat=ierr) 	at% Z(i), at% N(i), at% A(i), at% BE(i), at% Y(i) 
			if (ierr /=0) then
			   call alert(ierr,'ash_table_support: load_ash_table:: unable to read lines')
			   exit
		   end if
		end do
		close(iounit)
		call free_iounit(iounit)


		if (ierr == 0) ash_table_is_loaded = .TRUE.
	end subroutine load_ash_table

end module ash_table
