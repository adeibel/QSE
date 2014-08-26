module ash_table
	! storage for the mass table
	type ash_table_type
		integer :: Ntable	! number of table entries
		real, dimension(:), pointer :: Z_ash !average Z of distribution
		real, dimension(:), pointer :: A_ash !average A of distribution
		real, dimension(:), pointer :: BE_ash !binding energy of ashes
		real, dimension(:), pointer :: Y_ash !abundance fraction 
	end type ash_table_type
	logical, save :: ash_table_is_loaded = .FALSE.

	type(ash_table_type), target :: winvn_ash_table

contains

	subroutine eos_table_shutdown()
		type(eos_table_type), pointer :: et
		if (.not.eos_table_is_loaded) return
		et => winvn_eos_table
		deallocate(et% nb)
		deallocate(et% rho)
		deallocate(et% Ye)
		deallocate(et% Yn)
		deallocate(et% fr)
		deallocate(et% Z_bar)
		deallocate(et% A_bar)
		deallocate(et% Q)
		deallocate(et% fr_x)
		deallocate(et% a_ionic)
		deallocate(et% heat)
		deallocate(et% heat_full)
		deallocate(et% nn)
		deallocate(et% pr)
		deallocate(et% gb)
		deallocate(et% ne)
		deallocate(et% mun)
		deallocate(et% mue)		
		et% Ntable = 0
	end subroutine eos_table_shutdown

	subroutine load_eos_table(datadir,datafile,ierr)
		use alert_lib
		use utils_lib
		character(len=*), intent(in) :: datadir,datafile
		integer, intent(out) :: ierr
		character(len=*), parameter :: default_dir = 'crust_eos'
		character(len=256) :: filename
		type(eos_table_type), pointer :: et
		integer :: i, Ntab, Nele, iounit
		integer :: currentZ, id_el

		if (eos_table_is_loaded) then
			ierr = 1
			call alert(1,'eos_table_support: load_eos_table:: table is already loaded')
			return
		end if
		
		iounit = alloc_iounit(ierr)
		if (ierr /= 0) then
			call alert(ierr,'eos_table_support: load_eos_table:: unable to alloc iounit')
			return
		end if
		filename = trim(datadir)//'/'//default_dir//'/'//trim(datafile)
		open(unit=iounit, file=trim(filename), iostat=ierr, action='read')
		if (ierr /= 0) then
			call alert(ierr,'eos_table_support: load_eos_table:: unable to open '//trim(filename))
			return
		end if

		! read through the file to get number of entries
		Ntab = 0
		read(iounit,*,iostat=ierr)
		read(iounit,*,iostat=ierr)
		if (ierr /= 0) then
			call alert(ierr,'eos_table_support: load_eos_table:: unable to read any lines')
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
		et => winvn_eos_table
		et% Ntable = Ntab

		! allocate the tables		
		allocate(et% nb(Ntab), et% rho(Ntab), et% Ye(Ntab), et% Yn(Ntab), et% fr(Ntab), &
			& et% Z_bar(Ntab), et% A_bar(Ntab), et% Q(Ntab), et% fr_x(Ntab), et% a_ionic(Ntab), &
			& et% heat(Ntab), et% heat_full(Ntab), et% nn(Ntab), et% pr(Ntab), et% gb(Ntab), &
			& et% ne(Ntab), et% mun(Ntab), et% mue(Ntab))
		
		! now read in the table, skipping first two lines
		read(iounit,*,iostat=ierr)
		read(iounit,*,iostat=ierr)
		if (ierr /= 0) then
			call alert(ierr,'eos_table_support: load_eos_table:: unable to read any lines')
			return
		end if
		
		do i = 1, et% Ntable
			read(iounit,*,iostat=ierr) 	et% nb(i), et% rho(i), et% Ye(i), et% Yn(i), et% fr(i), &
			& et% Z_bar(i), et% A_bar(i), et% Q(i), et% fr_x(i), et% a_ionic(i), &
			& et% heat(i), et% heat_full(i), et% nn(i), et% pr(i), et% gb(i), &
			& et% ne(i), et% mun(i), et% mue(i)
			if (ierr /=0) then
			   call alert(ierr,'eos_table_support: load_eos_table:: unable to read lines')
			   exit
		   end if
		end do
		close(iounit)
		call free_iounit(iounit)

		if (ierr == 0) eos_table_is_loaded = .TRUE.
	end subroutine load_eos_table

!	function eos_table_index(Z, N, ierr) result(id)
!		integer, intent(in) :: Z, N
!		integer, intent(out) :: ierr
!		integer :: id
!		integer :: Zindex
!		type(eos_table_type), pointer :: et
!	
!		et => winvn_eos_table
!		ierr = 0
!		id = -1
!		if (Z < et% Zmin .or. Z > et% Zmax) then
!			ierr = -1
!			return
!		end if
!		Zindex = Z - et% Zmin + 1
!		if (N < et% Nmin(Zindex) .or. N > et% Nmax(Zindex)) then
!			ierr = -2
!			return
!		end if
!	
!		id = et% Zstart(Zindex) + N - et% Nmin(Zindex)
!	end function eos_table_index

end module eos_table
