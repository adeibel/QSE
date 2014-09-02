      program acc_crust

      use const_def
      use const_lib
      use num_def
      use mtx_lib
      use mtx_def
	  use outer_crust
	  use inner_crust
             
      implicit none
      logical, parameter :: show_all = .false.  ! false for releases
      integer ::  k, solver, omp_get_thread_num
      logical :: m_band, j_band, quiet

	  call follow_outer_crust_composition
      call follow_inner_crust_composition

      end program acc_crust