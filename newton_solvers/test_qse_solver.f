
      program test_qse_solver
      
      use qse_solver              
      use const_def
      use const_lib
      use num_def
      use mtx_lib
      use mtx_def
      
      implicit none
      logical, parameter :: show_all = .false.  ! false for releases
      integer ::  k, solver, omp_get_thread_num
      logical :: m_band, j_band, quiet

      call do_test_newton

      end program test_qse_solver