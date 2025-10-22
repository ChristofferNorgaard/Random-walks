program Hello

  Use omp_lib

  Implicit None

  INTEGER :: nthreads
  nthreads = 4

  CALL OMP_SET_NUM_THREADS(nthreads)

  write(*,*) omp_get_num_procs()
  write(*,*) omp_get_max_threads()
  write(*,*) omp_get_num_threads()

  !$OMP PARALLEL
    PRINT *, "Hello from process: ", OMP_GET_THREAD_NUM()
  !$OMP END PARALLEL

end program Hello
