module constants
    implicit none
    ! Define constants (parameters in Fortran)
    real(8), parameter :: PI = 4.d0*DATAN(1.D0)
end module constants

program pareto_sampler
    implicit none
    
    ! Variable declarations
    integer :: n_samples, i, seed_size
    real(8) :: alpha, xm
    real(8), allocatable :: samples(:)
    integer, allocatable :: seed(:)
    
    ! Pareto distribution parameters
    ! alpha: shape parameter (must be > 0)
    ! xm: scale parameter (minimum value, must be > 0)
    alpha = 2.0d0
    xm = 1.0d0
    n_samples = 10000
    
    ! Allocate array for samples
    allocate(samples(n_samples))
    
    ! Initialize random seed
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 12345  ! Fixed seed for reproducibility; remove for different results each run
    call random_seed(put=seed)
    
    ! Generate samples from Pareto distribution
    do i = 1, n_samples
        samples(i) = sample_pareto(alpha, xm)
    end do
    
    ! Display first 20 samples
    print *, "Pareto Distribution Sampler"
    print *, "Parameters: alpha =", alpha, ", xm =", xm
    print *, "Number of samples:", n_samples
    print *, ""
    print *, "First 20 samples:"
    do i = 1, min(20, n_samples)
        print '(I5, A, F12.6)', i, ":", samples(i)
    end do
    
    ! Calculate and display statistics
    print *, ""
    print *, "Sample Statistics:"
    print *, "Mean:", sum(samples) / n_samples
    print *, "Min:", minval(samples)
    print *, "Max:", maxval(samples)
    
    ! Theoretical mean for Pareto distribution (when alpha > 1)
    if (alpha > 1.0d0) then
        print *, "Theoretical mean:", alpha * xm / (alpha - 1.0d0)
    end if
    
    ! Clean up
    deallocate(samples)
    deallocate(seed)
    
contains
    
    ! Function to sample from Pareto distribution using inverse transform method
    ! Pareto PDF: f(x) = alpha * xm^alpha / x^(alpha+1) for x >= xm
    ! Inverse CDF: x = xm / u^(1/alpha) where u ~ Uniform(0,1)
    function p(mu) result(x)
        implicit none
        real(8), intent(in) :: mu
        real(8) :: x, u
        
        ! Generate uniform random number
        call random_number(u)
        
        ! Apply inverse transform
        x = 1 / (u**(1.0d0/(mu + 1)))
        
    end function p

    function intp(r1, t1, r2, t2, t) result(r)
        implicit none
        real(8), intent(in) :: r1, t1, r2, t2, t
        r = (t2 - t)/(t2 - t1)*r1 + (t - t1)/(t2 - t1)*r2
    end function intp

    subroutine register_time(ct, i, t, nt, r, t_array, store_t)
        implicit none
        integer :: ct, i, j
        real(8) :: t, nt, x, y
        real(8), intent(in) :: t_array(:)
        real(8) :: store_t(:, :)
        do while (ct <= size(store_t) )
            if (t <= store_t(ct) .and. store_t(ct) <= nt) then
                store_t(i, ct) = r
                ct += 1
            else
                exit
            end if
        end do
    end subroutine register_time
        
    function walk(n, mu, t_array) result(mad_array)
        use constants
        implicit none
        integer :: n
        real(8), intent(in) :: t_array(:)
        real(8), intent(in) :: mu
        real(8), dimension(n, size(t_array)) :: store_t
        real(8), dimension(n) :: mad_array, ct, x, y
        integer, dimension(n) :: tcurr
        integer :: i
        real(8) :: nx, ny, phi, l, nt
        logical :: done_something
        do
            done_something = .FALSE.
            do i = 1, n
                if (ct(i) > t(size(t))) cycle
                done_something = .TRUE.
                call random_number(phi)
                phi = 2d0*PI*phi
                l = p(mu)
                nx = x(i) + l*cos(phi)
                ny = y(i) + l*sin(phi)
                r = sqrt(nx*nx + ny*ny) 
                register_time(ct(i), i, t, t + r, r, t_array, store_t)
            end do
        end do
    end function walk
    
end program pareto_sampler
