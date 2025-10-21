module weighted_fit_module
    implicit none
    
contains
    subroutine combine_series(y_series, n, nseries, y_mean, y_sigma, status)
        ! Combines multiple series to calculate mean and standard deviation
        ! at each time point
        !
        ! Inputs:
        !   y_series(n, nseries) : array of series data
        !   n                    : number of data points per series
        !   nseries              : number of series
        !
        ! Outputs:
        !   y_mean(n)   : mean at each time point
        !   y_sigma(n)  : standard deviation at each time point
        !   status      : 0 if successful, -1 if nseries < 2
        
        implicit none
        
        integer, intent(in) :: n, nseries
        real(8), intent(in) :: y_series(n, nseries)
        real(8), intent(out) :: y_mean(n), y_sigma(n)
        integer, intent(out) :: status
        
        integer :: i, j
        real(8) :: sum_val, sum_sq, variance
        
        status = 0
        
        if (nseries < 2) then
            status = -1
            return
        end if
        
        do i = 1, n
            sum_val = 0.0d0
            sum_sq = 0.0d0
            
            ! Calculate mean
            do j = 1, nseries
                sum_val = sum_val + y_series(i, j)
            end do
            y_mean(i) = sum_val / real(nseries, 8)
            
            ! Calculate standard deviation
            do j = 1, nseries
                sum_sq = sum_sq + (y_series(i, j) - y_mean(i))**2
            end do
            
            ! Use sample standard deviation (divide by nseries-1)
            variance = sum_sq / real(nseries - 1, 8)
            y_sigma(i) = sqrt(variance)
            
            ! Handle case where all values are identical (sigma = 0)
            if (y_sigma(i) < 1.0d-15) then
                y_sigma(i) = 1.0d-15  ! Small non-zero value
            end if
        end do
        
    end subroutine combine_series

        subroutine weighted_linear_fit(x, y, sigma, n, intercept, slope, &
                                   intercept_err, slope_err, status)
        ! Performs weighted linear regression: y = intercept + slope * x
        !
        ! Inputs:
        !   x(n)     : independent variable data
        !   y(n)     : dependent variable data
        !   sigma(n) : errors on y values
        !   n        : number of data points
        !
        ! Outputs:
        !   intercept     : fitted intercept
        !   slope         : fitted slope
        !   intercept_err : error on intercept
        !   slope_err     : error on slope
        !   chi2          : chi-squared statistic
        !   status        : 0 if successful, non-zero otherwise
        
        implicit none
        
        ! Arguments
        integer, intent(in) :: n
        real(8), intent(in) :: x(n), y(n), sigma(n)
        real(8), intent(out) :: intercept, slope
        real(8), intent(out) :: intercept_err, slope_err
        integer, intent(out) :: status
        
        ! Local variables
        integer, parameter :: nparam = 2
        real(8), dimension(n, nparam) :: A_scaled, A_work
        real(8), dimension(n, 1) :: y_scaled_2d, y_work_2d
        real(8), dimension(n) :: y_scaled
        real(8), dimension(nparam, nparam) :: covar
        real(8), allocatable :: work(:)
        integer :: i, info, lwork
        real(8) :: work_query(1), y_fit
        
        ! External LAPACK routines
        external :: dgels, dpotrf, dpotri
        
        ! Initialize status
        status = 0
        
        ! Check for valid input
        if (n < nparam) then
            status = -1
            return
        end if
        
        ! Create scaled design matrix [1/sigma, x/sigma] and scaled y
        do i = 1, n
            if (sigma(i) <= 0.0d0) then
                status = -2  ! Invalid sigma
                return
            end if
            A_scaled(i, 1) = 1.0d0 / sigma(i)
            A_scaled(i, 2) = x(i) / sigma(i)
            y_scaled(i) = y(i) / sigma(i)
            y_scaled_2d(i, 1) = y(i) / sigma(i)
        end do
        
        ! Make working copies (DGELS overwrites inputs)
        A_work = A_scaled
        y_work_2d = y_scaled_2d
        
        ! Query optimal workspace size
        lwork = -1
        call dgels('N', n, nparam, 1, A_work, n, y_work_2d, n, &
                   work_query, lwork, info)
        
        if (info /= 0) then
            status = 1
            return
        end if
        
        lwork = int(work_query(1)) + 1
        allocate(work(lwork))
        
        ! Solve weighted least squares
        call dgels('N', n, nparam, 1, A_work, n, y_work_2d, n, &
                   work, lwork, info)
        
        if (info /= 0) then
            status = 2
            deallocate(work)
            return
        end if
        
        ! Extract coefficients
        intercept = y_work_2d(1, 1)
        slope = y_work_2d(2, 1)
        
        deallocate(work)
        
        ! Calculate covariance matrix: (A^T * A)^(-1)
        covar = 0.0d0
        do i = 1, n
            covar(1, 1) = covar(1, 1) + A_scaled(i, 1) * A_scaled(i, 1)
            covar(1, 2) = covar(1, 2) + A_scaled(i, 1) * A_scaled(i, 2)
            covar(2, 1) = covar(2, 1) + A_scaled(i, 2) * A_scaled(i, 1)
            covar(2, 2) = covar(2, 2) + A_scaled(i, 2) * A_scaled(i, 2)
        end do
        
        ! Invert using Cholesky decomposition
        call dpotrf('U', nparam, covar, nparam, info)
        if (info /= 0) then
            status = 3
            return
        end if
        
        call dpotri('U', nparam, covar, nparam, info)
        if (info /= 0) then
            status = 4
            return
        end if
        
        ! Copy upper triangle to lower (symmetric matrix)
        covar(2, 1) = covar(1, 2)
        
        ! Extract parameter errors
        intercept_err = sqrt(covar(1, 1))
        slope_err = sqrt(covar(2, 2))
        
    end subroutine weighted_linear_fit

end module weighted_fit_module

module constants
    implicit none
    real(8), parameter :: PI = 4.d0*DATAN(1.D0)
end module constants

program pareto_sampler
    use constants
    use weighted_fit_module
    implicit none
    
    ! Variable declarations
    integer :: n_samples, i, j, seed_size, nseries
    real(8) :: alpha, xm
    integer, allocatable :: seed(:)
    real(8), allocatable :: t_array(:), mad_array(:), y_mean(:), y_sigma(:), t_log_array(:)
    real(8), allocatable :: y_series(:, :)
    integer :: n_particles, n_times, io, status
    real(8) :: intercept, slope, intercept_err, slope_err, mu
    
    nseries = 100
    ! Distribution parameter
    ! p(mu) gives: x = 1 / u^(1/(mu+1))
    alpha = 2.0d0
    xm = 1.0d0
    
    ! Initialize random seed
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 12345  ! Fixed seed for reproducibility
    call random_seed(put=seed)
    
    ! Test Levy walk simulation
    print *, ""
    print *, "=== Levy Walk Simulation ==="
    n_particles = 10000
    n_times = 200
    allocate(y_series(n_times, nseries))
    allocate(y_mean(n_times))
    allocate(y_sigma(n_times))

    allocate(t_array(n_times))
    allocate(mad_array(n_times))
    do i = 1, n_times
        t_array(i) = real(i, 8) * 100d0
    end do

    t_log_array = log(t_array)

    open(newunit=io, file="slopes.txt")
    do j = 5, 100
        mu = 1d0 + 3d0* real(j, 8) / 100d0
        print *,"Mu now: ", mu
        do i = 1, nseries 
            y_series(:, i) = walk(n_particles, mu, t_array)
            y_series(:, i) = 2*log(y_series(:, i))
        end do

        call combine_series(y_series, n_times, nseries, y_mean, y_sigma, status)

        call weighted_linear_fit(t_log_array, y_mean, y_sigma, n_times, &
                             intercept, slope, &
                             intercept_err, slope_err, &
                             status)
        ! Clean up
        write(io, *) mu, slope, slope_err
    end do
    close(io)
    deallocate(seed, t_array)
    
contains
    
    ! Function to sample from power-law distribution
    ! Formula: x = 1 / u^(1/(mu+1)) where u ~ Uniform(0,1)
    function p(mu) result(x)
        implicit none
        real(8), intent(in) :: mu
        real(8) :: x, u
        
        ! Generate uniform random number
        call random_number(u)
        
        ! Apply transform
        x = 1.0d0 / (u**(1.0d0/(mu - 1.0d0)))
        
    end function p

    function mad(r) result (m)
        real(8), intent(in) :: r(:)
        real(8), allocatable :: c(:)
        real(8) :: m
        integer :: n, info
        
        n = size(r)
        
        ! Handle empty or invalid arrays
        if (n <= 0) then
            m = 0.0d0
            return
        end if
        
        ! Allocate and copy
        allocate(c(n))
        c = r
        
        ! Sort
        call DLASRT('I', n, c, info)
        
        if (info /= 0) then
            print *, "DLASRT error, info =", info
            m = 0.0d0
            deallocate(c)
            return
        end if
        
        ! Calculate median
        if (mod(n, 2) == 0) then
            m = (c(n/2) + c(n/2+1)) / 2.0d0
        else
            m = c((n+1)/2)
        end if
       
        c = abs(c - m)

        call DLASRT('I', n, c, info)
        
        if (info /= 0) then
            print *, "DLASRT error, info =", info
            m = 0.0d0
            deallocate(c)
            return
        end if
        
        ! Calculate median
        if (mod(n, 2) == 0) then
            m = (c(n/2) + c(n/2+1)) / 2.0d0
        else
            m = c((n+1)/2)
        end if

        ! Clean up
        deallocate(c)
        
    end function mad

    function intp(q1, t1, q2, t2, t) result(q)
        implicit none
        real(8), intent(in) :: q1, t1, q2, t2, t
        real(8) :: q
        q = (t2 - t)/(t2 - t1)*q1 + (t - t1)/(t2 - t1)*q2
    end function intp

    ! Levy walk simulation
    function walk(n, mu, t_array) result(mad_array)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: mu
        real(8), intent(in) :: t_array(:)
        real(8), dimension(n, size(t_array)) :: r_at_t
        real(8), dimension(n) :: x, y, t_current, mad_array
        integer, dimension(n) :: ct  ! Current time index for each particle
        integer :: i, j, step
        real(8) :: phi, step_length, t_next, r, xn, yn, xc, yc
        integer, parameter :: max_steps = 100000
        
        ! Initialize
        x = 0.0d0
        y = 0.0d0
        t_current = 0.0d0
        r_at_t = 0.0d0
        ct = 1  ! Start at first time index
        
        ! Simulate each particle
        do i = 1, n
            do step = 1, max_steps
                ! Exit if we've passed the last time point
                if (ct(i) > size(t_array)) exit
                
                ! Sample step length from distribution p(mu)
                step_length = p(mu)
                
                ! Sample random direction uniformly
                call random_number(phi)
                phi = 2.0d0 * PI * phi
                
                ! Calculate next time (ballistic motion: distance = time with unit velocity)
                t_next = t_current(i) + step_length
                
                ! Update position after completing the step
                xn = x(i) + step_length * cos(phi)
                yn = y(i) + step_length * sin(phi)
                
                ! Calculate distance from origin
                r = sqrt(x(i)**2 + y(i)**2)
                
                ! Record distance only for time points we haven't passed yet
                ! Start from ct(i) instead of 1 to avoid rechecking old times
                do while (ct(i) <= size(t_array))
                    if (t_current(i) <= t_array(ct(i)) .and. t_array(ct(i)) <= t_next) then
                        xc = intp(x(i), t_current(i), xn, t_next, t_array(ct(i)))
                        yc = intp(y(i), t_current(i), yn, t_next, t_array(ct(i)))
                        r_at_t(i, ct(i)) = sqrt(xc**2 + yc**2)
                        ct(i) = ct(i) + 1
                    else
                        exit
                    end if
                end do
                
                ! Move to next time
                t_current(i) = t_next
                y(i) = yn
                x(i) = xn
            end do
        end do
        do j = 1, size(t_array)
            mad_array(j) = mad(r_at_t(:, j))
        end do
    end function walk

    function flight(n, mu, t_array) result(mad_array)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: mu
        real(8), intent(in) :: t_array(:)
        real(8), dimension(n, size(t_array)) :: r_at_t
        real(8), dimension(n) :: x, y, t_current, mad_array
        integer, dimension(n) :: ct  ! Current time index for each particle
        integer :: i, j, step
        real(8) :: phi, step_length, t_next, r, yn, xn, xc, yc
        integer, parameter :: max_steps = 100000
        
        ! Initialize
        x = 0.0d0
        y = 0.0d0
        t_current = 0.0d0
        r_at_t = 0.0d0
        ct = 1  ! Start at first time index
        
        ! Simulate each particle
        do i = 1, n
            do step = 1, max_steps
                ! Exit if we've passed the last time point
                if (ct(i) > size(t_array)) exit
                
                ! Sample step length from distribution p(mu)
                step_length = p(mu)
                
                ! Sample random direction uniformly
                call random_number(phi)
                phi = 2.0d0 * PI * phi
                
                ! Calculate next time (flight motion: time step is constant)
                t_next = t_current(i) + 1 
                
                ! Update position after completing the step
                xn = x(i) + step_length * cos(phi)
                yn = y(i) + step_length * sin(phi)
                
                ! Calculate distance from origin
                r = sqrt(x(i)**2 + y(i)**2)
                
                ! Record distance only for time points we haven't passed yet
                ! Start from ct(i) instead of 1 to avoid rechecking old times
                do while (ct(i) <= size(t_array))
                    if (t_current(i) <= t_array(ct(i)) .and. t_array(ct(i)) <= t_next) then
                        xc = intp(x(i), t_current(i), xn, t_next, t_array(ct(i)))
                        yc = intp(y(i), t_current(i), yn, t_next, t_array(ct(i)))
                        r_at_t(i, ct(i)) = sqrt(xc**2 + yc**2)
                        ct(i) = ct(i) + 1
                    else
                        exit
                    end if
                end do
                
                ! Move to next time
                t_current(i) = t_next
                y(i) = yn
                x(i) = xn
            end do
        end do
        do j = 1, size(t_array)
            mad_array(j) = mad(r_at_t(:, j))
        end do
    end function flight
    

end program pareto_sampler
