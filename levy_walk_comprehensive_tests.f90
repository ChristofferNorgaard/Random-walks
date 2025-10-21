module constants
    implicit none
    real(8), parameter :: PI = 4.d0*DATAN(1.D0)
end module constants

module levy_walk_module
    use constants
    implicit none
contains
    
    ! Power-law distribution: x = 1 / u^(1/(mu+1))
    function p(mu) result(x)
        real(8), intent(in) :: mu
        real(8) :: x, u
        call random_number(u)
        x = 1.0d0 / (u**(1.0d0/(mu + 1.0d0)))
    end function p
    
    ! Levy walk simulation
    function walk(n, mu, t_array) result(mad_array)
        integer, intent(in) :: n
        real(8), intent(in) :: mu, t_array(:)
        real(8), dimension(size(t_array)) :: mad_array
        real(8), dimension(n, size(t_array)) :: r_at_t
        real(8), dimension(n) :: x, y, t_current
        integer :: i, j, step
        real(8) :: phi, step_length, t_next, r
        integer, parameter :: max_steps = 100000
        
        x = 0.0d0
        y = 0.0d0
        t_current = 0.0d0
        r_at_t = 0.0d0
        
        do i = 1, n
            do step = 1, max_steps
                if (t_current(i) > t_array(size(t_array))) exit
                
                step_length = p(mu)
                call random_number(phi)
                phi = 2.0d0 * PI * phi
                
                t_next = t_current(i) + step_length
                x(i) = x(i) + step_length * cos(phi)
                y(i) = y(i) + step_length * sin(phi)
                r = sqrt(x(i)**2 + y(i)**2)
                
                do j = 1, size(t_array)
                    if (t_current(i) <= t_array(j) .and. t_array(j) <= t_next) then
                        r_at_t(i, j) = r
                    end if
                end do
                
                t_current(i) = t_next
            end do
        end do
        
        do j = 1, size(t_array)
            mad_array(j) = sum(r_at_t(:, j)) / real(n, 8)
        end do
    end function walk
    
end module levy_walk_module

program test_levy_walk
    use constants
    use levy_walk_module
    implicit none
    
    call test_distribution()
    call test_walk_monotonicity()
    call test_walk_initial_conditions()
    call test_walk_mu_dependence()
    call test_distribution_moments()
    
contains
    
    subroutine test_distribution()
        integer, parameter :: n_samples = 50000
        real(8) :: mu, samples(n_samples), mean_val, variance, min_val, max_val
        integer :: i, count_below_10
        integer :: seed(10)
        
        print *, "==========================================="
        print *, "TEST 1: Distribution p(mu) Properties"
        print *, "==========================================="
        print *, "Formula: x = 1 / u^(1/(mu+1))"
        print *, ""
        
        seed = 54321
        call random_seed(put=seed)
        
        mu = 2.0d0
        
        do i = 1, n_samples
            samples(i) = p(mu)
        end do
        
        mean_val = sum(samples) / n_samples
        variance = sum((samples - mean_val)**2) / n_samples
        min_val = minval(samples)
        max_val = maxval(samples)
        count_below_10 = count(samples < 10.0d0)
        
        print *, "  mu =", mu
        print *, "  Number of samples:", n_samples
        print *, "  Mean:", mean_val
        print *, "  Std dev:", sqrt(variance)
        print *, "  Min:", min_val
        print *, "  Max:", max_val
        print *, "  Fraction < 10:", real(count_below_10)/n_samples
        
        ! Check that all samples >= 1 (theoretical minimum)
        if (min_val >= 0.999d0) then
            print *, "  ✓ PASS: All samples >= 1 (as expected)"
        else
            print *, "  ✗ FAIL: Found samples < 1"
        end if
        
        ! For mu=2, theoretical mean = (mu+1)/mu = 1.5
        if (abs(mean_val - 1.5d0) / 1.5d0 < 0.05d0) then
            print *, "  ✓ PASS: Mean within 5% of theoretical value (1.5)"
        else
            print *, "  ✗ WARNING: Mean differs from theoretical (1.5)"
            print *, "    Difference:", abs(mean_val - 1.5d0)
        end if
        
        print *, ""
    end subroutine test_distribution
    
    subroutine test_walk_monotonicity()
        integer :: n_walkers, n_times, i
        real(8) :: mu
        real(8), allocatable :: t_array(:), mad_array(:)
        logical :: is_monotonic
        
        print *, "==========================================="
        print *, "TEST 2: Walk MAD Monotonicity"
        print *, "==========================================="
        
        n_walkers = 100
        n_times = 30
        mu = 2.0d0
        
        allocate(t_array(n_times))
        do i = 1, n_times
            t_array(i) = real(i, 8) * 0.5d0
        end do
        
        mad_array = walk(n_walkers, mu, t_array)
        
        ! Check monotonicity (with small tolerance for numerical issues)
        is_monotonic = .true.
        do i = 2, n_times
            if (mad_array(i) < mad_array(i-1) - 1.0d-6) then
                is_monotonic = .false.
                exit
            end if
        end do
        
        print *, "  n_walkers =", n_walkers
        print *, "  mu =", mu
        print *, "  Time range: 0 to", t_array(n_times)
        print *, ""
        print *, "  MAD values (every 5th point):"
        do i = 1, n_times, 5
            print '(A, F8.2, A, F12.6)', "    t = ", t_array(i), " : MAD = ", mad_array(i)
        end do
        
        if (is_monotonic) then
            print *, "  ✓ PASS: MAD is monotonically increasing"
        else
            print *, "  ✗ WARNING: MAD not perfectly monotonic"
            print *, "    (Can occur with few walkers or statistical fluctuations)"
        end if
        
        deallocate(t_array, mad_array)
        print *, ""
    end subroutine test_walk_monotonicity
    
    subroutine test_walk_initial_conditions()
        integer :: n_walkers, n_times, i
        real(8) :: mu
        real(8), allocatable :: t_array(:), mad_array(:)
        
        print *, "==========================================="
        print *, "TEST 3: Walk Initial Conditions"
        print *, "==========================================="
        
        n_walkers = 50
        n_times = 10
        mu = 1.5d0
        
        allocate(t_array(n_times))
        t_array(1) = 0.0d0
        do i = 2, n_times
            t_array(i) = real(i-1, 8)
        end do
        
        mad_array = walk(n_walkers, mu, t_array)
        
        print *, "  n_walkers =", n_walkers
        print *, "  mu =", mu
        print *, "  Initial MAD (t=0):", mad_array(1)
        print *, "  MAD at t=1:", mad_array(2)
        print *, "  MAD at t=", t_array(n_times), ":", mad_array(n_times)
        
        if (mad_array(1) < 0.01d0) then
            print *, "  ✓ PASS: Initial MAD ≈ 0 (walkers start at origin)"
        else
            print *, "  ✗ FAIL: Initial MAD should be near 0"
        end if
        
        if (mad_array(n_times) > mad_array(1)) then
            print *, "  ✓ PASS: MAD increases over time"
        else
            print *, "  ✗ FAIL: MAD should increase"
        end if
        
        deallocate(t_array, mad_array)
        print *, ""
    end subroutine test_walk_initial_conditions
    
    subroutine test_walk_mu_dependence()
        integer :: n_walkers, n_times, i
        real(8) :: mu1, mu2, mu3
        real(8), allocatable :: t_array(:), mad1(:), mad2(:), mad3(:)
        
        print *, "==========================================="
        print *, "TEST 4: Walk Dependence on mu"
        print *, "==========================================="
        
        n_walkers = 150
        n_times = 15
        mu1 = 1.0d0
        mu2 = 2.0d0
        mu3 = 3.0d0
        
        allocate(t_array(n_times))
        do i = 1, n_times
            t_array(i) = real(i, 8)
        end do
        
        print *, "  Comparing mu = 1.0, 2.0, 3.0"
        print *, "  (Lower mu = heavier tail = longer steps)"
        print *, ""
        
        mad1 = walk(n_walkers, mu1, t_array)
        mad2 = walk(n_walkers, mu2, t_array)
        mad3 = walk(n_walkers, mu3, t_array)
        
        print *, "  MAD comparison:"
        print *, "      Time     mu=1.0     mu=2.0     mu=3.0"
        do i = 1, n_times, 3
            print '(A, F8.2, 3F12.4)', "    ", t_array(i), mad1(i), mad2(i), mad3(i)
        end do
        
        ! Lower mu should give larger steps on average
        if (mad1(n_times) > mad2(n_times) .and. mad2(n_times) > mad3(n_times)) then
            print *, "  ✓ PASS: MAD decreases with increasing mu (expected)"
        else
            print *, "  ✗ WARNING: Expected MAD(mu=1) > MAD(mu=2) > MAD(mu=3)"
            print *, "    (Statistical fluctuations may cause deviations)"
        end if
        
        deallocate(t_array, mad1, mad2, mad3)
        print *, ""
    end subroutine test_walk_mu_dependence
    
    subroutine test_distribution_moments()
        integer, parameter :: n_samples = 100000
        real(8) :: mu, samples(n_samples), theoretical_mean, empirical_mean
        integer :: i
        integer :: seed(10)
        
        print *, "==========================================="
        print *, "TEST 5: Distribution Moments"
        print *, "==========================================="
        
        seed = 99999
        call random_seed(put=seed)
        
        print *, "  Testing theoretical moments of p(mu)"
        print *, "  For x = 1/u^(1/(mu+1)), E[x] = (mu+1)/mu for mu > 0"
        print *, ""
        
        ! Test different mu values
        do i = 1, 3
            if (i == 1) mu = 1.0d0
            if (i == 2) mu = 2.0d0
            if (i == 3) mu = 4.0d0
            
            call sample_distribution(mu, n_samples, samples)
            empirical_mean = sum(samples) / n_samples
            theoretical_mean = (mu + 1.0d0) / mu
            
            print '(A, F5.2)', "  mu = ", mu
            print '(A, F10.6)', "    Theoretical mean: ", theoretical_mean
            print '(A, F10.6)', "    Empirical mean:   ", empirical_mean
            print '(A, F8.4, A)', "    Relative error:   ", &
                abs(empirical_mean - theoretical_mean)/theoretical_mean * 100, "%"
            
            if (abs(empirical_mean - theoretical_mean)/theoretical_mean < 0.02d0) then
                print *, "    ✓ PASS: Within 2% of theory"
            else
                print *, "    ✗ WARNING: Differs by >2% from theory"
            end if
            print *, ""
        end do
        
    end subroutine test_distribution_moments
    
    subroutine sample_distribution(mu, n, samples)
        real(8), intent(in) :: mu
        integer, intent(in) :: n
        real(8), intent(out) :: samples(n)
        integer :: i
        
        do i = 1, n
            samples(i) = p(mu)
        end do
    end subroutine sample_distribution
    
end program test_levy_walk
