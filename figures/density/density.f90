module constants
    implicit none
    real(8), parameter :: PI = 4.d0*DATAN(1.D0)
end module constants

program pareto_sampler
    use constants
    implicit none
    
    ! Variable declarations
    integer :: n_samples, i, seed_size
    real(8) :: alpha, xm
    real(8), allocatable :: samples(:)
    integer, allocatable :: seed(:)
    real(8), allocatable :: t_array(:), pos_array(:, :, :)
    integer :: n_particles, n_times, io
    
    ! Distribution parameter
    ! p(mu) gives: x = 1 / u^(1/(mu+1))
    alpha = 2.0d0
    xm = 1.0d0
    n_samples = 10000
    
    ! Allocate array for samples
    allocate(samples(n_samples))
    
    ! Initialize random seed
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 12345  ! Fixed seed for reproducibility
    call random_seed(put=seed)
    
    ! Test Levy walk simulation
    print *, ""
    print *, "=== Levy Walk Simulation ==="
    n_particles = 1000
    n_times = 20
    
    allocate(t_array(n_times))
    do i = 1, n_times
        t_array(i) = real(i-1, 8) * 100d0
    end do

    open(newunit=io, file="w15.dat")
    pos_array = walk(n_particles, 1.5d0, t_array)
    do i = 1, n_particles
        write(io, *) pos_array(i, 20, 1), pos_array(i, 20, 2)
    end do
    close(io) 
    open(newunit=io, file="w25.dat")
    pos_array = walk(n_particles, 2.5d0, t_array)
    do i = 1, n_particles
        write(io, *) pos_array(i, 20, 1), pos_array(i, 20, 2)
    end do
    close(io)
    open(newunit=io, file="w35.dat")
    pos_array = walk(n_particles, 3.5d0, t_array)
    do i = 1, n_particles
        write(io, *) pos_array(i, 20, 1), pos_array(i, 20, 2)
    end do
    close(io)

    open(newunit=io, file="f15.dat")
    pos_array = flight(n_particles, 1.5d0, t_array)
    do i = 1, n_particles
        write(io, *) pos_array(i, 20, 1), pos_array(i, 20, 2)
    end do
    close(io) 
    open(newunit=io, file="f25.dat")
    pos_array = flight(n_particles, 2d0, t_array)
    do i = 1, n_particles
        write(io, *) pos_array(i, 20, 1), pos_array(i, 20, 2)
    end do
    close(io)
    open(newunit=io, file="f35.dat")
    pos_array = flight(n_particles, 3.5d0, t_array)
    do i = 1, n_particles
        write(io, *) pos_array(i, 20, 1), pos_array(i, 20, 2)
    end do
    close(io)
    ! Clean up
    deallocate(samples, seed, t_array, pos_array)
    
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
    function walk(n, mu, t_array) result(pos_at_t)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: mu
        real(8), intent(in) :: t_array(:)
        real(8), dimension(n, size(t_array), 2) :: pos_at_t
        real(8), dimension(n) :: x, y, t_current
        real(8), dimension(n, 2) :: pos 
        integer, dimension(n) :: ct  ! Current time index for each particle
        integer :: i, j, step
        real(8) :: phi, step_length, t_next, r, xn, yn
        integer, parameter :: max_steps = 100000
        
        ! Initialize
        x = 0.0d0
        y = 0.0d0
        t_current = 0.0d0
        pos_at_t = 0.0d0
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
                        pos_at_t(i, ct(i), 1) = intp(x(i), t_current(i), xn, t_next, t_array(ct(i)))
                        pos_at_t(i, ct(i), 2) = intp(y(i), t_current(i), yn, t_next, t_array(ct(i)))
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
        
    end function walk

    function flight(n, mu, t_array) result(pos_at_t)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: mu
        real(8), intent(in) :: t_array(:)
        real(8), dimension(n, size(t_array), 2) :: pos_at_t
        real(8), dimension(n) :: x, y, t_current
        integer, dimension(n) :: ct  ! Current time index for each particle
        integer :: i, j, step
        real(8) :: phi, step_length, t_next, r, yn, xn
        integer, parameter :: max_steps = 100000
        
        ! Initialize
        x = 0.0d0
        y = 0.0d0
        t_current = 0.0d0
        pos_at_t = 0.0d0
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
                        pos_at_t(i, ct(i), 1) = intp(x(i), t_current(i), xn, t_next, t_array(ct(i)))
                        pos_at_t(i, ct(i), 2) = intp(y(i), t_current(i), yn, t_next, t_array(ct(i)))
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
    end function flight
    

end program pareto_sampler
