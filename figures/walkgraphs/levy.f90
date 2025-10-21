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
    real(8), allocatable :: t_array(:), mad_array(:)
    integer :: n_particles, n_times, j

    call flight(3d0, 100, "f3.dat")
    call flight(3d0, 100, "w3.dat")
    call flight(1d0, 100, "f1.dat")
    call flight(1d0, 100, "w1.dat")

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
        x = 1.0d0 / (u**(1.0d0/(mu + 1.0d0)))
        
    end function p
       
    ! Levy walk simulation
    subroutine walk(mu, t_end, file_name) 
        implicit none
        real(8), intent(in) :: mu, t_end
        real(8) :: x, y, t_current
        integer :: step, io
        real(8) :: phi, step_length, t_next, r
        character(len=*), intent(in) :: file_name
        integer, parameter :: max_steps = 100000
        ! Initialize
        x = 0.0d0
        y = 0.0d0
        t_current = 0.0d0
        open(newunit=io, file=file_name) 

        ! Simulate each particle
        do step = 1, max_steps
            if(t_current >= t_end) exit
            ! Exit if we've passed the last time point
            step_length = p(mu)
            
            ! Sample random direction uniformly
            call random_number(phi)
            phi = 2.0d0 * PI * phi
            
            ! Calculate next time (ballistic motion: distance = time with unit velocity)
            t_next = t_current + step_length
            
            ! Update position after completing the step
            x = x + step_length * cos(phi)
            y = y + step_length * sin(phi)
            
            ! Calculate distance from origin
            r = sqrt(x**2 + y**2)
            
            ! Move to next time
            t_current = t_next
            write(io, *) x, y, t_current
        end do
        close(io) 
    end subroutine walk

    subroutine flight(mu, max_steps, file_name) 
        implicit none
        real(8) :: x, y, t_current
        real(8), intent(in) :: mu
        integer :: step, io
        real(8) :: phi, step_length, t_next, r
        integer, intent(in) :: max_steps
        character(len=*), intent(in) :: file_name
        
        ! Initialize
        x = 0.0d0
        y = 0.0d0
        t_current = 0.0d0
        open(newunit=io, file=file_name) 
        ! Simulate each particle
        do step = 1, max_steps
            ! Sample step length from distribution p(mu)
            step_length = p(mu)
            
            ! Sample random direction uniformly
            call random_number(phi)
            phi = 2.0d0 * PI * phi
            
            ! Calculate next time (flight motion: time step is constant)
            t_next = t_current + 1 
            
            ! Update position after completing the step
            x = x + step_length * cos(phi)
            y = y + step_length * sin(phi)
            
            ! Calculate distance from origin
            r = sqrt(x**2 + y**2)
            
            ! Move to next time
            t_current = t_next
            write(io, *) x, y, t_current
        end do
    close(io)
    end subroutine flight
    

end program pareto_sampler
