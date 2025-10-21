program minimal_dlasrt_test
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none
    
    real(8) :: arr(5)
    integer(int32) :: n, info
    
    print *, "=== Minimal DLASRT Test ==="
    print *, ""
    
    ! Initialize array
    arr(1) = 5.0d0
    arr(2) = 2.0d0
    arr(3) = 8.0d0
    arr(4) = 1.0d0
    arr(5) = 6.0d0
    
    n = 5_int32
    info = 0_int32
    
    print *, "Array before:"
    print *, arr
    print *, ""
    print *, "n =", n
    print *, "info =", info
    print *, "size(arr) =", size(arr)
    print *, ""
    
    print *, "Calling DLASRT..."
    call DLASRT('I', n, arr, info)
    
    print *, ""
    print *, "After DLASRT:"
    print *, "info =", info
    
    if (info == 0) then
        print *, "Array after:"
        print *, arr
        print *, "✓ SUCCESS"
    else
        print *, "✗ FAILED with info =", info
    end if
    
end program minimal_dlasrt_test
