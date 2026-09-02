program test_radial_potential
    ! Standalone check of the new RADIAL_POTENTIAL radial normalization:
    ! confirms A(l) = R0^2/(l+1), i.e. RADIAL_SURFACE's R0^2/(l(l+1)) with the
    ! extra factor of l removed -- exactly the relationship requested when
    ! this option was added -- and that the LWS preset now carries it.
    use math_constants, only: prec
    use output_convention, only: radial_amplitude, RADIAL_SURFACE, RADIAL_POTENTIAL, LWS

    implicit none
    integer :: l
    real(prec) :: r0, A_surface, A_potential, ratio
    logical :: all_pass

    r0 = 6371.0e3_prec
    all_pass = .true.

    write(*,*) '=== test_radial_potential: A_potential(l) vs A_surface(l)*l ==='
    do l = 1, 5
        A_surface   = radial_amplitude(RADIAL_SURFACE,   l, r0)
        A_potential = radial_amplitude(RADIAL_POTENTIAL, l, r0)
        ratio = A_potential / (A_surface * real(l, prec))
        write(*,'(a,i2,a,es16.9,a,es16.9,a,f12.9)') &
            '  l=', l, '  A_potential=', A_potential, '  A_surface*l=', A_surface*real(l,prec), &
            '  ratio=', ratio
        if (abs(ratio - 1.0_prec) > 1.0e-12_prec) all_pass = .false.
        if (abs(A_potential - r0*r0/real(l+1,prec)) > 1.0e-6_prec) all_pass = .false.
    end do

    write(*,*)
    if (trim(LWS%radial_norm) == trim(RADIAL_POTENTIAL)) then
        write(*,*) 'LWS%radial_norm == RADIAL_POTENTIAL : PASS'
    else
        write(*,*) 'LWS%radial_norm == RADIAL_POTENTIAL : FAIL (got ', trim(LWS%radial_norm), ')'
        all_pass = .false.
    end if

    write(*,*)
    if (all_pass) then
        write(*,*) 'ALL CHECKS PASSED'
    else
        write(*,*) 'SOME CHECKS FAILED'
        stop 1
    end if

end program test_radial_potential
