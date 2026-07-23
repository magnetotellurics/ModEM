program test_Tnr_uniform_sphere
    ! Minimal standalone driver to print Tnr/Tnsp AND Tnrp/Tns (l=1) at the
    ! surface of a uniform 100 Ohm.m sphere (conductivity all the way to
    ! r=0, no separate "core" layer -- unlike FWD1D.f90, which hardcodes an
    ! extra 1e5 S/m core layer at the bottom of whatever model is given to
    ! it; we deliberately avoid that here to match the closed-form
    ! uniform-sphere check).
    !
    ! IMPORTANT: FWD1D.f90 hardcodes primary_grid='E', which makes
    ! H%gridType=FACE in sourceField1d -- and the FACE branch uses Tnrp/Tns,
    ! NOT Tnr/Tnsp (confirmed directly in field1d.f90's FACE-branch loops):
    !     H%x, H%y  ->  Tnrp   (cell-centre Rr; NOT Tnsp)
    !     H%z, E%x, E%y  ->  Tns   (face Rs; NOT Tnr)
    ! The first version of this test only checked Tnr/Tnsp (the EDGE-branch
    ! quantities), which don't actually validate what the l=1,m=-1
    ! Ex/Ey/Hz-vs-Hx/Hy comparison exercises. This version prints all four.
    !
    ! Compare against the independent Python predictions in
    ! uniform_sphere_Tnr_predict.py (same closed-form formulas, cross-checked
    ! against two unambiguous physical limits). Since Rr=Rs=r0 exactly in
    ! this single-point test, Tnr should equal Tns and Tnsp should equal
    ! Tnrp at the level of the underlying physics (T(r0), T'(r0)) -- if they
    ! DON'T match each other here, that pins the discrepancy on the
    ! Tnr/Tnsp-vs-Tnrp/Tns code paths inside sourcePotential specifically,
    ! not on the closed-form physics.

    use field1d
    implicit none

    type(conf1d_t)                     :: earth
    complex(8), dimension(1,1)         :: Tnr, Tnsp, Tnrp, Tns   ! (1 radius, l=1 only)
    real(8), dimension(1)              :: Rr, Rs
    real(8)                            :: period
    integer                            :: lmax, istat

    lmax = 1
    period = 1000.0d0

    earth%r0   = 6.371d6
    earth%rmax = earth%r0 + 1000.0d3   ! arbitrary air layer thickness, unused for this test
    earth%tol  = 1.0d-9
    earth%tau  = 0.0d0                  ! no thin sheet / surface conductance

    ! SINGLE layer, uniform 100 Ohm.m all the way to r=0 -- matches
    ! demo_general_lm.py's current test config (sigma=1e-2 S/m).
    allocate(earth%layer(1), earth%sigma(1), STAT=istat)
    earth%layer(1) = 0.0d0
    earth%sigma(1) = 1.0d-2             ! 100 Ohm.m

    ! Evaluate at the surface (r0) for both grids. sourcePotential internally
    ! adds 1m to r0 for numerical safety (field1d.f90:527); negligible here.
    Rr(1) = earth%r0
    Rs(1) = earth%r0

    call sourcePotential(earth, lmax, period, Rr, Rs, Tnr, Tnsp, Tnrp, Tns)

    write(*,'(a)') 'Tnr(surface, l=1) = '
    write(*,'(a,es20.12,a,es20.12,a)') '  (', real(Tnr(1,1)), ', ', aimag(Tnr(1,1)), ')'
    write(*,'(a)') 'Tnsp(surface, l=1) = '
    write(*,'(a,es20.12,a,es20.12,a)') '  (', real(Tnsp(1,1)), ', ', aimag(Tnsp(1,1)), ')'
    write(*,'(a)') 'Tnrp(surface, l=1) = '
    write(*,'(a,es20.12,a,es20.12,a)') '  (', real(Tnrp(1,1)), ', ', aimag(Tnrp(1,1)), ')'
    write(*,'(a)') 'Tns(surface, l=1) = '
    write(*,'(a,es20.12,a,es20.12,a)') '  (', real(Tns(1,1)), ', ', aimag(Tns(1,1)), ')'

    deallocate(earth%layer, earth%sigma, STAT=istat)

end program test_Tnr_uniform_sphere
