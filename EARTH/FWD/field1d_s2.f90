! *****************************************************************************
module field1d_s2
    ! Independent Fortran implementation of the layered-earth 1D EM induction
    ! solver derived in Sun & Egbert (2012), "A thin-sheet model for global
    ! electromagnetic induction", Geophys. J. Int. 189, 343-356, Section 2
    ! ("The homogeneous layered earth") + Appendix A.
    !
    ! Written independently of EARTH/FWD/field1d.f90 (S1 -- originally written
    ! using Kelbert (2006)'s conventions, later adjusted to match S2's for
    ! consistency; best published reference is Sun, Kelbert & Egbert (2015),
    ! J. Geophys. Res. Solid Earth, 120(10), 6771-6796) as a cross-check: no
    ! radial-solver or
    ! field-assembly code is shared between the two modules. Only the
    ! already-verified angular (Legendre / vector spherical harmonic) machinery
    ! is reused from field1d, whose Condon-Shortley phase convention has
    ! already been independently verified against sympy -- re-deriving that
    ! by hand a second time would only re-risk the kind of subtle phase bug
    ! this project has already spent a long time isolating and fixing.
    !
    ! TIME CONVENTION: e^{-i*omega*t}, matching the paper's own eq (1)-(2) and
    ! its k = sqrt(i*omega*mu0*sigma) (same formula/sign as field1d.f90's kl),
    ! so output is directly comparable to field1d.f90's output with no extra
    ! conjugation for time-convention reasons.
    !
    ! PHYSICS. Only the toroidal potential T is excited: the paper shows that
    ! under the quasi-static approximation the poloidal potential P must
    ! vanish identically in the source region (eq 21-22), so an external
    ! MT/ionospheric source only ever drives T -- exactly as in field1d.f90's
    ! Tnr-only approach. With P=0, eq (5)-(6) reduce to
    !
    !     E = -i*omega*mu0 * (r_hat/r) x grad_a(T)
    !     H = (1/r) grad_a(dT/dr) - (r_hat/r^2) laplacian_a(T)
    !
    ! Expanding with the standard identities r_hat x theta_hat = phi_hat,
    ! r_hat x phi_hat = -theta_hat, and laplacian_a(Y_l^m) = -l(l+1) Y_l^m
    ! (paper's own eq 9), and writing
    !     T(r,theta,phi) = sum_{l,m} coeff(l,m) * T_l(r) * Y_l^m(theta,phi)
    ! gives (independently re-derived twice -- directly, and via the
    ! T=R*Y/sqrt(l(l+1)) substitution used in pythonSolver -- and
    ! self-consistent both times):
    !
    !     E_theta = +(i*omega*mu0/r) * T_l(r)  * Yp(l,m)     Yp=(1/sin theta) dY/dphi
    !     E_phi   = -(i*omega*mu0/r) * T_l(r)  * Yt(l,m)     Yt=dY/dtheta
    !     E_r     = 0
    !     H_r     = +(l(l+1)/r^2)    * T_l(r)  * Y(l,m)
    !     H_theta = +(1/r)           * T_l'(r) * Yt(l,m)
    !     H_phi   = +(1/r)           * T_l'(r) * Yp(l,m)
    !
    ! No R0^2 or l(l+1)-division normalization is applied anywhere here (unlike
    ! field1d.f90's formulas) -- these are the literal, dimensionally complete
    ! eq (5)-(6) in SI units (r in metres, mu0 in H/m, omega in rad/s), nothing
    ! more is needed. Likewise, no final conjg() is applied to the assembled
    ! field, unlike field1d.f90's sourceField1d: that conjugation exists there
    ! specifically to compensate for how THEIR m,-m conjugate-pairing loop is
    ! built; here each Y_l^{-m} is obtained directly and independently via the
    ! standard identity Y_l^{-m}=(-1)^m*conjg(Y_l^m) (see Yval() below) and
    ! summed directly over the full m=-l..l range, which is already a
    ! standard, complete complex spherical-harmonic reconstruction with no
    ! further global conjugation required.
    !
    ! RADIAL SOLVE. T_l(r) satisfies the Helmholtz equation (paper's eq 7),
    ! general solution (eq 8) alpha*psi_l(kr) + beta*xi_l(kr) in each
    ! homogeneous layer, psi_l, xi_l the Riccati-Bessel functions of the first
    ! and third kind. These are evaluated via the rescaled / continued-
    ! fraction (Lentz, 1976, cited by the paper) scheme of Appendix A
    ! (eq A1-A10), needed for numerical stability at complex argument.
    ! Adjacent layers are connected via continuity of T, T' (eq 12) -- used
    ! here directly, in place of literally transcribing the explicit 2x2
    ! transfer matrix of eq (13)-(14), which is algebraically equivalent but a
    ! considerably more error-prone transcription target (many primes and
    ! subscripts); continuity of (T,T') plus a per-layer Wronskian projection
    ! onto that layer's own (psi_l, xi_l) basis is the same physics, more
    ! directly implemented. The thin-sheet jump (eq 17) is applied at the
    ! surface, and the (alpha_l^T, beta_l^T) split in the air is obtained by
    ! matching (eq 23). alpha_l^T is exactly the paper's own "external
    ! multipole moment" of the source (text following eq 23) -- the whole
    ! radial profile is normalized here to alpha_l^T=1 ("per unit external
    ! field", the same convention field1d.f90's sourcePotential uses), and
    ! sourceField1d_s2 multiplies by the user's coeff(l,m) afterward.
    !
    ! Numerically this reuses the same two robustness techniques already
    ! proven in field1d.f90's sourcePotential -- both are standard, paper-
    ! cited numerical-analysis techniques (Appendix A's rescaling; Lentz's
    ! continued fraction), not anything specific to field1d.f90's own
    ! derivation: per-layer sub-stepping to avoid exp() overflow in thick/
    ! conductive layers, and ratio-based (not T'/T) evaluation of the core
    ! impedance to avoid 0/0 at high l.

    use field1d, only: conf1d_t, legendre_norm, vsharm
    use griddef
    use sg_vector
    use utilities
    use math_constants
    implicit none

    private

    public :: sourcePotential_s2, sourceField1d_s2, surfaceField1d_s2

contains

    !===========================================================================
    ! Riccati-Bessel functions, rescaled per Appendix A, complex argument
    !===========================================================================

    subroutine riccati_ratio_psi(lmax, z, Rl)
        ! R_l^psi(z) = psi_l(z)/psi_{l-1}(z), l=1..lmax.
        !
        ! An earlier version of this routine seeded the l=lmax ratio via the
        ! alternating-sign continued-fraction form of eq A4 (Lentz, 1976) as
        ! literally transcribed from the paper's (OCR'd) Appendix A text.
        ! That transcription turned out to be wrong -- confirmed two ways:
        ! against an explicit closed-form psi_1(z)/psi_0(z), psi_2(z)/psi_1(z)
        ! check (sin/cos-based, unambiguous for complex z), and against this
        ! routine's own downward-recursion result. Rather than re-risk a
        ! second mistranscription of the same dense equation, this routine
        ! instead uses the more basic (if less elegant) Miller-type downward
        ! recursion: the EXACT 3-term relation R_l = 1/[(2l+1)/z - R_{l+1}]
        ! (paper's eq A3, itself independently re-derived from the standard
        ! spherical Bessel recursion j_{l+1}=(2l+1)/z*j_l-j_{l-1} and
        ! confirmed correct) is seeded with R_{l+1}=0 at a degree well above
        ! lmax (where the true ratio is negligibly small) and recursed
        ! downward -- downward recursion of this relation is numerically
        ! stable, unlike the forward direction. The seed degree is doubled
        ! until the l=1..lmax results stop changing, so no z- or l-dependent
        ! tuning is required.
        integer, intent(in)                      :: lmax
        complex(8), intent(in)                   :: z
        complex(8), dimension(:), intent(inout)  :: Rl
        complex(8), dimension(lmax) :: Rl_prev
        complex(8) :: Rnext, Rcur
        integer    :: l, Ltop, buffer, buffer_cap
        logical    :: converged

        buffer = max(50, 4*lmax + 4*ceiling(abs(z)))
        ! The FIRST pass below can never register as "converged" (Rl_prev
        ! starts at the 1e30 sentinel, guaranteeing a mismatch on purpose,
        ! so there are always >=2 real estimates to compare) -- so the loop
        ! always doubles buffer at least once even when the initial buffer
        ! was already more than sufficient. A FIXED cap (e.g. 200000) is
        ! therefore only safe when |z| is small enough that the initial
        ! buffer itself is well under the cap; for a large |z| (e.g. the
        ! Earth's core: earth%sigma for the core layer is hardcoded to 1e5
        ! S/m in FWD1D.f90, so |kl*r_core| ~ 1e5 at ordinary MT periods),
        ! the initial buffer can ALREADY exceed a fixed 200000 cap, so the
        ! very first mandatory doubling trips the "no convergence" warning
        ! even though the Miller recursion itself converges just fine --
        ! found 2026-07-27, `riccati_ratio_psi` warned at every call for a
        ! realistic layered Earth model (layered.prm) at T=1000s. Scale the
        ! cap off the INITIAL buffer instead, so at least a few real
        ! doublings are always possible regardless of |z|; this subroutine
        ! is called only once or twice per period (not per grid point,
        ! see sourcePotential_s2), so the extra iterations for
        ! large |z| are computationally cheap.
        buffer_cap = max(200000, 8*buffer)
        Rl_prev(:) = dcmplx(1.0d30, 0.0d0)

        do
            Ltop = lmax + buffer
            Rnext = dcmplx(0.0d0, 0.0d0)
            do l = Ltop,1,-1
                Rcur = 1.0d0/(dcmplx(dble(2*l+1),0.0d0)/z - Rnext)
                if (l <= lmax) Rl(l) = Rcur
                Rnext = Rcur
            end do

            converged = .true.
            do l = 1,lmax
                if (abs(Rl(l) - Rl_prev(l)) > 1.0d-13*max(1.0d0, abs(Rl(l)))) converged = .false.
            end do
            if (converged) exit

            Rl_prev(:) = Rl(:)
            buffer = buffer*2
            if (buffer > buffer_cap) then
                write(0,*) 'WARNING field1d_s2/riccati_ratio_psi: no convergence at z=',z
                exit
            end if
        end do

    end subroutine riccati_ratio_psi

    subroutine riccati_rescaled(lmax, z0, z, kind, F, Fp)
        ! kind=1: rescaled psi_bar_l(z|z0), psi_bar_l'(z|z0)  [eq A1,A5,A6,A9]
        ! kind=2: rescaled  xi_bar_l(z|z0),  xi_bar_l'(z|z0)  [eq A2,A7,A8,A10]
        ! F,Fp are the value and its z-derivative, for l=1..lmax.
        integer, intent(in)                     :: lmax, kind
        complex(8), intent(in)                  :: z0, z
        complex(8), dimension(:), intent(inout) :: F, Fp
        complex(8) :: F0, R1
        complex(8), dimension(lmax) :: Rl
        complex(8), parameter :: ci = dcmplx(0.0d0,1.0d0)
        integer :: l

        F(:) = dcmplx(0.0d0,0.0d0)
        Fp(:) = dcmplx(0.0d0,0.0d0)

        select case (kind)
        case (1)
            ! psi_bar_0(z|z0), eq A9
            F0 = (exp(ci*z - aimag(z0)) - exp(-ci*z - aimag(z0))) / (2.0d0*ci)
            call riccati_ratio_psi(lmax, z, Rl)
            do l = 1,lmax
                R1 = dcmplx(dble(2*l+1),0.0d0)/z0 * F0
                F0 = R1 * Rl(l)
                F(l) = F0
                Fp(l) = R1 - dcmplx(dble(l),0.0d0)*F0/z
            end do
        case (2)
            ! xi_bar_0(z|z0), eq A10. R_l^xi(z) via the FORWARD recursion
            ! (stable in this direction, unlike R_l^psi) from R_1^xi(z)=1/z-i
            ! (paper's eq A3, text following eq A4). NOTE the prefactor here
            ! is z0/(2l+1) (eq A7-A8) -- the RECIPROCAL of the psi case's
            ! (2l+1)/z0 (eq A5-A6) -- easy to invert by mistake, so flagged.
            F0 = -ci * exp(ci*z + aimag(z0))
            R1 = (1.0d0/z) - ci
            do l = 1,lmax
                Fp(l) = z0/dcmplx(dble(2*l+1),0.0d0) * F0   ! = z0/(2l+1) * xi_bar_{l-1}
                F0 = Fp(l) * R1
                F(l) = F0
                Fp(l) = Fp(l) - dcmplx(dble(l),0.0d0)*F0/z
                R1 = dcmplx(dble(2*l+1),0.0d0)/z - 1.0d0/R1
            end do
        end select

    end subroutine riccati_rescaled

    subroutine propagate_T(lmax, z0, T0, Tp0, z, T, Tp)
        ! Propagate the (T, dT/dz) pair from z0 to z within a single
        ! homogeneous layer (constant k), via Wronskian projection onto the
        ! layer's own (psi,xi) basis at z0 -- equivalent to, and used here in
        ! place of, forming the explicit 2x2 transfer matrix of eq (13)-(14).
        integer, intent(in)                    :: lmax
        complex(8), intent(in)                 :: z0, z
        complex(8), dimension(lmax), intent(in):: T0, Tp0
        complex(8), dimension(lmax), intent(inout) :: T, Tp
        complex(8), dimension(lmax) :: S0,Sp0,C0,Cp0,Sz,Spz,Cz,Cpz,A,B,T0h,Tp0h
        complex(8), parameter :: wr = dcmplx(0.0d0,-1.0d0)
        integer :: l

        call riccati_rescaled(lmax, z0, z0, 1, S0, Sp0)
        call riccati_rescaled(lmax, z0, z0, 2, C0, Cp0)

        do l = 1,lmax
            A(l) = (Cp0(l)*T0(l) - C0(l)*Tp0(l)) * wr
            T0h(l) = T0(l) - A(l)*S0(l)
            Tp0h(l) = Tp0(l) - A(l)*Sp0(l)
            B(l) = (S0(l)*Tp0h(l) - Sp0(l)*T0h(l)) * wr
        end do

        call riccati_rescaled(lmax, z0, z, 1, Sz, Spz)
        call riccati_rescaled(lmax, z0, z, 2, Cz, Cpz)

        do l = 1,lmax
            T(l) = Sz(l)*A(l) + Cz(l)*B(l)
            Tp(l) = Spz(l)*A(l) + Cpz(l)*B(l)
        end do

    end subroutine propagate_T

    !===========================================================================
    ! Radial solve: source potential T_l(r), T_l'(r), normalized to
    ! alpha_l^T = 1 (per unit external multipole amplitude)
    !===========================================================================

    subroutine sourcePotential_s2(earth, lmax, period, Rr, Rs, Tnr, Tnsp, Tnrp, Tns)

        type(conf1d_t), intent(in)          :: earth
        integer, intent(in)                 :: lmax
        real(8), intent(in)                 :: period
        real(8), dimension(:), intent(in)   :: Rr, Rs
        complex(8), dimension(:,:), intent(inout)           :: Tnr, Tnsp
        complex(8), dimension(:,:), intent(inout), optional :: Tnrp, Tns
        ! local
        complex(8), dimension(:), allocatable :: kl
        real(8), dimension(:), allocatable    :: rl
        complex(8), dimension(lmax) :: tn, tnp, Tval, Tpval, alpha_raw, beta_norm
        complex(8), dimension(lmax) :: rn0, rnp0, phn0, phnp0
        real(8) :: omega, rmax, overflow_arg, dr_sub, r_sub_bot, r_sub_top
        integer :: Nl, idr, idrmin, idrmax, ids, idsmin, idsmax
        integer :: j, l, nsub, isub, istat

        if (lmax <= 0) then
            write(0,*) 'Error in sourcePotential_s2: no potentials for degree 0'
            return
        end if

        Nl = size(earth%layer)
        allocate(rl(Nl), kl(Nl), STAT=istat)
        omega = 2.0d0*PI/period
        do j = 1,Nl
            rl(j) = earth%r0 + 1.0d0 - earth%layer(j)
            kl(j) = sqrt(dcmplx(0.0d0,1.0d0)*omega*MU_0*earth%sigma(j))
        end do
        rmax = earth%rmax + 1.0d0

        Tnr(:,:) = dcmplx(0.0d0,0.0d0)
        Tnsp(:,:) = dcmplx(0.0d0,0.0d0)
        if (present(Tnrp)) Tnrp(:,:) = dcmplx(0.0d0,0.0d0)
        if (present(Tns))  Tns(:,:)  = dcmplx(0.0d0,0.0d0)

        !--------------------------------------------------------------
        ! Core (layer Nl): regular solution T_l(r) = psi_bar_l(kl*r | kl*rl),
        ! self-referenced rescaling keeps the starting value O(1).
        !--------------------------------------------------------------
        call find_index(Rr, 0.0d0, rl(Nl), idrmin, idrmax)
        if ((idrmin > 0) .and. (idrmax > 0)) then
            do idr = idrmin,idrmax
                call riccati_rescaled(lmax, kl(Nl)*rl(Nl), kl(Nl)*Rr(idr), 1, Tval, Tpval)
                Tnr(idr,:) = Tval(:)
                if (present(Tnrp)) Tnrp(idr,:) = kl(Nl)*Tpval(:)
            end do
        end if

        call find_index(Rs, 0.0d0, rl(Nl), idsmin, idsmax)
        if ((idsmin > 0) .and. (idsmax > 0)) then
            do ids = idsmin,idsmax
                call riccati_rescaled(lmax, kl(Nl)*rl(Nl), kl(Nl)*Rs(ids), 1, Tval, Tpval)
                Tnsp(ids,:) = kl(Nl)*Tpval(:)
                if (present(Tns)) Tns(ids,:) = Tval(:)
            end do
        end if

        call riccati_rescaled(lmax, kl(Nl)*rl(Nl), kl(Nl)*rl(Nl), 1, tn, tnp)

        ! Robust core impedance: avoid 0/0 for l >> |kl(Nl)*rl(Nl)| (tn
        ! underflows) by computing kl*(1/R_l^psi - l/z0) directly, R_l^psi
        ! staying O(1) throughout (same trick documented for field1d.f90).
        block
            complex(8), dimension(lmax) :: Rl_core
            call riccati_ratio_psi(lmax, kl(Nl)*rl(Nl), Rl_core)
            do l = 1,lmax
                tnp(l) = kl(Nl)*(1.0d0/Rl_core(l) - dcmplx(dble(l),0.0d0)/(kl(Nl)*rl(Nl)))
            end do
        end block

        do l = 1,lmax
            if (abs(tn(l)) > 0.0d0) then
                Tnr(:,l) = Tnr(:,l)/tn(l)
                Tnsp(:,l) = Tnsp(:,l)/tn(l)
                if (present(Tnrp)) Tnrp(:,l) = Tnrp(:,l)/tn(l)
                if (present(Tns))  Tns(:,l)  = Tns(:,l)/tn(l)
            else
                Tnr(:,l) = dcmplx(0.0d0,0.0d0)
                Tnsp(:,l) = dcmplx(0.0d0,0.0d0)
                if (present(Tnrp)) Tnrp(:,l) = dcmplx(0.0d0,0.0d0)
                if (present(Tns))  Tns(:,l)  = dcmplx(0.0d0,0.0d0)
            end if
        end do

        !--------------------------------------------------------------
        ! Outward, layer by layer, using continuity of (T,T') at each
        ! interface (eq 12); sub-stepped within thick/conductive layers to
        ! avoid exp() overflow (same criterion as field1d.f90).
        !--------------------------------------------------------------
        do j = Nl-1,1,-1

            overflow_arg = abs(aimag(kl(j))) * (rl(j) - rl(j+1))
            nsub = max(1, ceiling(overflow_arg/200.0d0))
            dr_sub = (rl(j) - rl(j+1)) / nsub

            phn0(:) = dcmplx(1.0d0,0.0d0)
            phnp0(:) = tnp(:)/kl(j)
            r_sub_bot = rl(j+1)

            do isub = 1,nsub
                r_sub_top = rl(j+1) + isub*dr_sub

                call find_index(Rr, r_sub_bot, r_sub_top, idrmin, idrmax)
                if ((idrmin > 0) .and. (idrmax > 0)) then
                    do idr = idrmin,idrmax
                        call propagate_T(lmax, kl(j)*r_sub_bot, phn0, phnp0, kl(j)*Rr(idr), Tval, Tpval)
                        Tnr(idr,:) = Tval(:)
                        if (present(Tnrp)) Tnrp(idr,:) = kl(j)*Tpval(:)
                    end do
                end if

                call find_index(Rs, r_sub_bot, r_sub_top, idsmin, idsmax)
                if ((idsmin > 0) .and. (idsmax > 0)) then
                    do ids = idsmin,idsmax
                        call propagate_T(lmax, kl(j)*r_sub_bot, phn0, phnp0, kl(j)*Rs(ids), Tval, Tpval)
                        Tnsp(ids,:) = kl(j)*Tpval(:)
                        if (present(Tns)) Tns(ids,:) = Tval(:)
                    end do
                end if

                call propagate_T(lmax, kl(j)*r_sub_bot, phn0, phnp0, kl(j)*r_sub_top, tn, tnp)

                do l = 1,lmax
                    tnp(l) = kl(j)*tnp(l)/tn(l)
                end do
                do l = 1,lmax
                    if (abs(tn(l)) > 0.0d0) then
                        Tnr(:,l) = Tnr(:,l)/tn(l)
                        Tnsp(:,l) = Tnsp(:,l)/tn(l)
                        if (present(Tnrp)) Tnrp(:,l) = Tnrp(:,l)/tn(l)
                        if (present(Tns))  Tns(:,l)  = Tns(:,l)/tn(l)
                    end if
                end do

                phn0(:) = dcmplx(1.0d0,0.0d0)
                phnp0(:) = tnp(:)/kl(j)
                r_sub_bot = r_sub_top

            end do ! isub
        end do ! layers

        !--------------------------------------------------------------
        ! Thin-sheet jump at the surface (eq 17): T continuous, T' jumps.
        !--------------------------------------------------------------
        do l = 1,lmax
            rn0(l) = dcmplx(1.0d0,0.0d0)
            rnp0(l) = tnp(l) - dcmplx(0.0d0,1.0d0)*omega*MU_0*earth%tau*rn0(l)
        end do

        !--------------------------------------------------------------
        ! Air-layer matching (eq 23): [alpha_l^T; beta_l^T] from (T,T') at r1+.
        !
        ! NOTE: eq (23) as literally written involves r1^(-l-1), r1^l,
        ! r1^(l+1) directly, which overflows complex(8) range for even
        ! moderate l at Earth-radius scale (r1 ~ 6.4e6 m -- e.g. r1^15 is
        ! already ~1e97). field1d.f90's airprop avoids this by working with
        ! the DIMENSIONLESS ratio r/r0 throughout; re-deriving eq (23) in
        ! terms of Ahat = alpha_l^T * r1^(l+1), Bhat = beta_l^T * r1^(-l)
        ! (i.e. the coefficients of x^(l+1), x^(-l) with x=r/r1) removes all
        ! high powers of r1, leaving only a single factor of r1^1:
        !     Ahat = (l*T1 + r1*Tp1)/(2l+1)
        !     Bhat = ((l+1)*T1 - r1*Tp1)/(2l+1)
        ! This is algebraically identical to eq (23) (confirmed: matches
        ! airprop's A1,A2 exactly) and stays O(1) for any l.
        !--------------------------------------------------------------
        do l = 1,lmax
            alpha_raw(l) = (dble(l)*rn0(l) + rl(1)*rnp0(l)) / dble(2*l+1)
            beta_norm(l) = (dble(l+1)*rn0(l) - rl(1)*rnp0(l)) / dble(2*l+1)
        end do

        ! IMPORTANT: normalizing Tnr(:,l)/=alpha_raw(l) alone only sets the
        ! coefficient of x^(l+1) (x=r/r1) to 1 -- NOT the coefficient of the
        ! physical r^(l+1), which is what "alpha_l^T=1" is supposed to mean
        ! (T_l(r)=r^(l+1)+beta*r^-l with UNIT coefficient on r^(l+1), matching
        ! how coeff(l,m) is meant to represent the actual external source
        ! amplitude, same convention as field1d.f90/pythonSolver's K0). Since
        ! x^(l+1)=r^(l+1)/r1^(l+1), the coefficient-of-x^(l+1)=1 basis and the
        ! coefficient-of-r^(l+1)=1 basis differ by exactly r1^(l+1) (found by
        ! direct comparison against reference_unit_sphere_s2.py, 2026-07-23)
        ! -- dividing by alpha_raw(l)/rl(1)**(l+1) (instead of alpha_raw(l)
        ! alone) restores the intended convention. This factor is applied
        ! once, multiplicatively, here -- unlike the raw r1^(l+1) terms
        ! removed from the eq(23)/eq(18) algebra above, it cannot compound
        ! into overflow, though for lmax high enough that r1**(lmax+1)
        ! itself approaches complex(8) range (~1e308; irrelevant for
        ! realistic MT lmax) it would eventually need sub-stepping too.
        do l = 1,lmax
            if (abs(alpha_raw(l)) > 0.0d0) then
                block
                    complex(8) :: norm_factor
                    norm_factor = alpha_raw(l) / rl(1)**(l+1)
                    Tnr(:,l) = Tnr(:,l)/norm_factor
                    Tnsp(:,l) = Tnsp(:,l)/norm_factor
                    if (present(Tnrp)) Tnrp(:,l) = Tnrp(:,l)/norm_factor
                    if (present(Tns))  Tns(:,l)  = Tns(:,l)/norm_factor
                end block
                beta_norm(l) = beta_norm(l)/alpha_raw(l)
            else
                Tnr(:,l) = dcmplx(0.0d0,0.0d0)
                Tnsp(:,l) = dcmplx(0.0d0,0.0d0)
                if (present(Tnrp)) Tnrp(:,l) = dcmplx(0.0d0,0.0d0)
                if (present(Tns))  Tns(:,l)  = dcmplx(0.0d0,0.0d0)
                beta_norm(l) = dcmplx(0.0d0,0.0d0)
            end if
        end do

        !--------------------------------------------------------------
        ! Air region (eq 18, quasi-static limit), same r/r1-relative
        ! rewrite: with x=r/r1, T_l(r) = r1^(l+1)*[x^(l+1) + beta_norm*x^(-l)]
        ! (the r1^(l+1) restores the coefficient-of-r^(l+1)=1 convention, see
        ! note above), and, by the chain rule (dx/dr = 1/r1),
        ! T_l'(r) = r1^l*[(l+1)*x^l - l*beta_norm*x^-(l+1)].
        !--------------------------------------------------------------
        call find_index(Rr, rl(1), rmax, idrmin, idrmax)
        if ((idrmin > 0) .and. (idrmax > 0)) then
            do idr = idrmin,idrmax
                block
                    real(8) :: x
                    x = Rr(idr)/rl(1)
                    do l = 1,lmax
                        Tnr(idr,l) = rl(1)**(l+1) * (x**(l+1) + beta_norm(l)*x**(-l))
                        if (present(Tnrp)) Tnrp(idr,l) = rl(1)**l * (dble(l+1)*x**l &
                                                        - dble(l)*beta_norm(l)*x**(-l-1))
                    end do
                end block
            end do
        end if

        call find_index(Rs, rl(1), rmax, idsmin, idsmax)
        if ((idsmin > 0) .and. (idsmax > 0)) then
            do ids = idsmin,idsmax
                block
                    real(8) :: x
                    x = Rs(ids)/rl(1)
                    do l = 1,lmax
                        Tnsp(ids,l) = rl(1)**l * (dble(l+1)*x**l - dble(l)*beta_norm(l)*x**(-l-1))
                        if (present(Tns)) Tns(ids,l) = rl(1)**(l+1) * (x**(l+1) + beta_norm(l)*x**(-l))
                    end do
                end block
            end do
        end if

        deallocate(rl, kl, STAT=istat)

    end subroutine sourcePotential_s2

    !===========================================================================
    ! Angular helper: Y_l^m (or Yt_l^m, Yp_l^m) for negative m, via the
    ! identity Y_l^{-m} = (-1)^m * conjg(Y_l^m), applied to whichever of
    ! vsharm's Y/Yt/Yp arrays is passed in (all three follow the same rule,
    ! independently verified -- see plan/commit notes).
    !===========================================================================

    complex(8) function Yval(F, l, m, lmax)
        integer, intent(in) :: l, m, lmax
        complex(8), dimension(lmax,lmax+1), intent(in) :: F

        if (m >= 0) then
            Yval = F(l, m+1)
        else
            Yval = ((-1.0d0)**(-m)) * conjg(F(l, -m+1))
        end if

    end function Yval

    !===========================================================================
    ! Field assembly
    !===========================================================================

    subroutine sourceField1d_s2(earth, lmax, coeff, period, grid, H, E)

        type(conf1d_t), intent(in)              :: earth
        integer, intent(in)                     :: lmax
        complex(8), dimension(:), intent(in)    :: coeff
        real(8), intent(in)                     :: period
        type(grid_t), intent(in)                :: grid
        type(cvector), intent(inout)            :: H
        type(cvector), intent(inout), optional  :: E
        ! local
        real(8), dimension(:), allocatable       :: Rr, Rs
        real(8), dimension(lmax+1,lmax+1)        :: P_lm
        complex(8), dimension(lmax,lmax+1)       :: Y, Yt, Yp
        complex(8), dimension(:,:), allocatable  :: Tnr, Tnsp, Tnrp, Tns
        complex(8), dimension(:,:), allocatable  :: coeff2
        real(8)    :: omega
        complex(8) :: iomu0, csum
        integer    :: Np, Nt, Nr, Nrr, Nrs, ncoeff
        integer    :: i, j, k, l, m, istat
        integer    :: j1, j2
        real(8), parameter :: pole_tol = 1.0d-8 ! radians; true global poles only, not regional-grid edges

        ncoeff = 0
        do l = 0,lmax
            ncoeff = ncoeff + (2*l+1)
        end do
        if (size(coeff) /= ncoeff) then
            write(0,*) 'Error in sourceField1d_s2: bad sph. harm. coeffs vector (size ', &
                       size(coeff),'); must be size ',ncoeff
        end if

        ! Unpack the driver's flat, degree-blocked coeff(:) array (block for
        ! degree l starts at index l^2+1, size 2l+1, ordered m=0,1,-1,2,-2,...
        ! -- see FWD1D.f90/field1d.f90's read_modelParam pipeline) into a
        ! plain coeff2(l,m), m=-l..l array for direct summation below.
        allocate(coeff2(lmax,-lmax:lmax), STAT=istat)
        coeff2(:,:) = dcmplx(0.0d0,0.0d0)
        do l = 1,lmax
            coeff2(l,0) = coeff(l*l+1)
            do m = 1,l
                coeff2(l,m)  = coeff(l*l+2*m)
                coeff2(l,-m) = coeff(l*l+2*m+1)
            end do
        end do

        Np = grid%nx
        Nt = grid%ny
        Nr = grid%nz
        Nrs = Nr+1
        Nrr = Nr

        allocate(Rr(Nr), Rs(Nr+1), STAT=istat)
        Rs(1:Nr+1) = 1.0d3 * grid%r(1:Nr+1)
        Rr(1:Nr) = Rs(1:Nr) - 1.0d3 * grid%dr(1:Nr)/2.0d0

        allocate(Tnr(Nrr,lmax), Tnsp(Nrs,lmax), STAT=istat)
        if (H%gridType == FACE) then
            allocate(Tnrp(Nrr,lmax), Tns(Nrs,lmax), STAT=istat)
            call sourcePotential_s2(earth, lmax, period, Rr, Rs, Tnr, Tnsp, Tnrp, Tns)
        else
            call sourcePotential_s2(earth, lmax, period, Rr, Rs, Tnr, Tnsp)
        end if

        if (any(Tnr /= Tnr) .or. any(Tnsp /= Tnsp)) then
            write(0,*) 'Error in sourceField1d_s2: NaN in source potentials at T=',period,' s'
            deallocate(Tnr, Tnsp, STAT=istat)
            if (H%gridType == FACE) deallocate(Tnrp, Tns, STAT=istat)
            deallocate(Rr, Rs, coeff2, STAT=istat)
            return
        end if

        omega = 2.0d0*PI/period
        iomu0 = dcmplx(0.0d0,1.0d0) * omega * MU_0

        call zero_cvector(H)
        if (present(E)) call zero_cvector(E)

        if (H%gridType == EDGE) then

        ! EDGE H staggering (primary_grid='H'): H on primary edges, E on
        ! primary faces -- same grid pattern as field1d.f90's EDGE branch:
        !   H%x (phi,   EDGE x): node_th x mid_ph  x Rs -- Yp, Tnsp
        !   H%y (theta, EDGE y): mid_th  x node_ph x Rs -- Yt, Tnsp
        !   H%z (r,     EDGE z): node_th x node_ph x Rr -- Y,  Tnr
        !   E%y (theta, FACE y): node_th x mid_ph  x Rr -- Yp, Tnr   (if present)
        !   E%x (phi,   FACE x): mid_th  x node_ph x Rr -- Yt, Tnr   (if present)

        ! H%x -- phi component: H_phi = +(1/r) * T'(r) * Yp(l,m)
        ! node theta: when H%gridType==EDGE, H%x is undefined at BOTH the
        ! North pole (j=1, theta=0) and the South pole (j=Nt+1, theta=pi)
        ! -- skip each endpoint ONLY if this grid actually reaches it.
        ! The below logic allows for global vs regional grids;
        ! it's good practice to define regional grids away from the poles - for
        ! polar regions, will need to rotate to fake pole (supported)
        j1 = 1
        j2 = Nt+1
        if (grid%th(1) < pole_tol) j1 = 2
        if (grid%th(Nt+1) > PI - pole_tol) j2 = Nt
        do j = j1,j2
            call legendre_norm(lmax, cos(grid%th(j)), P_lm)
            do i = 1,Np
                call vsharm(lmax, cos(grid%th(j)), grid%ph(i)+grid%dp(i)/2, P_lm, Y, Yt, Yp)
                do k = 1,Nrs
                    csum = dcmplx(0.0d0,0.0d0)
                    do l = 1,lmax
                        do m = -l,l
                            csum = csum + coeff2(l,m)*Tnsp(k,l)/Rs(k)*Yval(Yp,l,m,lmax)
                        end do
                    end do
                    H%x(i,j,k) = csum
                end do
            end do
        end do

        ! H%y -- theta component: H_theta = +(1/r) * T'(r) * Yt(l,m)
        do j = 1,Nt
            call legendre_norm(lmax, cos(grid%th(j)+grid%dt(j)/2), P_lm)
            do i = 1,Np+1
                call vsharm(lmax, cos(grid%th(j)+grid%dt(j)/2), grid%ph(i), P_lm, Y, Yt)
                do k = 1,Nrs
                    csum = dcmplx(0.0d0,0.0d0)
                    do l = 1,lmax
                        do m = -l,l
                            csum = csum + coeff2(l,m)*Tnsp(k,l)/Rs(k)*Yval(Yt,l,m,lmax)
                        end do
                    end do
                    H%y(i,j,k) = csum
                end do
            end do
        end do

        ! H%z -- radial component: H_r = +(l(l+1)/r^2) * T(r) * Y(l,m)
        do j = 1,Nt+1
            call legendre_norm(lmax, cos(grid%th(j)), P_lm)
            do i = 1,Np+1
                call vsharm(lmax, cos(grid%th(j)), grid%ph(i), P_lm, Y)
                do k = 1,Nrr
                    csum = dcmplx(0.0d0,0.0d0)
                    do l = 1,lmax
                        do m = -l,l
                            csum = csum + coeff2(l,m)*Tnr(k,l)*dble(l*(l+1))/Rr(k)**2*Yval(Y,l,m,lmax)
                        end do
                    end do
                    H%z(i,j,k) = csum
                end do
            end do
        end do

        if (present(E)) then
        ! E%y -- theta component: E_theta = +(i*omega*mu0/r) * T(r) * Yp(l,m)
        ! node theta: when H%gridType==EDGE (so E%gridType==FACE), E%y is
        ! undefined at BOTH poles -- j1/j2 computed above for H%x, same
        ! node-theta range. The below logic allows for global vs regional grids;
        ! it's good practice to define regional grids away from the poles - for
        ! polar regions, will need to rotate to fake pole (supported)
        do j = j1,j2
            call legendre_norm(lmax, cos(grid%th(j)), P_lm)
            do i = 1,Np
                call vsharm(lmax, cos(grid%th(j)), grid%ph(i)+grid%dp(i)/2, P_lm, Y, Yt, Yp)
                do k = 1,Nrr
                    csum = dcmplx(0.0d0,0.0d0)
                    do l = 1,lmax
                        do m = -l,l
                            csum = csum + coeff2(l,m)*iomu0*Tnr(k,l)/Rr(k)*Yval(Yp,l,m,lmax)
                        end do
                    end do
                    E%y(i,j,k) = csum
                end do
            end do
        end do

        ! E%x -- phi component: E_phi = -(i*omega*mu0/r) * T(r) * Yt(l,m)
        ! node phi (i=1..Np+1) -- FIX (2026-07-26): loop previously stopped at
        ! i=Np, leaving E%x(Np+1,:,:) at its zero_cvector default; see the
        ! identical fix + rationale in field1d.f90's corresponding E%x block.
        do j = 1,Nt
            call legendre_norm(lmax, cos(grid%th(j)+grid%dt(j)/2), P_lm)
            do i = 1,Np+1
                call vsharm(lmax, cos(grid%th(j)+grid%dt(j)/2), grid%ph(i), P_lm, Y, Yt)
                do k = 1,Nrr
                    csum = dcmplx(0.0d0,0.0d0)
                    do l = 1,lmax
                        do m = -l,l
                            csum = csum - coeff2(l,m)*iomu0*Tnr(k,l)/Rr(k)*Yval(Yt,l,m,lmax)
                        end do
                    end do
                    E%x(i,j,k) = csum
                end do
            end do
        end do
        end if ! present(E)

        else if (H%gridType == FACE) then

        ! FACE H staggering (primary_grid='E'): H on primary faces, E on
        ! primary edges -- mirrors the EDGE staggering with x<->y swapped:
        !   H%y (theta, FACE y): node_th x mid_ph  x Rr -- Yt, Tnrp
        !   H%x (phi,   FACE x): mid_th  x node_ph x Rr -- Yp, Tnrp
        !   H%z (r,     FACE z): mid_th  x mid_ph  x Rs -- Y,  Tns
        !   E%x (phi,   EDGE x): node_th x mid_ph  x Rs -- Yt, Tns   (if present)
        !   E%y (theta, EDGE y): mid_th  x node_ph x Rs -- Yp, Tns   (if present)

        ! Node-theta components below (H%y on FACE, E%x on EDGE) use Yt and
        ! are undefined at BOTH the North pole (j=1, theta=0) and the South
        ! pole (j=Nt+1, theta=pi) -- skip each endpoint ONLY if this grid
        ! actually reaches it. The below logic allows for global vs regional grids;
        ! it's good practice to define regional grids away from the poles - for
        ! polar regions, will need to rotate to fake pole (supported)
        j1 = 1
        j2 = Nt+1
        if (grid%th(1) < pole_tol) j1 = 2
        if (grid%th(Nt+1) > PI - pole_tol) j2 = Nt

        ! H%y -- theta component
        do j = j1,j2
            call legendre_norm(lmax, cos(grid%th(j)), P_lm)
            do i = 1,Np
                call vsharm(lmax, cos(grid%th(j)), grid%ph(i)+grid%dp(i)/2, P_lm, Y, Yt)
                do k = 1,Nrr
                    csum = dcmplx(0.0d0,0.0d0)
                    do l = 1,lmax
                        do m = -l,l
                            csum = csum + coeff2(l,m)*Tnrp(k,l)/Rr(k)*Yval(Yt,l,m,lmax)
                        end do
                    end do
                    H%y(i,j,k) = csum
                end do
            end do
        end do

        ! H%x -- phi component
        do j = 1,Nt
            call legendre_norm(lmax, cos(grid%th(j)+grid%dt(j)/2), P_lm)
            do i = 1,Np+1
                call vsharm(lmax, cos(grid%th(j)+grid%dt(j)/2), grid%ph(i), P_lm, Y, Yt, Yp)
                do k = 1,Nrr
                    csum = dcmplx(0.0d0,0.0d0)
                    do l = 1,lmax
                        do m = -l,l
                            csum = csum + coeff2(l,m)*Tnrp(k,l)/Rr(k)*Yval(Yp,l,m,lmax)
                        end do
                    end do
                    H%x(i,j,k) = csum
                end do
            end do
        end do

        ! H%z -- radial component
        do j = 1,Nt
            call legendre_norm(lmax, cos(grid%th(j)+grid%dt(j)/2), P_lm)
            do i = 1,Np
                call vsharm(lmax, cos(grid%th(j)+grid%dt(j)/2), grid%ph(i)+grid%dp(i)/2, P_lm, Y)
                do k = 1,Nrs
                    csum = dcmplx(0.0d0,0.0d0)
                    do l = 1,lmax
                        do m = -l,l
                            csum = csum + coeff2(l,m)*Tns(k,l)*dble(l*(l+1))/Rs(k)**2*Yval(Y,l,m,lmax)
                        end do
                    end do
                    H%z(i,j,k) = csum
                end do
            end do
        end do

        if (present(E)) then
        ! E%x -- phi component (node theta -- j1/j2 computed above, same
        ! node-theta range as H%y). The below logic allows for global vs regional grids;
        ! it's good practice to define regional grids away from the poles - for
        ! polar regions, will need to rotate to fake pole (supported)
        do j = j1,j2
            call legendre_norm(lmax, cos(grid%th(j)), P_lm)
            do i = 1,Np
                call vsharm(lmax, cos(grid%th(j)), grid%ph(i)+grid%dp(i)/2, P_lm, Y, Yt)
                do k = 1,Nrs
                    csum = dcmplx(0.0d0,0.0d0)
                    do l = 1,lmax
                        do m = -l,l
                            csum = csum - coeff2(l,m)*iomu0*Tns(k,l)/Rs(k)*Yval(Yt,l,m,lmax)
                        end do
                    end do
                    E%x(i,j,k) = csum
                end do
            end do
        end do

        ! E%y -- theta component
        do j = 1,Nt
            call legendre_norm(lmax, cos(grid%th(j)+grid%dt(j)/2), P_lm)
            do i = 1,Np+1
                call vsharm(lmax, cos(grid%th(j)+grid%dt(j)/2), grid%ph(i), P_lm, Y, Yt, Yp)
                do k = 1,Nrs
                    csum = dcmplx(0.0d0,0.0d0)
                    do l = 1,lmax
                        do m = -l,l
                            csum = csum + coeff2(l,m)*iomu0*Tns(k,l)/Rs(k)*Yval(Yp,l,m,lmax)
                        end do
                    end do
                    E%y(i,j,k) = csum
                end do
            end do
        end do
        end if ! present(E)

        else
            write(0,*) 'Error in sourceField1d_s2: unknown H%gridType: ', trim(H%gridType)
        end if

        deallocate(Tnr, Tnsp, STAT=istat)
        if (H%gridType == FACE) deallocate(Tnrp, Tns, STAT=istat)
        deallocate(Rr, Rs, coeff2, STAT=istat)

    end subroutine sourceField1d_s2

    !=======================================================================
    ! surfaceField1d_s2 (S2): evaluate the full EM field at the Earth's
    ! SURFACE (r = earth%r0 = R_earth) at the LATERAL CELL CENTRES of the
    ! grid (theta = grid%th(j)+grid%dt(j)/2, phi = grid%ph(i)+grid%dp(i)/2),
    ! for validation and plotting. S2 counterpart of field1d.f90's
    ! surfaceField1d -- same Sun & Egbert (2012) eq (5)-(6) assembly as
    ! sourceField1d_s2's EDGE branch, specialized to the single radius r=R0.
    ! Returns all six components co-located at one common set of surface
    ! points, native e^{-i*omega*t}, in Hsurf/Esurf(Np,Nt,3), component
    ! index 1=radial (r), 2=colatitude (theta), 3=longitude (phi); E_r=0.
    !=======================================================================
    subroutine surfaceField1d_s2(earth,lmax,coeff,period,grid,Hsurf,Esurf)

        type(conf1d_t), intent(in)                  :: earth
        integer, intent(in)                         :: lmax
        complex(8), dimension(:), intent(in)        :: coeff
        real(8), intent(in)                         :: period
        type(grid_t), intent(in)                    :: grid
        complex(8), dimension(:,:,:), intent(out)   :: Hsurf, Esurf ! (Np,Nt,3)
        ! local
        real(8), dimension(1)                       :: Rr, Rs
        complex(8), dimension(1,lmax)               :: Tnr, Tnsp
        real(8), dimension(lmax+1,lmax+1)           :: P_lm
        complex(8), dimension(lmax,lmax+1)          :: Y, Yt, Yp
        complex(8), dimension(:,:), allocatable     :: coeff2
        real(8)     :: R0, th_mid, ph_mid, omega
        complex(8)  :: iomu0, csr, cst, csp
        integer     :: Np, Nt, i, j, l, m, istat

        Np = grid%nx
        Nt = grid%ny
        R0 = earth%r0
        omega = 2.0d0*PI/period
        iomu0 = dcmplx(0.0d0,1.0d0) * omega * MU_0

        ! unpack flat, degree-blocked coeff(:) into coeff2(l,-l:l) (as in
        ! sourceField1d_s2): degree-l block starts at l^2+1, size 2l+1,
        ! ordered m=0,1,-1,2,-2,...
        allocate(coeff2(lmax,-lmax:lmax), STAT=istat)
        coeff2(:,:) = dcmplx(0.0d0,0.0d0)
        do l = 1,lmax
            coeff2(l,0) = coeff(l*l+1)
            do m = 1,l
                coeff2(l,m)  = coeff(l*l+2*m)
                coeff2(l,-m) = coeff(l*l+2*m+1)
            end do
        end do

        ! source potentials at the single surface radius r = R_earth
        Rr(1) = R0
        Rs(1) = R0
        call sourcePotential_s2(earth,lmax,period,Rr,Rs,Tnr,Tnsp)

        if (any(Tnr /= Tnr) .or. any(Tnsp /= Tnsp)) then
            deallocate(coeff2, STAT=istat)
            call errStop('surfaceField1d_s2: NaN in source potentials at the surface')
        end if

        Hsurf = dcmplx(0.0d0,0.0d0)
        Esurf = dcmplx(0.0d0,0.0d0)

        do j = 1,Nt
            th_mid = grid%th(j) + grid%dt(j)/2
            call legendre_norm(lmax,cos(th_mid),P_lm)
            do i = 1,Np
                ph_mid = grid%ph(i) + grid%dp(i)/2
                call vsharm(lmax,cos(th_mid),ph_mid,P_lm,Y,Yt,Yp)

                csr = dcmplx(0.0d0,0.0d0)
                cst = dcmplx(0.0d0,0.0d0)
                csp = dcmplx(0.0d0,0.0d0)
                do l = 1,lmax
                    do m = -l,l
                        ! Hr = (l(l+1)/r^2) T(r) Y ; at r=R0 -> l(l+1)/R0^2
                        csr = csr + coeff2(l,m)*Tnr(1,l)*dble(l*(l+1))/R0**2*Yval(Y,l,m,lmax)
                        ! Htheta = (1/r) T'(r) Yt ; at r=R0 -> Tnsp/R0
                        cst = cst + coeff2(l,m)*Tnsp(1,l)/R0*Yval(Yt,l,m,lmax)
                        ! Hphi = (1/r) T'(r) Yp
                        csp = csp + coeff2(l,m)*Tnsp(1,l)/R0*Yval(Yp,l,m,lmax)
                    end do
                end do
                Hsurf(i,j,1) = csr
                Hsurf(i,j,2) = cst
                Hsurf(i,j,3) = csp

                cst = dcmplx(0.0d0,0.0d0)
                csp = dcmplx(0.0d0,0.0d0)
                do l = 1,lmax
                    do m = -l,l
                        ! Etheta = (i*omega*mu0/r) T(r) Yp ; at r=R0 -> iomu0*Tnr/R0
                        cst = cst + coeff2(l,m)*iomu0*Tnr(1,l)/R0*Yval(Yp,l,m,lmax)
                        ! Ephi = -(i*omega*mu0/r) T(r) Yt
                        csp = csp - coeff2(l,m)*iomu0*Tnr(1,l)/R0*Yval(Yt,l,m,lmax)
                    end do
                end do
                Esurf(i,j,2) = cst
                Esurf(i,j,3) = csp
                ! Esurf(i,j,1) = 0 (toroidal field, E_r=0)
            end do ! phi
        end do ! theta

        deallocate(coeff2, STAT=istat)

    end subroutine surfaceField1d_s2

end module field1d_s2
