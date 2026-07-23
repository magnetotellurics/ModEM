"""
Check whether Faraday's law, curl(E) = i*omega*mu0*H, holds on global1d's
ACTUAL .hfield/.efield output files, via numerical (finite-difference)
differentiation on the grid -- companion to test_pythonsolver_faraday.py,
which did the same self-consistency check on spherical_em_induction.py's
analytic fields instead. Together these two tests localize whether the
l=1,m=-1 "(-1i) on Hz,Ex,Ey but not Hx,Hy" discrepancy (see CLAUDE.md,
2026-07-22) is a field1d.f90 bug, a pythonSolver bug, or neither.

Only the RADIAL component of Faraday's law is checked:

    (curl E)_r = (1/(r sin(theta))) * [ d(sin(theta)*Ephi)/dtheta - dEtheta/dphi ]
               = i*omega*mu0*Hr

IMPORTANT -- grid staggering (this is a REVISION of a first, WRONG attempt
that naively compared Etheta/Ephi/Hr at the same raw (i,j,k) file index,
which gave meaningless noise): with primary_grid='E' (hardcoded in
FWD1D.f90), E%gridType=EDGE and H%gridType=FACE, so Etheta, Ephi, Hr live on
DIFFERENT, staggered sub-grids, not the same physical (theta,phi) points.
Per field1d.f90's FACE-branch loop headers (see CLAUDE.md):
    H%z (r,     FACE z): mid_theta  x mid_phi  x Rs
    E%x (phi,   EDGE x): node_theta x mid_phi  x Rs
    E%y (theta, EDGE y): mid_theta  x node_phi x Rs
and confirmed against write_cvector's embedding rule (EDGE: x(1:Nx,:,:)=E%x,
y(:,1:Ny,:)=E%y; FACE: z(1:Nx,1:Ny,:)=E%z), the raw file arrays (0-indexed
here) correspond to:
    Hr[i,j,k]     at (phi_mid(i),  theta_mid(j))     i=0..nx-1, j=0..ny-1
    Ephi[i,j,k]   at (phi_mid(i),  theta_node(j))    i=0..nx-1, j=0..ny
    Etheta[i,j,k] at (phi_node(i), theta_mid(j))     i=0..nx,   j=0..ny-1
where theta_node(j)=j*dtheta, theta_mid(j)=theta_node(j)+dtheta/2 (same
pattern for phi). This means, for a chosen Hr point (i,j): Ephi[i,j,k] and
Ephi[i,j+1,k] sit EXACTLY at theta_mid(j)-+dtheta/2 (bracketing Hr's theta
with no interpolation needed -- the whole point of Yee-grid staggering), and
similarly Etheta[i,j,k] / Etheta[i+1,j,k] bracket Hr's phi. So the curl at
Hr's own location is a simple, DIRECT forward-type difference using exactly
these neighbors -- no interpolation, no assumed sign flips.
"""
import numpy as np

EFIELD_FILE = r"C:\Users\Anna Kelbert\Developer\ModEM-global1d\MTsource\MT.1000sec.test_m-ve.sunegbert2012.E-grid.T01.efield"
HFIELD_FILE = r"C:\Users\Anna Kelbert\Developer\ModEM-global1d\MTsource\MT.1000sec.test_m-ve.sunegbert2012.E-grid.T01.hfield"

PERIOD = 1000.0
omega = 2 * np.pi / PERIOD
mu0 = 1.256637e-6   # matches the constant hardcoded in field1d.f90 / TSModel.m

R_INDEX = 13         # 1-indexed; k=13 = first layer at/below the surface (see CLAUDE.md)

# test point: comfortably interior (away from poles, away from the phi=0/360
# wraparound seam), 0-indexed into Hr's own (mid_phi, mid_theta) grid
I_H = 89     # phi_mid index  (0..nx-1)
J_H = 44     # theta_mid index (0..ny-1)


def read_field(fname):
    with open(fname) as f:
        f.readline()  # header comment line
        f.readline()  # "  1"
        nx, ny, nz, gridType = f.readline().split()
        nx, ny, nz = int(nx), int(ny), int(nz)
        gridType = gridType.strip()

        n = (nx + 1) * (ny + 1) * (nz + 1)
        data = np.loadtxt(f, max_rows=n)

    Cphi = data[:, 3] + 1j * data[:, 4]     # phi component
    Ctheta = data[:, 5] + 1j * data[:, 6]   # theta component
    Cr = data[:, 7] + 1j * data[:, 8]       # r component

    shape = (nx + 1, ny + 1, nz + 1)
    Cphi = Cphi.reshape(shape, order='F')
    Ctheta = Ctheta.reshape(shape, order='F')
    Cr = Cr.reshape(shape, order='F')

    return dict(phi=Cphi, theta=Ctheta, r=Cr, nx=nx, ny=ny, nz=nz, gridType=gridType)


Efield = read_field(EFIELD_FILE)
Hfield = read_field(HFIELD_FILE)
print(f"E file: nx={Efield['nx']} ny={Efield['ny']} nz={Efield['nz']} gridType={Efield['gridType']}")
print(f"H file: nx={Hfield['nx']} ny={Hfield['ny']} nz={Hfield['nz']} gridType={Hfield['gridType']}")
assert Efield['gridType'] == 'EDGE' and Hfield['gridType'] == 'FACE', \
    "staggering assumptions below are specific to E=EDGE, H=FACE (primary_grid='E')"

nx, ny = Efield['nx'], Efield['ny']
dtheta = np.pi / ny
dphi = 2 * np.pi / nx
k = R_INDEX - 1

theta_node = lambda j: j * dtheta
theta_mid = lambda j: theta_node(j) + dtheta / 2
phi_mid_val = (I_H + 0.5) * dphi
theta_mid_val = theta_mid(J_H)
r0 = 6.371e6   # approx radius at this near-surface layer, only for the 1/(r sin theta) prefactor

print(f"\nHr test point: i={I_H} j={J_H}  ->  phi_mid={np.degrees(phi_mid_val):.2f} deg, "
      f"theta_mid={np.degrees(theta_mid_val):.2f} deg, k={R_INDEX}")

Etheta = Efield['theta']   # shape (nx+1, ny, nz+1)  -- node_phi x mid_theta
Ephi = Efield['phi']       # shape (nx, ny+1, nz+1)  -- mid_phi  x node_theta
Hr = Hfield['r']           # shape (nx, ny, nz+1)    -- mid_phi  x mid_theta

# d(sin(theta)*Ephi)/dtheta at (phi_mid(I_H), theta_mid(J_H)): Ephi[I_H,J_H,k]
# and Ephi[I_H,J_H+1,k] sit EXACTLY at theta_node(J_H) and theta_node(J_H+1),
# which bracket theta_mid(J_H) -- no interpolation needed.
sinEphi_lo = np.sin(theta_node(J_H)) * Ephi[I_H, J_H, k]
sinEphi_hi = np.sin(theta_node(J_H + 1)) * Ephi[I_H, J_H + 1, k]
dEphi_dtheta = (sinEphi_hi - sinEphi_lo) / dtheta

# dEtheta/dphi at the same point: Etheta[I_H,J_H,k] and Etheta[I_H+1,J_H,k]
# sit EXACTLY at phi_node(I_H), phi_node(I_H+1), bracketing phi_mid(I_H).
dEtheta_dphi = (Etheta[I_H + 1, J_H, k] - Etheta[I_H, J_H, k]) / dphi

curlE_r = (dEphi_dtheta - dEtheta_dphi) / (r0 * np.sin(theta_mid_val))
rhs = 1j * omega * mu0 * Hr[I_H, J_H, k]

print(f"\n(curl E)_r (numerical, staggering-aware) = {curlE_r}")
print(f"i*omega*mu0*Hr (raw Hr from file)         = {rhs}")
print(f"ratio (curl E)_r / (i*omega*mu0*Hr)        = {curlE_r / rhs}")
print(f"  (expect ~+1.0 if raw file components directly satisfy Faraday's law;")
print(f"   ~-1 would indicate a sign-convention mismatch (e.g. the 'down' Hr/Er")
print(f"   relabeling found earlier for MT-convention output); ~+-1j would")
print(f"   indicate a genuine phase bug)")

relerr = abs(curlE_r - rhs) / abs(rhs)
print(f"\nrelative error (vs ratio=+1) = {relerr:.3e}")
