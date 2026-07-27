"""
Reader for ModEM's binary .esoln E-field solution file (written by
ioAscii.f90's write_solnVectorMTX -- despite the module's name, this is
Fortran UNFORMATTED BINARY, not ASCII). Format (from direct inspection of
C:\\Users\\Anna Kelbert\\Developer\\ModEM\\f90\\3D_MT\\ioMod\\ioAscii.f90,
FileWriteInit lines 188-215 and EfileWrite lines 376-392):

Header (4 records):
  1: version(char20), nPer(int4), nMode(int4), nx(int4), ny(int4), nz(int4),
     nzAir(int4), ox(f8), oy(f8), oz(f8), rotdeg(f8)
  2: dx(1:nx)  f8 array
  3: dy(1:ny)  f8 array
  4: dz(1:nz)  f8 array

Then for j=1..nPer (outer), k=1..nMode (inner), 4 records per block:
  a: Omega(f8), iFreq(int4), iMode(int4), ModeName(char20)
  b: Ex  complex16, shape (nx,   ny+1, nz+1)  [Fortran order]
  c: Ey  complex16, shape (nx+1, ny,   nz+1)
  d: Ez  complex16, shape (nx+1, ny+1, nz)

iMode=1 -> ModeName='Ey' (source polarization), iMode=2 -> ModeName='Ex' --
per boundary_ws[S].f90's actual physics (iMode=1 sets only E0%y, iMode=2 sets
only E0%x; SolnSpace.f90's Pol_name labels were fixed 2026-07 to match, see
ModEM's own history). Matches this project's own Mode1(l=1,m=0,zonal,
Zyx=Ey/Hx)/Mode2(l=1,m=+-1,Zxy=Ex/Hy) naming (confirmed via
small_predicted.dat's ZXY/ZYX column correspondence).
"""
import sys
import numpy as np
from scipy.io import FortranFile


def read_esoln(fname):
    f = FortranFile(fname, 'r')

    rec1 = f.read_record('S20', 'i4', 'i4', 'i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8')
    version = rec1[0][0].decode(errors='replace').strip()
    nPer, nMode, nx, ny, nz, nzAir = (int(rec1[i][0]) for i in range(1, 7))
    ox, oy, oz, rotdeg = (float(rec1[i][0]) for i in range(7, 11))

    dx = f.read_record('f8')
    dy = f.read_record('f8')
    dz = f.read_record('f8')

    blocks = []
    for j in range(nPer):
        for k in range(nMode):
            hdr = f.read_record('f8', 'i4', 'i4', 'S20')
            omega = float(hdr[0][0])
            iFreq = int(hdr[1][0])
            iMode = int(hdr[2][0])
            modeName = hdr[3][0].decode(errors='replace').strip()

            Ex = f.read_record(np.complex128).reshape((nx, ny + 1, nz + 1), order='F')
            Ey = f.read_record(np.complex128).reshape((nx + 1, ny, nz + 1), order='F')
            Ez = f.read_record(np.complex128).reshape((nx + 1, ny + 1, nz), order='F')

            blocks.append(dict(omega=omega, period=2 * np.pi / omega, iFreq=iFreq,
                                iMode=iMode, modeName=modeName, Ex=Ex, Ey=Ey, Ez=Ez))

    f.close()
    return dict(version=version, nPer=nPer, nMode=nMode, nx=nx, ny=ny, nz=nz, nzAir=nzAir,
                ox=ox, oy=oy, oz=oz, rotdeg=rotdeg, dx=dx, dy=dy, dz=dz, blocks=blocks)


if __name__ == "__main__":
    fname = sys.argv[1] if len(sys.argv) > 1 else "out.esoln"
    d = read_esoln(fname)
    print(f"version={d['version']!r} nPer={d['nPer']} nMode={d['nMode']} "
          f"nx={d['nx']} ny={d['ny']} nz={d['nz']} nzAir={d['nzAir']}")
    print(f"ox={d['ox']} oy={d['oy']} oz={d['oz']} rotdeg={d['rotdeg']}")
    print(f"dx[:5]={d['dx'][:5]} ... dx[-5:]={d['dx'][-5:]}")
    print(f"dy[:5]={d['dy'][:5]} ... dy[-5:]={d['dy'][-5:]}")
    print(f"dz[:8]={d['dz'][:8]}")
    for b in d['blocks']:
        print(f"  block: iFreq={b['iFreq']} iMode={b['iMode']} modeName={b['modeName']!r} "
              f"period={b['period']:.3f}s  Ex.shape={b['Ex'].shape}")
