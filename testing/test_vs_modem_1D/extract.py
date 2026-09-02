import sys

def read_field(fname):
    with open(fname) as f:
        header = f.readline()
        f.readline()  # "  1"
        dims = f.readline().split()
        nx, ny, nz = int(dims[0]), int(dims[1]), int(dims[2])
        gridType = dims[3]
        data = {}
        for line in f:
            parts = line.split()
            i, j, k = int(parts[0]), int(parts[1]), int(parts[2])
            rex, imx, rey, imy, rez, imz = map(float, parts[3:9])
            data[(i, j, k)] = (rex, imx, rey, imy, rez, imz)
    return header.strip(), nx, ny, nz, gridType, data

for fname in sys.argv[1:]:
    header, nx, ny, nz, gridType, data = read_field(fname)
    # near-surface (small k = near top of air->just below surface varies;
    # pick a k a few layers below where field is well-established), interior i,j
    i0, j0 = nx // 2, ny // 2
    # find a k near the surface -- scan for a non-tiny value at (i0,j0)
    print(f"--- {fname} ---  nx={nx} ny={ny} nz={nz} gridType={gridType}")
    for k in (10, 13, 15, 20, 30):
        if (i0, j0, k) in data:
            rex, imx, rey, imy, rez, imz = data[(i0, j0, k)]
            print(f"  i={i0} j={j0} k={k}:  Ex(theta,'north')=Re,Im(y)= {rey:+.6e} {imy:+.6e}   "
                  f"Ey(phi,'east')=Re,Im(x)= {rex:+.6e} {imx:+.6e}")
