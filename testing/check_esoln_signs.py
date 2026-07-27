from read_esoln import read_esoln

d = read_esoln("out.esoln")
nx, ny, nz, nzAir = d['nx'], d['ny'], d['nz'], d['nzAir']
print(f"nx={nx} ny={ny} nz={nz} nzAir={nzAir}  dx={d['dx'][:3]}... dy={d['dy'][:3]}...")

# Station at X=4.000 Y=4.000 Z=0.000 (small_test.dat) -- with dx=dy=1 (8 cells,
# origin 0), node j/i=5 (1-indexed) sits at exactly y or x = 4*dy = 4, matching
# the station exactly. k=13 (1-indexed) is the surface node (12 air layers
# above, matching this whole project's R_INDEX=13 convention). Model is
# laterally uniform (1D) so the exact horizontal cell-center choice for the
# "other" axis shouldn't matter.
i_node, j_node, k_node = 5, 5, 13   # 1-indexed node position
i_mid, j_mid = 4, 4                  # 1-indexed cell-center position (any reasonable choice)

for b in d['blocks']:
    if b['iMode'] == 1:   # Mode1 (Ey source polarization, per boundary_ws[S].f90) -- want Ey
        # Ey shape (nx+1, ny, nz+1): x=node(i_node), y=cell-center(j_mid), z=node(k_node)
        val = b['Ey'][i_node - 1, j_mid - 1, k_node - 1]
        print(f"T={b['period']:7.1f}s  Mode1 (Ey-pol)  Ey(i={i_node},j={j_mid},k={k_node}) = "
              f"{val.real:+.6e} {val.imag:+.6e}j")
    else:                  # Mode2 (Ex source polarization, per boundary_ws[S].f90) -- want Ex
        # Ex shape (nx, ny+1, nz+1): x=cell-center(i_mid), y=node(j_node), z=node(k_node)
        val = b['Ex'][i_mid - 1, j_node - 1, k_node - 1]
        print(f"T={b['period']:7.1f}s  Mode2 (Ex-pol)  Ex(i={i_mid},j={j_node},k={k_node}) = "
              f"{val.real:+.6e} {val.imag:+.6e}j")
