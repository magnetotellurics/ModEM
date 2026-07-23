from scipy.special import lpmv, factorial
import numpy as np
from spherical_em_induction import Y_lm

x = 0.3
for m in [0,1,2]:
    for l in [1,2,3]:
        if m>l: continue
        print(f'l={l} m={m}  lpmv={lpmv(m,l,x):.6f}')
print()
# Standard (Condon-Shortley) closed forms for comparison, m=1:
# P_1^1(x) = -sqrt(1-x^2)
# P_2^1(x) = -3x*sqrt(1-x^2)
print('Expected WITH Condon-Shortley phase (physics-standard, (-1)^m built into P_l^m for m>0):')
print('P_1^1(x) =', -np.sqrt(1-x**2))
print('P_2^1(x) =', -3*x*np.sqrt(1-x**2))
print()
print('Expected WITHOUT Condon-Shortley phase (sign flipped):')
print('P_1^1(x) =', np.sqrt(1-x**2))
print('P_2^1(x) =', 3*x*np.sqrt(1-x**2))
print()

theta, phi = 1.1, 0.7   # arbitrary generic angles

for l in [1, 2, 3]:
    for m in range(1, l+1):
        Yp = Y_lm(l, m, theta, phi)
        Yn = Y_lm(l, -m, theta, phi)
        expected = (-1)**m * np.conj(Yp)
        match = np.isclose(Yn, expected)
        print(f'l={l} m={m}:  Y(-m) = {Yn:.6f}   (-1)^m*conj(Y(m)) = {expected:.6f}   match={match}')
