import cantera as ct
import numpy as np

np.set_printoptions(precision=2, linewidth=144)

surf = ct.Interface('dummysurf2.yml', 'surf')
surf2 = ct.Interface('dummysurf2.yml', 'surf2')
gas = surf.adjacent['gas']
gas.TPX =  1000, 2e5, 'A:0.6, B:0.3, C:0.2, D:0.1'
surf.coverages = 'A(S):0.1, B(S):0.2, C(S):0.3, D(S):0.2, (S):0.2'
surf2.coverages = 'A(S):0.1, D(S):0.2, (S):0.2'

surf.derivative_settings = {'skip-coverage-dependence': True, 'skip-electrochemistry': True}
surf2.derivative_settings = {'skip-coverage-dependence': True, 'skip-electrochemistry': True}
r = ct.IdealGasMoleReactor(gas)
r.volume = 3
rsurf2 = ct.ReactorSurface(surf2, r, A=9e-4)
rsurf1 = ct.ReactorSurface(surf, r, A=5e-4)
net = ct.ReactorNet([r])
net.step()
print('Analytical Jacobian:')
print(r.jacobian)
print()
print('Finite difference Jacobian')
print(r.finite_difference_jacobian)
