import cantera as ct
import numpy as np
import math

np.set_printoptions(precision=1, linewidth=np.inf)

model = "h2o2.yaml"
# initial conditions
T0 = 1500
P0 = ct.one_atm
equiv_ratio = 1
surf_area = 0.527
fuel = "H2"
air = "O2:1.0, N2:3.773"
# reactor 1
gas1 = ct.Solution(model, transport_model=None)
gas1.TP = T0, P0
gas1.set_equivalence_ratio(equiv_ratio, fuel, air)
r1 = ct.IdealGasMoleReactor(gas1)
# # surf 1
# surf1 = ct.Interface(model, "h_surf", [gas1])
# surf1.TP = T0, P0
# surf1.coverages = {"PT(S)":1}
# rsurf1 = ct.ReactorSurface(surf1, r1, A=surf_area)
# reactor network setup
net1 = ct.ReactorNet([r1,])
net1.initialize()

om = lambda number: math.floor(math.log(number, 10))
jac = r1.jacobian
fd_jac = r1.finite_difference_jacobian
for i in range(r1.n_vars):
    for j in range(r1.n_vars):
        if jac[i, j] != 0:
            jac[i, j] = np.sign(jac[i, j]) * om(abs(jac[i, j]))
        if fd_jac[i, j] != 0:
            fd_jac[i, j] = np.sign(fd_jac[i, j]) * om(abs(fd_jac[i, j]))
print(jac)
# print(fd_jac)
ctr = 0
for i in range(r1.n_vars):
    for j in range(r1.n_vars):
        if jac[i, j] == fd_jac[i, j]:# and jac[i, j] != 0:
            ctr += 1
print(f"{ctr} / {r1.n_vars ** 2}")
        # print(jac[i, j], fd_jac[i, j])
# print(np.sign(r1.jacobian) == np.sign(r1.finite_difference_jacobian))
