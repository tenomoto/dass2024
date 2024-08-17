import numpy as np

# l63
ns = 3
xt0 = np.array([1, 3, 5])
xb0 = np.array([1.1, 3.3, 5.5])
#xb0 = np.array([1, -1, 1])
p = 10
r = 48
b = 8 / 3
dt = 0.01
nt = 200
model_q = 0.1
obs_int = 60
obs_r = np.array([0.1, 0.3, 0.5])**2
nobs = nt // obs_int

# var
bgd_b = 1
ni = 10
a = 5e-4

# ens
ne = 10
fil = "enkf"
