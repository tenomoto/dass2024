import numpy as np


def analysis(zf, yo, r, loc_inf=1.0):
  k, ne = zf.shape
  zf_mean = zf.mean(axis=1)
  zf_anom = zf - zf_mean[:, None]
  s = loc_inf * (zf_anom @ zf_anom[-1,]) / (ne - 1)
  return np.outer(s, (yo - zf[-1, :])) / (r + s[-1])

