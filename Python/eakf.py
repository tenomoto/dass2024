import numpy as np


def analysis(zf, yo, r, loc_inf=1.0):
  k, ne = zf.shape
  zf_mean = zf.mean(axis=1)
  zf_anom = zf - zf_mean[:, None]
  s = loc_inf * (zf_anom @ zf_anom[-1,]) / (ne - 1)
  alpha = np.sqrt(r / (r + s[-1]))
  dz_mean = s / (r + s[-1]) * (yo - zf_mean[-1])
  return dz_mean[:, None] + (alpha - 1) / s[-1] * np.outer(s, zf_anom[-1,])

