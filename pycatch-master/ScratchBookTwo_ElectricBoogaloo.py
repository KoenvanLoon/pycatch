import skgstat as skg
import gstools as gs
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('ggplot')

dim = 3
# rotation around z, y ,x
angles = [np.deg2rad(90), np.deg2rad(45), np.deg2rad(22.5)]
model = gs.Gaussian(dim=3, len_scale=[16, 8, 4], angles=angles)
x = y = z = range(50)
pos = (x, y , z)
srf = gs.SRF(model, seed=1001)
field = srf.structured(pos)

# All 3 axes of the rotated coordinate-system
main_axes = gs.rotated_main_axes(dim, angles)
axis1, axis2, axis3 = main_axes

# in generalfunctions - R
# a = variogram(v ~ x, ~ x + y + z, dataFrame, beta=0,tol.ver=0.1, alpha=90,tol.hor=0.1, boundaries=boundariesVector)
## v ~ x or v ~ 1 == object, ~1 in case of absence of regressors; formula defining the response vector and (possible) regressors
## ~ x + y + z == coords/position # CHECK
## dataFrame == data/field # CHECK
## alpha=90 == direction in plane (x,y), in positive degrees clockwise from positive y (North): alpha=0 for direction North (increasing y), alpha=90 for direction East (increasing x); optional a vector of directions in (x,y)
## beta=0 == direction in z, in positive degrees up from the (x,y) plane; optional a vector of directions
## tol.hor=0.1 == horizontal tolerance angle in degrees # CHECK (angles_tol GS; hor and ver?)
## boundaries=boundariesVector == numerical vector with distance interval upper boundaries

boundariesVector = (30.5, 40.5)

# GS method
bin_center, dir_vario, counts = gs.vario_estimate(
    pos,
    field,
    direction=main_axes,
    bandwidth=10,
    sampling_size=2000,
    sampling_seed=1001,
    mesh_type="structured",
    return_counts=True
)
# returns the bin centers and the estimated variogram values at the bin centers

# bin_edges (numpy.ndarray, optional) - the bins on which the variogram will be calculated
# --- size of the bins; not needed
# sampling_size (int or none, optional) - for large input data, this method can take a long time to compute the variogram, therefore this argument specifies the number of data points to sample randomly
# --- useful to speed up model; tradeoff with return
# direction(list of numpy.ndarray, optional) - directions to evaluate a directional variogram. Anglular tolerance is given by angles_tol. bandwidth to cut off how wide the search for point pairs should be is given by bandwidth. You can provide multiple directions at once to get one variogram for each direction. For a single direction you can also use the angles parameter, to provide the direction by its spherical coordianates.
# --- 1 direction only
# angles (numpy.ndarray, optional) - the angles of the main axis to calculate the variogram for in radians angle definitions from ISO standard 80000-2:2009 for 1d this parameter will have no effect at all for 2d supply one angle which is azimuth \varphi (ccw from +x in xy plane) for 3d supply two angles which are azimuth \varphi (ccw from +x in xy plane) and inclination \theta (cw from +z). Can be used instead of direction.
# --- replaces R alpha/beta?
# angles_tol (class: float, optional) - the tolerance around the variogram angle to count a point as being within this direction from another point (the angular tolerance around the directional vector given by angles) Default: np.pi/8 = 22.5°
# --- replaces tol.hor (and tol.ver?)
# bandwith (class: float, optional) - bandwidth to cut off the angular tolerance for directional variograms. If None is given, only the angles_tol parameter will control the point selection.

print("Original:")
print(model)
model.fit_variogram(bin_center, dir_vario)
print("Fitted:")
print(model)

# Plots
fig = plt.figure(figsize=[10, 5])
ax1 = fig.add_subplot(121, projection=Axes3D.name)
ax2 = fig.add_subplot(122)

ax1.plot([0, axis1[0]], [0, axis1[1]], [0, axis1[2]], label="0.")
ax1.plot([0, axis2[0]], [0, axis2[1]], [0, axis2[2]], label="1.")
ax1.plot([0, axis3[0]], [0, axis3[1]], [0, axis3[2]], label="2.")
ax1.set_xlim(-1, 1)
ax1.set_ylim(-1, 1)
ax1.set_zlim(-1, 1)
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.set_zlabel("Z")
ax1.set_title("Tait-Bryan main axis")
ax1.legend(loc="lower left")

x_max = max(bin_center)
ax2.scatter(bin_center, dir_vario[0], label="0. axis")
ax2.scatter(bin_center, dir_vario[1], label="1. axis")
ax2.scatter(bin_center, dir_vario[2], label="2. axis")
model.plot("vario_axis", axis=0, ax=ax2, x_max=x_max, label="fit on axis 0")
model.plot("vario_axis", axis=1, ax=ax2, x_max=x_max, label="fit on axis 1")
model.plot("vario_axis", axis=2, ax=ax2, x_max=x_max, label="fit on axis 2")
ax2.set_title("Fitting an anisotropic model")
ax2.legend()

plt.show()

srf.plot()
plt.show()
