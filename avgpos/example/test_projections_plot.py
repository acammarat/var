#!/usr/bin/env python3
"""
Matplotlib script to visualize plane projection data as a smooth interpolated heatmap.
Generated automatically by avgpos tool.

Usage: python3 test_projections_plot.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import Rbf
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Read data from file
data = np.loadtxt('test_projections.dat', dtype=str)

# Extract coordinates and g values
e_orig = data[:, 0].astype(float)
f_orig = data[:, 1].astype(float)
g_orig = data[:, 2].astype(float)

# Extract labels if present (4th column)
has_labels = data.shape[1] > 3
if has_labels:
    labels_orig = data[:, 3]

# Replication parameters
ne_rep, nf_rep = 2.0, 2.0

# Replicate data along e and f axes
e_list, f_list, g_list, labels_list = [], [], [], []

# Use lattice-based shifts if provided, otherwise use data range (legacy)
lattice_shifts = (np.float64(1.6616715304787615), np.float64(2.878099516279951))
if lattice_shifts is not None:
    e_shift_unit = lattice_shifts[0]
    f_shift_unit = lattice_shifts[1]
else:
    # Legacy behavior: use data range
    e_shift_unit = e_orig.max() - e_orig.min() if len(e_orig) > 1 else 1.0
    f_shift_unit = f_orig.max() - f_orig.min() if len(f_orig) > 1 else 1.0

# Determine how many full and partial replications to make
ne_full = int(np.ceil(ne_rep))
nf_full = int(np.ceil(nf_rep))

for ie in range(ne_full):
    for jf in range(nf_full):
        # Calculate if this replica is fully or partially included
        e_factor = min(1.0, ne_rep - ie) if ie < ne_full - 1 else (ne_rep - ie)
        f_factor = min(1.0, nf_rep - jf) if jf < nf_full - 1 else (nf_rep - jf)
        
        # Include this replica if it has non-zero contribution
        if e_factor > 0 and f_factor > 0:
            e_shift = ie * e_shift_unit
            f_shift = jf * f_shift_unit
            e_list.append(e_orig + e_shift)
            f_list.append(f_orig + f_shift)
            g_list.append(g_orig)
            if has_labels:
                labels_list.append(labels_orig)

e = np.concatenate(e_list)
f = np.concatenate(f_list)
g = np.concatenate(g_list)
if has_labels:
    labels = np.concatenate(labels_list)

# Create a regular grid for interpolation
# Determine the range of e and f
e_min, e_max = e.min(), e.max()
f_min, f_max = f.min(), f.max()

e_grid = np.linspace(e_min, e_max, 200)
f_grid = np.linspace(f_min, f_max, 200)
e_mesh, f_mesh = np.meshgrid(e_grid, f_grid)

# Use Radial Basis Function interpolation with very small smoothing
# This ensures interpolation passes extremely close to data points while handling duplicates
# smooth value is set to a tiny value to get nearly exact values at data points
try:
    rbf = Rbf(e, f, g, function='thin_plate', smooth=1e-10)
    g_interp = rbf(e_mesh, f_mesh)
except:
    # Fall back to multiquadric if thin_plate fails
    try:
        rbf = Rbf(e, f, g, function='multiquadric', smooth=1e-10)
        g_interp = rbf(e_mesh, f_mesh)
    except:
        # Last resort: use linear with small smoothing
        rbf = Rbf(e, f, g, function='linear', smooth=1e-8)
        g_interp = rbf(e_mesh, f_mesh)

# Determine color range from actual data values (not interpolated)
vmin, vmax = g.min(), g.max()

# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 8))

# Create smooth heatmap using pcolormesh with RGB gradient (jet colormap)
# Use the same vmin/vmax as the scatter plot for consistent colors
heatmap = ax.pcolormesh(e_mesh, f_mesh, g_interp, cmap='jet', shading='auto', vmin=vmin, vmax=vmax)

# Overlay the original data points with their EXACT g values colored
# This ensures atomic positions correspond to the real g value from the data file
scatter = ax.scatter(e, f, c=g, cmap='jet', s=150, edgecolors='black', linewidths=2, zorder=10, vmin=vmin, vmax=vmax)

# Add colorbar with height matching the plot
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(heatmap, cax=cax)
cbar.set_label('Distance from plane (Å)', fontsize=12)

# Set axis labels with units
ax.set_xlabel('x (Å)', fontsize=12)
ax.set_ylabel('y (Å)', fontsize=12)

# Add grid
ax.grid(True, alpha=0.3)

# Set equal aspect ratio
ax.set_aspect('equal', adjustable='box')

# Save figure
plt.tight_layout()
plt.savefig('test_projections_heatmap.png', dpi=150, bbox_inches='tight')
print(f"Plot saved to test_projections_heatmap.png")

# Optionally display the plot (comment out if running headless)
# plt.show()
