#!/usr/bin/env python3
"""
Matplotlib script to visualize plane projection data as a smooth interpolated heatmap.
Generated automatically by avgpos tool.

Usage: python3 mod135m_s1_plot.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import Rbf
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Read data from file
data = np.loadtxt('mod135m_s1.dat', dtype=str)

# Extract coordinates and g values
e_orig = data[:, 0].astype(float)
f_orig = data[:, 1].astype(float)
g_orig = data[:, 2].astype(float)

# Extract labels if present (4th column)
has_labels = data.shape[1] > 3
if has_labels:
    labels_orig = data[:, 3]

# Replication parameters
ne_rep, nf_rep = 6.0, 6.0

# Replicate data along e and f axes
e_list, f_list, g_list, labels_list = [], [], [], []

# Use lattice-based shifts if provided, otherwise use data range (legacy)
lattice_shifts = (np.float64(9.45996171912), np.float64(3.7256864482))
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
# Always interpolate on the full data range
e_data_min, e_data_max = e.min(), e.max()
f_data_min, f_data_max = f.min(), f.max()

# Create grid for full data range
e_grid_full = np.linspace(e_data_min, e_data_max, 200)
f_grid_full = np.linspace(f_data_min, f_data_max, 200)
e_mesh_full, f_mesh_full = np.meshgrid(e_grid_full, f_grid_full)

# Use Radial Basis Function interpolation with very small smoothing
# This ensures interpolation passes extremely close to data points while handling duplicates
# smooth value is set to a tiny value to get nearly exact values at data points
try:
    rbf = Rbf(e, f, g, function='thin_plate', smooth=1e-10)
    g_interp_full = rbf(e_mesh_full, f_mesh_full)
except:
    # Fall back to multiquadric if thin_plate fails
    try:
        rbf = Rbf(e, f, g, function='multiquadric', smooth=1e-10)
        g_interp_full = rbf(e_mesh_full, f_mesh_full)
    except:
        # Last resort: use linear with small smoothing
        rbf = Rbf(e, f, g, function='linear', smooth=1e-8)
        g_interp_full = rbf(e_mesh_full, f_mesh_full)

# Determine the range to display (either erange or full data range)
erange = (10.0, 30.0, 5.0, 10.0)
if erange is not None:
    e_min, e_max, f_min, f_max = erange[0], erange[1], erange[2], erange[3]
    # Calculate shift amounts to make labels start from 0
    e_display_shift = e_min
    f_display_shift = f_min
    
    # Extract the portion of the interpolated data within erange
    # Find indices that fall within the erange window
    e_mask = (e_grid_full >= e_min) & (e_grid_full <= e_max)
    f_mask = (f_grid_full >= f_min) & (f_grid_full <= f_max)
    
    # Extract the sub-grid
    e_grid = e_grid_full[e_mask]
    f_grid = f_grid_full[f_mask]
    
    # Create meshgrid for the extracted region
    e_mesh, f_mesh = np.meshgrid(e_grid, f_grid)
    
    # Extract the corresponding portion of g_interp
    # We need to slice the 2D array using the masks
    g_interp = g_interp_full[np.ix_(f_mask, e_mask)]
else:
    e_min, e_max = e_data_min, e_data_max
    f_min, f_max = f_data_min, f_data_max
    # No shift when erange is not specified
    e_display_shift = 0
    f_display_shift = 0
    
    # Use the full interpolation
    e_mesh = e_mesh_full
    f_mesh = f_mesh_full
    g_interp = g_interp_full

# Apply coordinate shift for display (labels start from 0 when erange is used)
e_mesh_display = e_mesh - e_display_shift
f_mesh_display = f_mesh - f_display_shift
e_display = e - e_display_shift
f_display = f - f_display_shift

# Determine color range from actual data values (not interpolated) or use custom range
vrange = (-0.3, 0.3)
if vrange is not None:
    vmin, vmax = vrange[0], vrange[1]
else:
    vmin, vmax = g.min(), g.max()

# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 8))

# Create smooth heatmap using pcolormesh with RGB gradient (jet colormap)
# Use the same vmin/vmax as the scatter plot for consistent colors
heatmap = ax.pcolormesh(e_mesh_display, f_mesh_display, g_interp, cmap='jet', shading='auto', vmin=vmin, vmax=vmax)

# Add colorbar with height matching the plot
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(heatmap, cax=cax)
cbar.set_label('Distance from plane (Å)', fontsize=12)

# Set axis labels with units
ax.set_xlabel('x (Å)', fontsize=12)
ax.set_ylabel('y (Å)', fontsize=12)

# Set axis limits to match the display range (shifted to start from 0 when erange is used)
ax.set_xlim(e_min - e_display_shift, e_max - e_display_shift)
ax.set_ylim(f_min - f_display_shift, f_max - f_display_shift)

# Add atom labels at exact projection coordinates without background box
if has_labels:
    for i in range(len(e)):
        ax.annotate(labels[i], (e_display[i], f_display[i]), 
                    ha='center', va='center',
                    fontsize=10, fontweight='bold', color='black')

# Add grid
ax.grid(True, alpha=0.3)

# Set equal aspect ratio
ax.set_aspect('equal', adjustable='box')

# Save figure
plt.tight_layout()
plt.savefig('mod135m_s1_heatmap.png', dpi=150, bbox_inches='tight')
print(f"Plot saved to mod135m_s1_heatmap.png")

# Optionally display the plot (comment out if running headless)
# plt.show()

# Write ASCII matrix file for Gwyddion
# This uses the same interpolated data (g_interp) that was used for the heatmap
gwyddion_file = 'mod135m_s1_matrix.dat'
print(f"Writing ASCII matrix to {gwyddion_file}...")

# Get the grid dimensions
yres, xres = g_interp.shape

# Write ASCII data matrix file
with open(gwyddion_file, 'w') as gwy_f:
    # Write header with metadata
    gwy_f.write(f"# Plane projection data from avgpos - ASCII matrix format\n")
    gwy_f.write(f"# Grid resolution: {xres} x {yres}\n")
    gwy_f.write(f"# Original data X range (e): {e_min:.10e} to {e_max:.10e} (width: {e_max - e_min:.10e} Å)\n")
    gwy_f.write(f"# Original data Y range (f): {f_min:.10e} to {f_max:.10e} (width: {f_max - f_min:.10e} Å)\n")
    gwy_f.write(f"# Display X range (e): {e_min - e_display_shift:.10e} to {e_max - e_display_shift:.10e} (shifted by {e_display_shift:.10e} Å)\n")
    gwy_f.write(f"# Display Y range (f): {f_min - f_display_shift:.10e} to {f_max - f_display_shift:.10e} (shifted by {f_display_shift:.10e} Å)\n")
    gwy_f.write(f"# Z values (g): signed distance from plane (Å)\n")
    gwy_f.write(f"# Data format: {yres} rows x {xres} columns\n")
    gwy_f.write(f"# Interpolation: RBF (thin_plate, smooth=1e-10)\n")
    gwy_f.write(f"#\n")
    
    # Write data matrix in row-major order
    # Each row of the matrix is one line in the file
    for i in range(yres):
        row_values = [f"{g_interp[i, j]:.10e}" for j in range(xres)]
        gwy_f.write(" ".join(row_values) + "\n")

print(f"ASCII matrix written to {gwyddion_file}")
