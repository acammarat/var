#!/usr/bin/env python3
"""
avgpos - Calculate average position and standard deviation of selected atoms
along a crystallographic direction from a POSCAR file.

Copyright (C) 2025
This file is part of phtools.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
"""

import sys
import numpy as np
import argparse


def read_poscar(filename):
    """
    Read a POSCAR file and return atomic structure information.
    
    Parameters:
    -----------
    filename : str
        Path to the POSCAR file
        
    Returns:
    --------
    dict : Dictionary containing structure information
        - 'lattice': 3x3 numpy array with lattice vectors
        - 'elements': list of element symbols
        - 'atom_counts': list of atom counts per element
        - 'positions': Nx3 numpy array with atomic positions
        - 'coordinate_type': 'Direct' or 'Cartesian'
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Read comment line
    comment = lines[0].strip()
    
    # Read scale factor
    scale = float(lines[1].strip())
    if abs(scale - 1.0) > 1e-6:
        print(f"Warning: Scale factor is {scale}, not 1.0. Applying scaling.")
    
    # Read lattice vectors
    lattice = np.zeros((3, 3))
    for i in range(3):
        lattice[i] = [float(x) for x in lines[2+i].split()]
    lattice *= scale
    
    # Read element symbols
    elements = lines[5].split()
    
    # Read atom counts
    atom_counts = [int(x) for x in lines[6].split()]
    total_atoms = sum(atom_counts)
    
    # Check for selective dynamics
    line_idx = 7
    if lines[line_idx].strip()[0].upper() in ['S']:
        line_idx += 1
    
    # Read coordinate type
    coord_type = lines[line_idx].strip()
    coordinate_type = 'Direct' if coord_type[0].upper() in ['D'] else 'Cartesian'
    
    # Read atomic positions
    positions = np.zeros((total_atoms, 3))
    for i in range(total_atoms):
        pos_line = lines[line_idx + 1 + i].split()
        positions[i] = [float(x) for x in pos_line[:3]]
    
    # Convert direct coordinates to Cartesian if needed
    if coordinate_type == 'Direct':
        positions = np.dot(positions, lattice)
    
    return {
        'lattice': lattice,
        'elements': elements,
        'atom_counts': atom_counts,
        'positions': positions,
        'coordinate_type': coordinate_type
    }


def select_atoms(structure, selection):
    """
    Select atoms based on element type or indices.
    
    Parameters:
    -----------
    structure : dict
        Structure dictionary from read_poscar
    selection : str or list
        Either element symbol(s) or atom indices (1-based)
        
    Returns:
    --------
    numpy.ndarray : Indices of selected atoms (0-based)
    """
    if isinstance(selection, str):
        # Selection by element
        selected_indices = []
        idx = 0
        for i, element in enumerate(structure['elements']):
            count = structure['atom_counts'][i]
            if element in selection.split(','):
                selected_indices.extend(range(idx, idx + count))
            idx += count
        return np.array(selected_indices)
    else:
        # Selection by indices (convert from 1-based to 0-based)
        return np.array([i-1 for i in selection])


def get_direction_vector(structure, direction):
    """
    Get the unit vector for the specified crystallographic direction.
    
    Parameters:
    -----------
    structure : dict
        Structure dictionary from read_poscar
    direction : str
        Direction specification: 'x', 'y', 'z', 'a', 'b', 'c', or custom [h,k,l]
        
    Returns:
    --------
    numpy.ndarray : Unit vector in Cartesian coordinates
    """
    lattice = structure['lattice']
    
    if direction.lower() == 'x':
        return np.array([1.0, 0.0, 0.0])
    elif direction.lower() == 'y':
        return np.array([0.0, 1.0, 0.0])
    elif direction.lower() == 'z':
        return np.array([0.0, 0.0, 1.0])
    elif direction.lower() == 'a':
        vec = lattice[0]
    elif direction.lower() == 'b':
        vec = lattice[1]
    elif direction.lower() == 'c':
        vec = lattice[2]
    else:
        # Parse custom direction [h,k,l]
        try:
            direction = direction.strip('[]')
            h, k, l = [float(x) for x in direction.split(',')]
            # Convert Miller indices to Cartesian
            vec = h * lattice[0] + k * lattice[1] + l * lattice[2]
        except:
            raise ValueError(f"Invalid direction specification: {direction}")
    
    # Normalize to unit vector
    return vec / np.linalg.norm(vec)


def calculate_average_position(structure, atom_indices, direction_vector):
    """
    Calculate average position and standard deviation along a direction.
    
    Parameters:
    -----------
    structure : dict
        Structure dictionary from read_poscar
    atom_indices : numpy.ndarray
        Indices of atoms to include in calculation
    direction_vector : numpy.ndarray
        Unit vector defining the direction
        
    Returns:
    --------
    tuple : (average, std_dev, positions_along_dir)
        - average: mean position along the direction
        - std_dev: standard deviation
        - positions_along_dir: array of positions projected onto direction
    """
    positions = structure['positions'][atom_indices]
    
    # Project positions onto the direction vector
    positions_along_dir = np.dot(positions, direction_vector)
    
    # Calculate statistics
    average = np.mean(positions_along_dir)
    std_dev = np.std(positions_along_dir)
    
    return average, std_dev, positions_along_dir


def get_atom_labels(structure, atom_indices, label_format='both'):
    """
    Generate atom labels for selected atoms based on the specified format.
    
    Parameters:
    -----------
    structure : dict
        Structure dictionary from read_poscar
    atom_indices : numpy.ndarray
        Indices of atoms (0-based)
    label_format : str, optional
        Format of labels: 'type' (atomic type only), 'id' (atom ID only), 
        'both' (atomic type + ID, default)
        
    Returns:
    --------
    list : List of atom labels based on format
        - 'type': ['Ti', 'Ti', 'O', 'O']
        - 'id': ['1', '2', '3', '4']
        - 'both': ['Ti1', 'Ti2', 'O3', 'O4']
        where the number corresponds to the atom's position in the POSCAR file (1-based)
    """
    labels = []
    
    # Create a mapping of atom index to element type
    idx = 0
    atom_to_element = []
    for i, element in enumerate(structure['elements']):
        count = structure['atom_counts'][i]
        atom_to_element.extend([element] * count)
        idx += count
    
    # Generate labels based on the specified format
    for atom_idx in atom_indices:
        element = atom_to_element[atom_idx]
        atom_id = atom_idx + 1  # atom_idx is 0-based, so add 1 to get POSCAR file position
        
        if label_format == 'type':
            labels.append(element)
        elif label_format == 'id':
            labels.append(str(atom_id))
        else:  # 'both'
            labels.append(f"{element}{atom_id}")
    
    return labels


def generate_plot_script(data_file, script_file, output_image='heatmap.png', with_labels=False, replicate=(1, 1), no_circles=False, lattice_shifts=None, label_no_box=False, vrange=None, label_at_projection=False):
    """
    Generate a Python script using matplotlib to plot the plane projection data as a heatmap.
    
    Parameters:
    -----------
    data_file : str
        Path to the data file containing projection data
    script_file : str
        Path where the Python plotting script will be written
    output_image : str
        Name of the output image file (default: 'heatmap.png')
    with_labels : bool
        Whether to include atom labels in the plot
    replicate : tuple
        Number of replications along e and f axes (ne, nf)
    no_circles : bool
        Whether to hide the circles representing atom positions (only when labels are not used)
    lattice_shifts : tuple or None
        Tuple of (e_shift, f_shift) representing the periodic shift distances based on lattice vectors.
        If None, uses the data range for shifts (legacy behavior).
    label_no_box : bool
        Whether to show labels without a background box (only when with_labels is True)
    vrange : tuple or None
        Tuple of (vmin, vmax) to specify the color map range. If None, uses data min/max (default).
    label_at_projection : bool
        Whether to position labels at the exact projection coordinates (no offset) and hide circles.
        Only effective when with_labels is True.
    """
    script_content = f"""#!/usr/bin/env python3
\"\"\"
Matplotlib script to visualize plane projection data as a smooth interpolated heatmap.
Generated automatically by avgpos tool.

Usage: python3 {script_file}
\"\"\"

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import Rbf
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Read data from file
data = np.loadtxt('{data_file}', dtype=str)

# Extract coordinates and g values
e_orig = data[:, 0].astype(float)
f_orig = data[:, 1].astype(float)
g_orig = data[:, 2].astype(float)

# Extract labels if present (4th column)
has_labels = data.shape[1] > 3
if has_labels:
    labels_orig = data[:, 3]

# Replication parameters
ne_rep, nf_rep = {replicate[0]}, {replicate[1]}

# Replicate data along e and f axes
e_list, f_list, g_list, labels_list = [], [], [], []

# Use lattice-based shifts if provided, otherwise use data range (legacy)
lattice_shifts = {lattice_shifts}
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

# Determine color range from actual data values (not interpolated) or use custom range
vrange = {vrange}
if vrange is not None:
    vmin, vmax = vrange[0], vrange[1]
else:
    vmin, vmax = g.min(), g.max()

# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 8))

# Create smooth heatmap using pcolormesh with RGB gradient (jet colormap)
# Use the same vmin/vmax as the scatter plot for consistent colors
heatmap = ax.pcolormesh(e_mesh, f_mesh, g_interp, cmap='jet', shading='auto', vmin=vmin, vmax=vmax)
"""
    
    # Add scatter points conditionally
    # Hide circles only if label_at_projection is True AND with_labels is True (labels replace circles)
    if (with_labels or not no_circles) and not (label_at_projection and with_labels):
        script_content += f"""
# Overlay the original data points with their EXACT g values colored
# This ensures atomic positions correspond to the real g value from the data file
scatter = ax.scatter(e, f, c=g, cmap='jet', s=150, edgecolors='black', linewidths=2, zorder=10, vmin=vmin, vmax=vmax)
"""
    
    script_content += f"""
# Add colorbar with height matching the plot
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(heatmap, cax=cax)
cbar.set_label('Distance from plane (Å)', fontsize=12)

# Set axis labels with units
ax.set_xlabel('x (Å)', fontsize=12)
ax.set_ylabel('y (Å)', fontsize=12)
"""
    
    if with_labels:
        if label_at_projection:
            # Position labels at exact projection coordinates (no offset)
            if label_no_box:
                script_content += f"""
# Add atom labels at exact projection coordinates without background box
if has_labels:
    for i in range(len(e)):
        ax.annotate(labels[i], (e[i], f[i]), 
                    ha='center', va='center',
                    fontsize=10, fontweight='bold', color='black')
"""
            else:
                script_content += f"""
# Add atom labels at exact projection coordinates
if has_labels:
    for i in range(len(e)):
        ax.annotate(labels[i], (e[i], f[i]), 
                    ha='center', va='center',
                    fontsize=10, fontweight='bold', color='black',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', alpha=0.7))
"""
        elif label_no_box:
            script_content += f"""
# Add atom labels without background box (already replicated above if needed)
if has_labels:
    for i in range(len(e)):
        ax.annotate(labels[i], (e[i], f[i]), 
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=10, fontweight='bold', color='black')
"""
        else:
            script_content += f"""
# Add atom labels (already replicated above if needed)
if has_labels:
    for i in range(len(e)):
        ax.annotate(labels[i], (e[i], f[i]), 
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=10, fontweight='bold', color='black',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', alpha=0.7))
"""
    
    script_content += f"""
# Add grid
ax.grid(True, alpha=0.3)

# Set equal aspect ratio
ax.set_aspect('equal', adjustable='box')

# Save figure
plt.tight_layout()
plt.savefig('{output_image}', dpi=150, bbox_inches='tight')
print(f"Plot saved to {output_image}")

# Optionally display the plot (comment out if running headless)
# plt.show()
"""
    
    with open(script_file, 'w') as f:
        f.write(script_content)
    
    # Make the script executable (Unix-like systems)
    import os
    try:
        os.chmod(script_file, 0o755)
    except:
        pass  # Ignore if chmod fails (e.g., on Windows)


def calculate_plane_projections(structure, atom_indices, direction_vector, average_position, flip_g=False):
    """
    Calculate orthogonal projections of atoms onto a plane perpendicular to 
    the direction vector and passing through the average position.
    
    Parameters:
    -----------
    structure : dict
        Structure dictionary from read_poscar
    atom_indices : numpy.ndarray
        Indices of atoms to include in calculation
    direction_vector : numpy.ndarray
        Unit vector defining the direction (normal to the plane)
    average_position : float
        Average position along the direction (defines plane location)
    flip_g : bool, optional
        If True, flip the sign of g (default: False)
        
    Returns:
    --------
    tuple : (projections, basis1, basis2)
        - projections: Nx3 array where each row contains [e, f, g]
            - e, f: 2D coordinates of the projection on the plane
            - g: average_position minus the distance of the atom from the plane
                 (or the opposite if flip_g=True)
        - basis1: First basis vector of the plane (unit vector)
        - basis2: Second basis vector of the plane (unit vector)
    """
    positions = structure['positions'][atom_indices]
    
    # Calculate the distance of each atom along the direction vector
    distances_along_dir = np.dot(positions, direction_vector)
    
    # Calculate the signed distance from each atom to the plane
    # (positive if atom is on the side of the direction vector, negative otherwise)
    signed_distances = distances_along_dir - average_position
    
    # Project each atom onto the plane
    # projection = position - (signed_distance * normal_vector)
    projections_3d = positions - np.outer(signed_distances, direction_vector)
    
    # Create an orthonormal basis for the plane
    # Find two orthogonal vectors in the plane
    # Start with an arbitrary vector not parallel to direction_vector
    if abs(direction_vector[2]) < 0.9:
        arbitrary = np.array([0.0, 0.0, 1.0])
    else:
        arbitrary = np.array([1.0, 0.0, 0.0])
    
    # First basis vector in the plane (orthogonal to direction_vector)
    basis1 = arbitrary - np.dot(arbitrary, direction_vector) * direction_vector
    basis1 = basis1 / np.linalg.norm(basis1)
    
    # Second basis vector (orthogonal to both direction_vector and basis1)
    basis2 = np.cross(direction_vector, basis1)
    basis2 = basis2 / np.linalg.norm(basis2)
    
    # Project the 3D projections onto the 2D plane coordinate system
    e_coords = np.dot(projections_3d, basis1)
    f_coords = np.dot(projections_3d, basis2)
    
    # Calculate g = average_position - distance_from_plane
    # The distance from the plane is the absolute value of signed_distances
    # But we want: average_position - distance_of_atom_from_plane
    # Since distance_along_dir = average_position + signed_distance
    # we have: g = average_position - abs(signed_distance) if we want actual distance
    # But the requirement says: "average_position minus the distance of the atom from the plane"
    # which could mean: average_position - distance_along_dir = -signed_distances
    if flip_g:
        g_coords = signed_distances
    else:
        g_coords = -signed_distances
    
    # Combine into Nx3 array
    result = np.column_stack((e_coords, f_coords, g_coords))
    
    return result, basis1, basis2


def extract_profile_along_x(projections, y_value, num_points=200):
    """
    Extract a 1D profile along x (e-axis) at a specified y (f-axis) coordinate.
    Uses RBF interpolation to get g values along the line.
    
    Parameters:
    -----------
    projections : numpy.ndarray
        Nx3 array where each row contains [e, f, g]
    y_value : float
        The y (f-axis) coordinate at which to extract the profile
    num_points : int, optional
        Number of points to sample along x (default: 200)
        
    Returns:
    --------
    tuple : (x_coords, g_values)
        - x_coords: array of x (e-axis) coordinates
        - g_values: array of interpolated g values at each x coordinate
    """
    from scipy.interpolate import Rbf
    
    e_coords = projections[:, 0]
    f_coords = projections[:, 1]
    g_coords = projections[:, 2]
    
    # Determine the range of x coordinates
    e_min, e_max = e_coords.min(), e_coords.max()
    
    # Create x coordinates along the profile line
    x_profile = np.linspace(e_min, e_max, num_points)
    
    # Create y coordinates (all at the specified y_value)
    y_profile = np.full(num_points, y_value)
    
    # Use Radial Basis Function interpolation to get g values
    # Try thin_plate first, fall back to multiquadric, then linear
    try:
        rbf = Rbf(e_coords, f_coords, g_coords, function='thin_plate', smooth=1e-10)
        g_profile = rbf(x_profile, y_profile)
    except:
        try:
            rbf = Rbf(e_coords, f_coords, g_coords, function='multiquadric', smooth=1e-10)
            g_profile = rbf(x_profile, y_profile)
        except:
            rbf = Rbf(e_coords, f_coords, g_coords, function='linear', smooth=1e-8)
            g_profile = rbf(x_profile, y_profile)
    
    return x_profile, g_profile


def main():
    parser = argparse.ArgumentParser(
        description='Calculate average position and standard deviation of selected atoms '
                    'along a crystallographic direction from a POSCAR file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Average position of all Se atoms along z-axis
  %(prog)s POSCAR -s Se -d z
  
  # Average position of atoms 2,3,4 along c lattice vector
  %(prog)s POSCAR -i 2,3,4 -d c
  
  # Average position along custom direction [1,1,0]
  %(prog)s POSCAR -s W,Mo -d [1,1,0]
        """
    )
    
    parser.add_argument('poscar', help='Path to POSCAR file')
    parser.add_argument('-s', '--select', type=str, 
                        help='Select atoms by element symbol(s), comma-separated (e.g., "Se" or "W,Mo")')
    parser.add_argument('-i', '--indices', type=str,
                        help='Select atoms by indices (1-based), comma-separated (e.g., "1,2,3")')
    parser.add_argument('-d', '--direction', type=str, required=True,
                        help='Direction: x, y, z (Cartesian) or a, b, c (lattice vectors) '
                             'or [h,k,l] (Miller indices)')
    parser.add_argument('-o', '--output', type=str,
                        help='Output file for plane projection data (3 columns: e, f, g)')
    parser.add_argument('--plot', action='store_true',
                        help='Generate Python matplotlib script for heatmap visualization (requires -o)')
    parser.add_argument('--labels', type=str, nargs='?', const='both', default=None, 
                        choices=['type', 'id', 'both'],
                        help='Include atom labels in output and plot (requires -o). '
                             'Options: "type" (atomic type only), "id" (atom ID only), '
                             '"both" (atomic type + ID, default if flag is used without value)')
    parser.add_argument('--replicate', type=str, default='1,1',
                        help='Replicate the plot along e and f axes (format: "ne,nf", e.g., "2.5,3" for 2.5x3 replication)')
    parser.add_argument('--no-circles', action='store_true',
                        help='Hide circles representing atom positions in the plot (only when --labels is not used)')
    parser.add_argument('--label-no-box', action='store_true',
                        help='Show labels without background box in the plot (only when --labels is used)')
    parser.add_argument('--label-at-projection', action='store_true',
                        help='Position labels at the exact atom projection coordinates instead of offset from circles. '
                             'When used, circles are hidden and labels are centered at the projection position. '
                             'Only effective when --labels is used. Requires both -o and --plot.')
    parser.add_argument('--vrange', type=str, default=None,
                        help='Specify color map range as "vmin,vmax" (e.g., --vrange=-5,5 or --vrange=0,10). '
                             'If not specified, uses data min/max. Requires both -o and --plot. '
                             'Useful for comparing multiple plots with consistent color scales.')
    parser.add_argument('--flip-g', action='store_true',
                        help='Flip the sign of g values in the plane projection output. '
                             'By default, g = average_position - distance_along_direction. '
                             'With this flag, g = distance_along_direction - average_position.')
    parser.add_argument('--profile-y', type=float, default=None,
                        help='Extract 1D profile along x (e-axis) at the specified y (f-axis) coordinate. '
                             'Creates a data file with x and g values. Requires both -o and --plot.')
    
    args = parser.parse_args()
    
    # Validate input
    if not args.select and not args.indices:
        parser.error("Must specify either --select or --indices")
    if args.select and args.indices:
        parser.error("Cannot specify both --select and --indices")
    
    # Read POSCAR file
    print(f"Reading POSCAR file: {args.poscar}")
    try:
        structure = read_poscar(args.poscar)
    except Exception as e:
        print(f"Error reading POSCAR file: {e}")
        sys.exit(1)
    
    # Print structure information
    total_atoms = sum(structure['atom_counts'])
    print(f"Structure contains {total_atoms} atoms:")
    for elem, count in zip(structure['elements'], structure['atom_counts']):
        print(f"  {elem}: {count}")
    print()
    
    # Select atoms
    if args.select:
        atom_indices = select_atoms(structure, args.select)
        print(f"Selected {len(atom_indices)} atom(s) of type: {args.select}")
    else:
        indices = [int(x) for x in args.indices.split(',')]
        atom_indices = select_atoms(structure, indices)
        print(f"Selected {len(atom_indices)} atom(s) by indices: {args.indices}")
    
    if len(atom_indices) == 0:
        print("Error: No atoms selected!")
        sys.exit(1)
    
    # Get direction vector
    try:
        direction_vector = get_direction_vector(structure, args.direction)
        print(f"Direction vector (Cartesian): [{direction_vector[0]:.6f}, "
              f"{direction_vector[1]:.6f}, {direction_vector[2]:.6f}]")
    except Exception as e:
        print(f"Error parsing direction: {e}")
        sys.exit(1)
    
    # Calculate average position
    average, std_dev, positions = calculate_average_position(
        structure, atom_indices, direction_vector
    )
    
    # Print results
    print()
    print("=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"Number of atoms: {len(atom_indices)}")
    print(f"Average position: {average:.6f} Å")
    print(f"Standard deviation: {std_dev:.6f} Å")
    print()
    
    # Print individual positions if not too many
    if len(atom_indices) <= 20:
        print("Individual positions along direction:")
        for i, (idx, pos) in enumerate(zip(atom_indices, positions)):
            print(f"  Atom {idx+1}: {pos:.6f} Å")
    
    # Calculate and write plane projections if output file is specified
    if args.output:
        projections, basis1, basis2 = calculate_plane_projections(
            structure, atom_indices, direction_vector, average, flip_g=args.flip_g
        )
        
        # Get atom labels if requested
        g_description = 'distance_from_plane - average_position' if args.flip_g else 'average_position - distance_from_plane'
        if args.labels:
            labels = get_atom_labels(structure, atom_indices, args.labels)
            # Create a combined array with projections and labels
            # Save with labels as 4th column
            with open(args.output, 'w') as f:
                f.write('# e f g label\n')
                f.write('# Projections onto plane perpendicular to direction vector\n')
                f.write('# e, f: 2D coordinates on plane\n')
                f.write(f'# g: {g_description}\n')
                f.write('# label: atom type and ID (e.g., Ti1, O2)\n')
                for i, label in enumerate(labels):
                    f.write(f"{projections[i, 0]:.6f} {projections[i, 1]:.6f} {projections[i, 2]:.6f} {label}\n")
            
            print()
            print(f"Plane projection data with labels written to: {args.output}")
            print(f"  Columns: e, f, g, label")
        else:
            # Write to file without labels
            np.savetxt(args.output, projections, fmt='%.6f', 
                       header=f'e f g\nProjections onto plane perpendicular to direction vector\n'
                              f'e, f: 2D coordinates on plane\n'
                              f'g: {g_description}',
                       comments='# ')
            
            print()
            print(f"Plane projection data written to: {args.output}")
            print(f"  Columns: e, f, g")
        
        print(f"  e, f: 2D coordinates of atom projection on plane")
        if args.flip_g:
            print(f"  g: signed distance from plane (atom_distance - average_position) [flipped]")
        else:
            print(f"  g: signed distance from plane (average_position - atom_distance)")
        
        # Generate matplotlib plot script if requested
        if args.plot:
            import os
            # Determine output paths
            base_name = os.path.splitext(args.output)[0]
            script_file = f"{base_name}_plot.py"
            image_file = f"{base_name}_heatmap.png"
            
            # Parse replication argument
            try:
                replicate_parts = args.replicate.split(',')
                ne_rep = float(replicate_parts[0])
                nf_rep = float(replicate_parts[1]) if len(replicate_parts) > 1 else ne_rep
                replicate = (ne_rep, nf_rep)
            except (ValueError, IndexError):
                print(f"Warning: Invalid replicate format '{args.replicate}'. Using default 1,1")
                replicate = (1, 1)
            
            # Parse vrange argument
            vrange = None
            if args.vrange:
                try:
                    vrange_parts = args.vrange.split(',')
                    if len(vrange_parts) != 2:
                        raise ValueError("vrange must have exactly 2 values")
                    vmin = float(vrange_parts[0])
                    vmax = float(vrange_parts[1])
                    if vmin >= vmax:
                        raise ValueError("vmin must be less than vmax")
                    vrange = (vmin, vmax)
                except (ValueError, IndexError) as e:
                    print(f"Warning: Invalid vrange format '{args.vrange}': {e}. Using data min/max")
                    vrange = None
            
            # Calculate lattice-based shift distances
            # Project each lattice vector onto the plane basis to find which provides the proper periodicity
            # For a general plane, we need to find the lattice vectors that give the minimal periodic shifts
            lattice = structure['lattice']
            
            # Project all three lattice vectors onto the plane basis
            lattice_proj_e = np.array([np.dot(lattice[i], basis1) for i in range(3)])
            lattice_proj_f = np.array([np.dot(lattice[i], basis2) for i in range(3)])
            
            # Find the lattice vector(s) that provide the minimal periodic shift along e and f
            # We want non-zero projections that represent the periodicity in the plane
            # The proper shift is the smallest non-zero projection magnitude for each axis
            
            # For e-axis: find the smallest non-zero absolute projection
            e_shifts = [abs(proj) for proj in lattice_proj_e if abs(proj) > 1e-6]
            e_shift = min(e_shifts) if e_shifts else (projections[:, 0].max() - projections[:, 0].min())
            
            # For f-axis: find the smallest non-zero absolute projection
            f_shifts = [abs(proj) for proj in lattice_proj_f if abs(proj) > 1e-6]
            f_shift = min(f_shifts) if f_shifts else (projections[:, 1].max() - projections[:, 1].min())
            
            lattice_shifts = (e_shift, f_shift)
            
            # Generate the plotting script
            generate_plot_script(args.output, script_file, image_file, args.labels, replicate, args.no_circles, lattice_shifts, args.label_no_box, vrange, args.label_at_projection)
            
            print()
            print(f"Matplotlib plotting script generated: {script_file}")
            print(f"To create the heatmap, run: python3 {script_file}")
            print(f"Output image will be: {image_file}")
            if vrange:
                print(f"  (with custom color range: {vrange[0]} to {vrange[1]})")
            if args.labels:
                if args.label_at_projection:
                    label_style = " at projection coordinates (no circles)" + (" without box" if args.label_no_box else "")
                else:
                    label_style = " without box" if args.label_no_box else ""
                print(f"  (with atom labels{label_style})")
            if replicate != (1, 1):
                print(f"  (with {replicate[0]}x{replicate[1]} replication)")
            if args.no_circles and not args.labels:
                print(f"  (without atom position circles)")
            
            # Extract profile along x at specified y if requested
            if args.profile_y is not None:
                profile_file = f"{base_name}_profile.dat"
                
                # Extract the profile
                x_profile, g_profile = extract_profile_along_x(projections, args.profile_y)
                
                # Save profile to file
                with open(profile_file, 'w') as f:
                    f.write('# x g\n')
                    f.write('# Profile along x (e-axis) at specified y (f-axis) coordinate\n')
                    f.write(f'# y coordinate: {args.profile_y:.6f}\n')
                    f.write('# x: coordinate along e-axis\n')
                    f.write('# g: interpolated g value from plane projection\n')
                    for x, g in zip(x_profile, g_profile):
                        f.write(f"{x:.6f} {g:.6f}\n")
                
                print()
                print(f"Profile data written to: {profile_file}")
                print(f"  Profile extracted along x at y = {args.profile_y:.6f}")
                print(f"  Columns: x, g")
        else:
            # --plot not specified but -o is specified
            if args.profile_y is not None:
                print()
                print("Warning: --profile-y flag requires both -o/--output and --plot to be specified. Ignoring.")
    elif args.plot or args.labels or args.vrange or args.label_at_projection or args.profile_y is not None:
        print()
        if args.plot or args.labels:
            print("Warning: --plot and --labels flags require -o/--output to be specified. Ignoring.")
        if args.vrange:
            print("Warning: --vrange flag requires both -o/--output and --plot to be specified. Ignoring.")
        if args.label_at_projection:
            print("Warning: --label-at-projection flag requires -o/--output, --plot, and --labels to be specified. Ignoring.")
        if args.profile_y is not None:
            print("Warning: --profile-y flag requires both -o/--output and --plot to be specified. Ignoring.")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
