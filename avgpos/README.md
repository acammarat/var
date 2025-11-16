# avgpos

A tool to calculate the average position and standard deviation of selected atoms along a specified crystallographic direction from a POSCAR file.

## Features

- Read POSCAR files (VASP structure format)
- Select atoms by element type or by indices
- Calculate average position along:
  - Cartesian directions (x, y, z)
  - Crystallographic lattice vectors (a, b, c)
  - Custom Miller indices [h,k,l]
- Calculate standard deviation of positions
- Calculate and export orthogonal projections onto a plane perpendicular to the direction vector and passing through the average position

## Requirements

- Python 3.6 or higher
- NumPy
- Matplotlib (for generating heatmap visualizations)
- SciPy (for interpolation in smooth heatmaps)

## Installation

No installation required. Simply make the script executable:

```bash
chmod +x avgpos.py
```

Or run it with Python:

```bash
python3 avgpos.py
```

## Usage

### Basic syntax

```bash
./avgpos.py POSCAR -s <elements> -d <direction>
./avgpos.py POSCAR -i <indices> -d <direction>
```

### Options

- `POSCAR`: Path to the POSCAR file (required)
- `-s, --select`: Select atoms by element symbol(s), comma-separated (e.g., "Se" or "W,Mo")
- `-i, --indices`: Select atoms by indices (1-based), comma-separated (e.g., "1,2,3")
- `-d, --direction`: Direction specification (required):
  - Cartesian: `x`, `y`, `z`
  - Lattice vectors: `a`, `b`, `c`
  - Miller indices: `[h,k,l]` (e.g., `[1,1,0]`)
- `-o, --output`: Output file for plane projection data (optional)
- `--plot`: Generate Python matplotlib script for heatmap visualization (requires `-o`)
- `--labels [type|id|both]`: Include atom labels in output and plot (requires `-o`)
  - `type`: Show only atomic type (e.g., Se, Ti)
  - `id`: Show only atom ID from POSCAR file (e.g., 2, 4)
  - `both` or no value: Show both type and ID (e.g., Se2, Ti4) - default
- `--replicate`: Replicate the plot along e and f axes (format: "ne,nf", default: "1,1")
  - Supports non-integer replication (e.g., "2.5,3" for 2.5x3 replication)
- `--no-circles`: Hide circles representing atom positions in the plot (only effective when `--labels` is not used)
- `--label-no-box`: Show labels without background box in the plot (only when `--labels` is used)
- `--label-at-projection`: Position labels at the exact atom projection coordinates instead of offset from circles. When used, circles are hidden and labels are centered at the projection position. Only effective when `--labels` is used. Requires both `-o` and `--plot`.
- `--vrange`: Specify color map range as "vmin,vmax" (e.g., `--vrange=-5,5` or `--vrange=0,10`). If not specified, uses data min/max. Requires both `-o` and `--plot`.
- `--flip-g`: Flip the sign of g values in the plane projection output. By default, g = average_position - distance_along_direction. With this flag, g = distance_along_direction - average_position.
- `--profile-y`: Extract 1D profile along x (e-axis) at the specified y (f-axis) coordinate. Creates a data file with x and g values. Requires both `-o` and `--plot`.

### Examples

Calculate average position of all Se atoms along the z-axis:
```bash
./avgpos.py POSCAR -s Se -d z
```

Calculate average position of atoms 2, 3, and 4 along the c lattice vector:
```bash
./avgpos.py POSCAR -i 2,3,4 -d c
```

Calculate average position of W and Mo atoms along the [1,1,0] direction:
```bash
./avgpos.py POSCAR -s W,Mo -d [1,1,0]
```

Calculate average position of all atoms of multiple elements along x-axis:
```bash
./avgpos.py POSCAR -s Se,Mo -d x
```

Calculate average position and export plane projection data:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat
```

Calculate average position and generate matplotlib heatmap script:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot
# Then run: python3 projections_plot.py
```

Calculate average position with atom labels and generate labeled heatmap:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot --labels
# Then run: python3 projections_plot.py
# Labels will show atom type and POSCAR file ID (e.g., Se2, Se3, Ti1)
```

Use different label formats:
```bash
# Show only atomic type in labels
./avgpos.py POSCAR -s Se -d z -o projections.dat --labels type

# Show only atom IDs in labels  
./avgpos.py POSCAR -s Se -d z -o projections.dat --labels id

# Show both (same as --labels without argument)
./avgpos.py POSCAR -s Se -d z -o projections.dat --labels both
```

Generate heatmap with 2.5x3 replication along e and f axes:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot --replicate 2.5,3
# Then run: python3 projections_plot.py
# Plot will show 2.5 replications along e-axis and 3 along f-axis
```

Generate smooth heatmap without atom position circles:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot --no-circles
# Then run: python3 projections_plot.py
# Plot will show only the smooth interpolated surface without circles
```

Generate heatmap with labels positioned at projection coordinates:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot --labels --label-at-projection
# Then run: python3 projections_plot.py
# Labels will be centered at the exact projection coordinates, replacing the circles
# Useful for cleaner plots where atom identities are more important than circles
```

Generate heatmap with labels at projection coordinates without background box:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot --labels --label-at-projection --label-no-box
# Then run: python3 projections_plot.py
# Labels will be centered at projection coordinates without background boxes
```

Generate heatmap with custom color map range:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot --vrange=-5,5
# Then run: python3 projections_plot.py
# Color map will range from -5 to 5 instead of using data min/max
# Useful for comparing multiple plots with consistent color scales
```

Generate heatmap with positive color range:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot --vrange=0,10
# Then run: python3 projections_plot.py
# Color map will range from 0 to 10
```

Calculate average position with flipped g values:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --flip-g
# g values will be calculated as: distance_along_direction - average_position
# instead of the default: average_position - distance_along_direction
```

Extract 1D profile along x at a specific y coordinate:
```bash
./avgpos.py POSCAR -s Se -d z -o projections.dat --plot --profile-y=0.5
# Then run: python3 projections_plot.py
# Creates projections_profile.dat with x and g values at y=0.5
# Useful for analyzing variation along a specific line in the plane
```

## Output

### Standard Output

The tool displays:
- Structure information (number and types of atoms)
- Selected atoms
- Direction vector in Cartesian coordinates
- Average position along the direction (in Ångströms)
- Standard deviation (in Ångströms)
- Individual atomic positions (if 20 or fewer atoms are selected)

### Plane Projection Output File (optional)

When the `-o` option is specified, the tool generates a data file containing:
- **Column 1 (e)**: First coordinate of the atom's orthogonal projection onto the plane
- **Column 2 (f)**: Second coordinate of the atom's orthogonal projection onto the plane  
- **Column 3 (g)**: Signed distance from the plane (average_position - atom_distance_along_direction by default, or atom_distance_along_direction - average_position with `--flip-g`)
- **Column 4 (label)**: Atom label - only when `--labels` is used. Format depends on the label option:
  - `--labels type`: Atomic type only (e.g., Ti, O)
  - `--labels id`: POSCAR atom ID only (e.g., 1, 2)
  - `--labels both`: Type and ID (e.g., Ti1, O2)

The plane is perpendicular to the specified direction vector and passes through the calculated average position. The e and f coordinates form an orthonormal 2D coordinate system in the plane.

### Matplotlib Plotting Script (optional)

When the `--plot` flag is used along with `-o`, the tool generates a Python script using matplotlib that creates a smooth interpolated heatmap visualization of the plane projection data:
- **Script file**: Named as `<output_basename>_plot.py`
- **Output image**: Named as `<output_basename>_heatmap.png`
- The heatmap uses Radial Basis Function (RBF) interpolation to create a smooth surface covering the entire e,f range
- g values are represented by a color gradient (coolwarm colormap)
- Original data points are overlaid as black dots for reference
- When `--labels` is also used, atom labels (element+POSCAR file ID) are annotated on the plot
- To generate the plot, run: `python3 <script_file>`

**Requirements**: Matplotlib and SciPy must be installed on your system to generate the visualization.

### Profile Data File (optional)

When the `--profile-y` flag is used along with `-o` and `--plot`, the tool generates a 1D profile data file containing:
- **Profile file**: Named as `<output_basename>_profile.dat`
- **Column 1 (x)**: x coordinate along the e-axis
- **Column 2 (g)**: Interpolated g value at the specified y coordinate
- The profile is extracted using the same RBF interpolation as the heatmap visualization
- This is useful for analyzing variation along a specific line in the plane projection

## Example Output

```
Reading POSCAR file: POSCAR
Structure contains 6 atoms:
  W: 1
  Se: 4
  Mo: 1

Selected 4 atom(s) of type: Se
Direction vector (Cartesian): [0.000000, 0.000000, 1.000000]

============================================================
RESULTS
============================================================
Number of atoms: 4
Average position: 35.123456 Å
Standard deviation: 2.345678 Å

Individual positions along direction:
  Atom 2: 38.343795 Å
  Atom 3: 28.534648 Å
  Atom 4: 31.870770 Å
  Atom 5: 35.012611 Å
```

## POSCAR File Format

The tool supports standard VASP POSCAR format with:
- Comment line
- Scale factor (preferably 1.0)
- Lattice vectors (3 lines)
- Element symbols
- Atom counts per element
- Coordinate type (Direct/Cartesian)
- Atomic positions

Both Direct (fractional) and Cartesian coordinates are supported.

## License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

## Author

Part of the phtools collection: https://github.com/acammarat/phtools
