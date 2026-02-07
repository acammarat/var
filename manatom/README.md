# Manatom - VASP POSCAR Manipulation Tools

A collection of Python scripts for manipulating VASP POSCAR files.

## Scripts

### displace_atoms.py

Displace atomic positions in a VASP POSCAR file along specified directions and calculate surface distances.

#### Features

- **Flexible atom selection:**
  - Select individual atoms by ID (1-based index in POSCAR)
  - Select all atoms of a specific element type
  - Select multiple atoms using comma-separated lists
  - Mix atom IDs and element symbols in the same command

- **Displacement directions:**
  - **Cartesian coordinates:** `x`, `y`, `z`
  - **Crystallographic coordinates:** `a`, `b`, `c` (along lattice vectors)

- **Displacement in Angstroms:** All displacements are specified in Ångström units

- **Surface distance calculation:**
  - Calculate the distance between the center of mass of selected atoms and the closest atom surface
  - Surface identification using plane fitting through the three atoms closest to the center of mass
  - Plane optimization to minimize distance to neighboring atoms
  - Configurable neighbor distance threshold

- **Compatibility:**
  - Works with both Direct and Cartesian coordinate formats
  - Preserves selective dynamics flags
  - Preserves velocities and additional data
  - Maintains POSCAR formatting and scaling factors

#### Requirements

- Python 3.x
- NumPy
- SciPy (for surface distance calculation)

Install requirements:
```bash
pip install numpy scipy
```

#### Usage

##### Displacement Mode
```bash
python displace_atoms.py POSCAR --atom <selector> --direction <dir> --amount <value> [--output <filename>]
```

##### Surface Distance Calculation Mode
```bash
python displace_atoms.py POSCAR --atom <selector> --calculate-surface-distance [--neighbor-threshold <value>]
```

**Arguments:**
- `POSCAR`: Input POSCAR file
- `--atom`: Atom selector (required)
  - Single atom ID: `3` (selects 3rd atom)
  - Element symbol: `Fe` (selects all Fe atoms)
  - Multiple atoms: `1,3,5` (selects atoms 1, 3, and 5)
  - Mixed: `2,4,Fe,O` (selects atoms 2, 4, and all Fe and O atoms)
- `--direction`: Displacement direction (required for displacement mode)
  - Cartesian: `x`, `y`, or `z`
  - Crystallographic: `a`, `b`, or `c`
- `--amount`: Displacement amount in Ångström (required for displacement mode)
- `--output`: Output filename (optional, default: `POSCAR_disp.vasp`)
- `--calculate-surface-distance`: Enable surface distance calculation mode
- `--neighbor-threshold`: Distance threshold for finding neighbors in surface calculation (optional, default: 5.0 Å)

#### Examples

##### Displacement Examples

##### Example 1: Displace a single atom by ID
```bash
python displace_atoms.py POSCAR --atom 3 --direction x --amount 0.5
```
Displaces the 3rd atom in the POSCAR by 0.5 Å along the Cartesian x-axis.

##### Example 2: Displace all atoms of a specific element
```bash
python displace_atoms.py POSCAR --atom Fe --direction z --amount -0.3
```
Displaces all Fe atoms by -0.3 Å along the Cartesian z-axis.

##### Example 3: Displace multiple specific atoms
```bash
python displace_atoms.py POSCAR --atom 1,3,5 --direction a --amount 0.2
```
Displaces atoms 1, 3, and 5 by 0.2 Å along the crystallographic a-direction.

##### Example 4: Displace atoms and element types together
```bash
python displace_atoms.py POSCAR --atom 2,4,Fe --direction b --amount 0.1
```
Displaces atoms 2, 4, and all Fe atoms by 0.1 Å along the crystallographic b-direction.

##### Example 5: Displace multiple element types
```bash
python displace_atoms.py POSCAR --atom Ti,O --direction c --amount 0.05
```
Displaces all Ti and O atoms by 0.05 Å along the crystallographic c-direction.

##### Example 6: Custom output filename
```bash
python displace_atoms.py POSCAR --atom 1 --direction y --amount 0.25 --output POSCAR_custom.vasp
```
Displaces atom 1 and saves the result to `POSCAR_custom.vasp`.

##### Surface Distance Calculation Examples

##### Example 7: Calculate surface distance for all Fe atoms
```bash
python displace_atoms.py POSCAR --atom Fe --calculate-surface-distance
```
Calculates the distance from the center of mass of all Fe atoms to the closest surface defined by neighboring atoms.

##### Example 8: Calculate surface distance with custom neighbor threshold
```bash
python displace_atoms.py POSCAR --atom 1,2,3,4,5 --calculate-surface-distance --neighbor-threshold 6.0
```
Calculates the surface distance for atoms 1-5, using a 6.0 Å threshold for finding neighboring atoms.

##### Example 9: Calculate surface distance for multiple element types
```bash
python displace_atoms.py POSCAR --atom Ti,O --calculate-surface-distance --neighbor-threshold 4.5
```
Calculates the surface distance for all Ti and O atoms with a 4.5 Å neighbor threshold.

#### Understanding Surface Distance Calculation

The surface distance calculation feature determines the distance between the center of mass of selected atoms and the closest atom surface using the following algorithm:

1. **Calculate center of mass**: Computes the geometric center of all selected atoms
2. **Find closest atoms**: Identifies the 3 atoms from the selection that are closest to the center of mass
3. **Initial plane fitting**: Creates a plane that passes through these 3 atoms using singular value decomposition (SVD)
4. **Find neighbors**: Identifies atoms within the specified distance threshold that are neighbors of the 3 plane-defining atoms
5. **Plane optimization**: Adjusts the plane position and orientation to minimize the sum of squared distances to the neighboring atoms using the BFGS optimization algorithm
6. **Distance calculation**: Computes the perpendicular distance from the center of mass to the optimized plane

This approach provides a robust estimate of the distance to the atomic surface, taking into account the local atomic environment.

#### Understanding Coordinate Systems

**Cartesian coordinates (x, y, z):**
- Standard 3D Cartesian axes
- Independent of the crystal lattice
- Good for absolute directional displacements

**Crystallographic coordinates (a, b, c):**
- Along the lattice vectors defined in the POSCAR
- Useful for displacements relative to the crystal structure
- The script automatically converts Ångström displacements to fractional coordinates based on lattice vector lengths

#### Output

**Displacement mode:**
The script creates a new POSCAR file (default: `POSCAR_disp.vasp`) with the displaced atomic positions. The output file maintains:
- Original comment line
- Scaling factor
- Lattice vectors
- Element names and counts
- Selective dynamics flags (if present)
- Coordinate type (Direct/Cartesian)
- All additional data (velocities, etc.)

**Surface distance calculation mode:**
The script outputs detailed information to the console including:
- Selected atoms and their indices
- Center of mass coordinates
- The 3 closest atoms used for initial plane fitting
- Initial plane normal and position
- Number and indices of neighboring atoms
- Optimized plane parameters
- Final distance from center of mass to surface (in Ångström)

#### Notes

- Atom IDs are 1-based (first atom is 1, not 0)
- If an atom is selected multiple times (e.g., `--atom 1,Fe` where atom 1 is Fe), it will only be counted once
- Negative displacement values are allowed
- The script automatically handles coordinate system conversions
- For surface distance calculation, at least 3 atoms must be selected
- The neighbor threshold affects the quality of the surface fit - typical values range from 3.0 to 7.0 Å depending on the structure

## License

MIT License
