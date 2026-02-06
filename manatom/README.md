# Manatom - VASP POSCAR Manipulation Tools

A collection of Python scripts for manipulating VASP POSCAR files.

## Scripts

### displace_atoms.py

Displace atomic positions in a VASP POSCAR file along specified directions.

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

- **Compatibility:**
  - Works with both Direct and Cartesian coordinate formats
  - Preserves selective dynamics flags
  - Preserves velocities and additional data
  - Maintains POSCAR formatting and scaling factors

#### Requirements

- Python 3.x
- NumPy

Install requirements:
```bash
pip install numpy
```

#### Usage

```bash
python displace_atoms.py POSCAR --atom <selector> --direction <dir> --amount <value> [--output <filename>]
```

**Arguments:**
- `POSCAR`: Input POSCAR file
- `--atom`: Atom selector (required)
  - Single atom ID: `3` (displaces 3rd atom)
  - Element symbol: `Fe` (displaces all Fe atoms)
  - Multiple atoms: `1,3,5` (displaces atoms 1, 3, and 5)
  - Mixed: `2,4,Fe,O` (displaces atoms 2, 4, and all Fe and O atoms)
- `--direction`: Displacement direction (required)
  - Cartesian: `x`, `y`, or `z`
  - Crystallographic: `a`, `b`, or `c`
- `--amount`: Displacement amount in Ångström (required)
- `--output`: Output filename (optional, default: `POSCAR_disp.vasp`)

#### Examples

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

The script creates a new POSCAR file (default: `POSCAR_disp.vasp`) with the displaced atomic positions. The output file maintains:
- Original comment line
- Scaling factor
- Lattice vectors
- Element names and counts
- Selective dynamics flags (if present)
- Coordinate type (Direct/Cartesian)
- All additional data (velocities, etc.)

#### Notes

- Atom IDs are 1-based (first atom is 1, not 0)
- If an atom is selected multiple times (e.g., `--atom 1,Fe` where atom 1 is Fe), it will only be displaced once
- Negative displacement values are allowed
- The script automatically handles coordinate system conversions

## License

MIT License
