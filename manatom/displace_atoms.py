#!/usr/bin/env python3
"""
Script to displace atomic positions in a VASP POSCAR file.

Usage:
    python displace_atoms.py POSCAR --atom <atom_id_or_symbol> --direction <x|y|z|a|b|c> --amount <displacement_in_angstrom>

Examples:
    # Displace atom ID 3 (third atom in POSCAR) along cartesian x by 0.5 Angstrom
    python displace_atoms.py POSCAR --atom 3 --direction x --amount 0.5
    
    # Displace all Fe atoms along crystallographic a by 0.2 Angstrom
    python displace_atoms.py POSCAR --atom Fe --direction a --amount 0.2
    
    # Displace atoms 1, 3, and 5 along cartesian z by 0.3 Angstrom
    python displace_atoms.py POSCAR --atom 1,3,5 --direction z --amount 0.3
    
    # Displace atoms 2, 4 and all Fe atoms along crystallographic b by 0.1 Angstrom
    python displace_atoms.py POSCAR --atom 2,4,Fe --direction b --amount 0.1
"""

import argparse
import numpy as np
import sys


class POSCAR:
    def __init__(self, filename):
        """Read and parse a POSCAR file."""
        self.filename = filename
        self.read_poscar()
    
    def read_poscar(self):
        """Read POSCAR file and store its contents."""
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        
        self.comment = lines[0].strip()
        self.scaling_factor = float(lines[1].strip())
        
        # Lattice vectors
        self.lattice = np.array([
            [float(x) for x in lines[2].split()],
            [float(x) for x in lines[3].split()],
            [float(x) for x in lines[4].split()]
        ])
        
        # Element names and counts
        self.elements = lines[5].split()
        self.element_counts = [int(x) for x in lines[6].split()]
        self.total_atoms = sum(self.element_counts)
        
        # Selective dynamics (optional)
        self.selective_dynamics = False
        coord_line_idx = 7
        if lines[7].strip()[0] in ['S', 's']:
            self.selective_dynamics = True
            self.sd_flags = []
            coord_line_idx = 8
        
        # Coordinate type (Direct/Cartesian)
        self.coord_type = lines[coord_line_idx].strip()
        self.is_direct = self.coord_type[0] in ['D', 'd']
        
        # Read atomic positions
        self.positions = []
        start_idx = coord_line_idx + 1
        
        for i in range(self.total_atoms):
            line_parts = lines[start_idx + i].split()
            pos = [float(x) for x in line_parts[:3]]
            self.positions.append(pos)
            
            if self.selective_dynamics:
                sd = line_parts[3:6] if len(line_parts) >= 6 else ['T', 'T', 'T']
                self.sd_flags.append(sd)
        
        self.positions = np.array(self.positions)
        
        # Store any additional lines (velocities, etc.)
        self.additional_lines = lines[start_idx + self.total_atoms:]
    
    def get_atom_indices(self, atom_selector):
        """
        Get indices of atoms matching the selector.
        
        Args:
            atom_selector: Can be:
                - Integer: atom ID (1-based index in POSCAR list)
                - String: element symbol (e.g., "Fe", "O") for all atoms of that type
        
        Returns:
            List of atom indices (0-based)
        """
        indices = []
        
        try:
            # Case 1: Integer input - atom ID (1-based)
            atom_id = int(atom_selector)
            if atom_id < 1 or atom_id > self.total_atoms:
                raise ValueError(f"Atom ID {atom_id} out of range (1-{self.total_atoms})")
            
            # Convert to 0-based index
            indices = [atom_id - 1]
            
        except ValueError:
            # Case 2: String input - element symbol
            element = atom_selector.strip()
            
            if element not in self.elements:
                raise ValueError(f"Element '{element}' not found in POSCAR. Available: {self.elements}")
            
            # Find all atoms of this element type
            atom_idx = 0
            for elem, count in zip(self.elements, self.element_counts):
                if elem == element:
                    # Add all atoms of this element type
                    indices.extend(range(atom_idx, atom_idx + count))
                atom_idx += count
        
        return indices
    
    def parse_atom_selectors(self, atom_string):
        """
        Parse comma-separated atom selectors.
        
        Args:
            atom_string: Comma-separated list of atom IDs and/or element symbols
                        e.g., "1,3,5", "Fe", "1,Fe,O", "2,4,6,Fe"
        
        Returns:
            List of unique atom indices (0-based), sorted
        """
        selectors = [s.strip() for s in atom_string.split(',')]
        all_indices = set()  # Use set to avoid duplicates
        
        for selector in selectors:
            if not selector:
                continue
            indices = self.get_atom_indices(selector)
            all_indices.update(indices)
        
        return sorted(list(all_indices))
    
    def displace_atoms(self, atom_string, direction, amount):
        """
        Displace selected atoms.
        
        Args:
            atom_string: Comma-separated atom IDs and/or element symbols
            direction: 'x', 'y', 'z' (cartesian) or 'a', 'b', 'c' (crystallographic)
            amount: Displacement in Angstrom
        """
        indices = self.parse_atom_selectors(atom_string)
        
        if not indices:
            raise ValueError(f"No atoms found matching selector: {atom_string}")
        
        # Print information about what's being displaced
        print(f"Displacing {len(indices)} atom(s): {atom_string}")
        print(f"Atom indices (1-based): {[i+1 for i in indices]}")
        
        # Determine displacement vector
        if direction.lower() in ['x', 'y', 'z']:
            # Cartesian displacement
            cart_disp = np.zeros(3)
            axis_map = {'x': 0, 'y': 1, 'z': 2}
            cart_disp[axis_map[direction.lower()]] = amount / self.scaling_factor
            
            # Convert to direct coordinates if needed
            if self.is_direct:
                # Transform cartesian displacement to fractional
                disp = np.linalg.solve(self.lattice.T, cart_disp)
            else:
                disp = cart_disp
                
        elif direction.lower() in ['a', 'b', 'c']:
            # Crystallographic displacement
            axis_map = {'a': 0, 'b': 1, 'c': 2}
            axis_idx = axis_map[direction.lower()]
            
            # Get the length of the lattice vector
            lattice_vector_length = np.linalg.norm(self.lattice[axis_idx]) * self.scaling_factor
            
            # Convert Angstrom displacement to fractional
            frac_disp = amount / lattice_vector_length
            
            if self.is_direct:
                disp = np.zeros(3)
                disp[axis_idx] = frac_disp
            else:
                # Convert fractional to cartesian
                disp = frac_disp * self.lattice[axis_idx]
        else:
            raise ValueError(f"Invalid direction: {direction}. Use x/y/z (cartesian) or a/b/c (crystallographic)")
        
        # Apply displacement
        for idx in indices:
            self.positions[idx] += disp
        
        print(f"Applied displacement: {amount} Å along {direction.upper()}")
    
    def write_poscar(self, filename):
        """Write POSCAR to file."""
        with open(filename, 'w') as f:
            # Header
            f.write(f"{self.comment}\n")
            f.write(f"  {self.scaling_factor:.16f}\n")
            
            # Lattice vectors
            for vec in self.lattice:
                f.write(f"  {vec[0]:22.16f} {vec[1]:22.16f} {vec[2]:22.16f}\n")
            
            # Elements
            f.write("  " + "  ".join(self.elements) + "\n")
            f.write("  " + "  ".join(str(c) for c in self.element_counts) + "\n")
            
            # Selective dynamics
            if self.selective_dynamics:
                f.write("Selective dynamics\n")
            
            # Coordinate type
            f.write(f"{self.coord_type}\n")
            
            # Positions
            for i, pos in enumerate(self.positions):
                line = f"  {pos[0]:20.16f} {pos[1]:20.16f} {pos[2]:20.16f}"
                if self.selective_dynamics:
                    line += f"  {self.sd_flags[i][0]}  {self.sd_flags[i][1]}  {self.sd_flags[i][2]}"
                f.write(line + "\n")
            
            # Additional lines (velocities, etc.)
            for line in self.additional_lines:
                f.write(line)


def main():
    parser = argparse.ArgumentParser(
        description='Displace atoms in a VASP POSCAR file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s POSCAR --atom 3 --direction x --amount 0.5
      Displace atom ID 3 (third atom in POSCAR) by 0.5 Å along cartesian x
  
  %(prog)s POSCAR --atom Fe --direction z --amount -0.3
      Displace all Fe atoms by -0.3 Å along cartesian z
  
  %(prog)s POSCAR --atom 1,3,5 --direction a --amount 0.2
      Displace atoms 1, 3, and 5 by 0.2 Å along crystallographic a
  
  %(prog)s POSCAR --atom 2,4,Fe --direction y --amount 0.1
      Displace atoms 2, 4, and all Fe atoms by 0.1 Å along cartesian y
  
  %(prog)s POSCAR --atom Fe,O --direction b --amount -0.15
      Displace all Fe and O atoms by -0.15 Å along crystallographic b
        """
    )
    
    parser.add_argument('input_file', help='Input POSCAR file')
    parser.add_argument('--atom', required=True, 
                       help='Comma-separated atom selectors: atom IDs (1-based integers) and/or element symbols (e.g., "1,3,5" or "Fe" or "2,4,6,Fe")')
    parser.add_argument('--direction', required=True, choices=['x', 'y', 'z', 'a', 'b', 'c', 'X', 'Y', 'Z', 'A', 'B', 'C'],
                       help='Displacement direction: x/y/z (cartesian) or a/b/c (crystallographic)')
    parser.add_argument('--amount', type=float, required=True,
                       help='Displacement amount in Angstrom')
    parser.add_argument('--output', default='POSCAR_disp.vasp',
                       help='Output filename (default: POSCAR_disp.vasp)')
    
    args = parser.parse_args()
    
    try:
        # Read POSCAR
        print(f"Reading {args.input_file}...")
        poscar = POSCAR(args.input_file)
        print(f"Found {poscar.total_atoms} atoms: {dict(zip(poscar.elements, poscar.element_counts))}")
        
        # Displace atoms
        poscar.displace_atoms(args.atom, args.direction, args.amount)
        
        # Write output
        poscar.write_poscar(args.output)
        print(f"Written displaced structure to {args.output}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()