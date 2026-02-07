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
from scipy.optimize import minimize


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
    
    def get_cartesian_positions(self):
        """Get atomic positions in Cartesian coordinates."""
        if self.is_direct:
            # Convert from fractional to Cartesian
            cart_positions = np.dot(self.positions, self.lattice) * self.scaling_factor
        else:
            cart_positions = self.positions * self.scaling_factor
        return cart_positions
    
    def calculate_center_of_mass(self, atom_indices):
        """
        Calculate center of mass of selected atoms.
        
        Args:
            atom_indices: List of atom indices (0-based)
        
        Returns:
            Center of mass position in Cartesian coordinates (Angstrom)
        """
        cart_positions = self.get_cartesian_positions()
        selected_positions = cart_positions[atom_indices]
        center_of_mass = np.mean(selected_positions, axis=0)
        return center_of_mass
    
    def find_closest_atoms(self, reference_point, atom_indices, n=3):
        """
        Find the n closest atoms to a reference point.
        
        Args:
            reference_point: Reference position in Cartesian coordinates
            atom_indices: List of atom indices to search from (0-based)
            n: Number of closest atoms to find
        
        Returns:
            List of atom indices (0-based) of the n closest atoms
        """
        cart_positions = self.get_cartesian_positions()
        selected_positions = cart_positions[atom_indices]
        
        # Calculate distances
        distances = np.linalg.norm(selected_positions - reference_point, axis=1)
        
        # Get indices of n smallest distances
        closest_local_indices = np.argsort(distances)[:n]
        
        # Map back to original atom indices
        closest_atom_indices = [atom_indices[i] for i in closest_local_indices]
        
        return closest_atom_indices
    
    def fit_plane(self, atom_indices):
        """
        Fit a plane through the given atoms.
        
        Args:
            atom_indices: List of at least 3 atom indices (0-based)
        
        Returns:
            Tuple of (normal_vector, point_on_plane)
            - normal_vector: Unit normal vector of the plane
            - point_on_plane: A point on the plane (centroid of the atoms)
        """
        cart_positions = self.get_cartesian_positions()
        points = cart_positions[atom_indices]
        
        # Calculate centroid
        centroid = np.mean(points, axis=0)
        
        # Center the points
        centered_points = points - centroid
        
        # Use SVD to find the normal vector
        # The normal is the singular vector corresponding to the smallest singular value
        _, _, vh = np.linalg.svd(centered_points)
        normal = vh[-1, :]  # Last row of V^T
        
        # Normalize the normal vector
        normal = normal / np.linalg.norm(normal)
        
        return normal, centroid
    
    def find_neighbors(self, atom_indices, distance_threshold=5.0):
        """
        Find neighbors of the given atoms within a distance threshold.
        
        Args:
            atom_indices: List of atom indices (0-based)
            distance_threshold: Maximum distance to consider as neighbor (Angstrom)
        
        Returns:
            List of neighbor atom indices (0-based), excluding the input atoms
        """
        cart_positions = self.get_cartesian_positions()
        reference_positions = cart_positions[atom_indices]
        
        neighbors = set()
        for ref_pos in reference_positions:
            distances = np.linalg.norm(cart_positions - ref_pos, axis=1)
            neighbor_mask = (distances < distance_threshold) & (distances > 0.01)
            neighbor_indices = np.where(neighbor_mask)[0]
            neighbors.update(neighbor_indices.tolist())
        
        # Remove the original atoms from the neighbor set
        for idx in atom_indices:
            neighbors.discard(idx)
        
        return sorted(list(neighbors))
    
    def optimize_plane(self, initial_normal, initial_point, neighbor_indices):
        """
        Optimize plane position and orientation to minimize distance to neighbors.
        
        Args:
            initial_normal: Initial normal vector
            initial_point: Initial point on plane
            neighbor_indices: List of neighbor atom indices (0-based)
        
        Returns:
            Tuple of (optimized_normal, optimized_point)
        """
        cart_positions = self.get_cartesian_positions()
        neighbor_positions = cart_positions[neighbor_indices]
        
        def plane_distance_sum(params):
            """
            Calculate sum of squared distances from neighbors to plane.
            params: [nx, ny, nz, px, py, pz] where n is normal and p is point on plane
            """
            normal = params[:3]
            point = params[3:]
            
            # Normalize the normal vector
            normal_norm = np.linalg.norm(normal)
            if normal_norm < 1e-10:
                return 1e10
            normal = normal / normal_norm
            
            # Calculate distances from each neighbor to the plane
            distances = np.abs(np.dot(neighbor_positions - point, normal))
            
            return np.sum(distances**2)
        
        # Initial parameters
        initial_params = np.concatenate([initial_normal, initial_point])
        
        # Optimize
        result = minimize(plane_distance_sum, initial_params, method='BFGS')
        
        # Extract optimized parameters
        optimized_normal = result.x[:3]
        optimized_normal = optimized_normal / np.linalg.norm(optimized_normal)
        optimized_point = result.x[3:]
        
        return optimized_normal, optimized_point
    
    def point_to_plane_distance(self, point, plane_normal, plane_point):
        """
        Calculate distance from a point to a plane.
        
        Args:
            point: Point coordinates
            plane_normal: Plane normal vector (should be normalized)
            plane_point: A point on the plane
        
        Returns:
            Distance from point to plane
        """
        return np.abs(np.dot(point - plane_point, plane_normal))
    
    def calculate_surface_distance(self, atom_string, neighbor_threshold=5.0):
        """
        Calculate distance from center of mass of selected atoms to closest surface.
        
        Args:
            atom_string: Comma-separated atom IDs and/or element symbols
            neighbor_threshold: Distance threshold for finding neighbors (Angstrom)
        
        Returns:
            Distance from center of mass to the optimized surface (Angstrom)
        """
        # Parse atom selection
        indices = self.parse_atom_selectors(atom_string)
        
        if not indices:
            raise ValueError(f"No atoms found matching selector: {atom_string}")
        
        if len(indices) < 3:
            raise ValueError("Need at least 3 atoms to define a surface")
        
        print(f"\n{'='*60}")
        print(f"SURFACE DISTANCE CALCULATION")
        print(f"{'='*60}")
        print(f"Selected atoms: {len(indices)} atom(s)")
        print(f"Atom indices (1-based): {[i+1 for i in indices]}")
        
        # Step 1: Calculate center of mass
        center_of_mass = self.calculate_center_of_mass(indices)
        print(f"\nCenter of mass (Cartesian, Å): [{center_of_mass[0]:.4f}, {center_of_mass[1]:.4f}, {center_of_mass[2]:.4f}]")
        
        # Step 2: Find 3 closest atoms to center of mass
        closest_atoms = self.find_closest_atoms(center_of_mass, indices, n=3)
        print(f"\n3 closest atoms to center of mass (1-based): {[i+1 for i in closest_atoms]}")
        
        # Step 3: Fit initial plane through these 3 atoms
        initial_normal, initial_point = self.fit_plane(closest_atoms)
        print(f"Initial plane normal: [{initial_normal[0]:.4f}, {initial_normal[1]:.4f}, {initial_normal[2]:.4f}]")
        print(f"Initial plane point (Å): [{initial_point[0]:.4f}, {initial_point[1]:.4f}, {initial_point[2]:.4f}]")
        
        # Step 4: Find neighbors of the 3 atoms
        neighbors = self.find_neighbors(closest_atoms, distance_threshold=neighbor_threshold)
        print(f"\nFound {len(neighbors)} neighbor atoms (within {neighbor_threshold} Å)")
        if neighbors:
            print(f"Neighbor indices (1-based): {[i+1 for i in neighbors[:10]]}{' ...' if len(neighbors) > 10 else ''}")
        
        # Step 5: Optimize plane if neighbors exist
        if neighbors:
            optimized_normal, optimized_point = self.optimize_plane(initial_normal, initial_point, neighbors)
            print(f"\nOptimized plane normal: [{optimized_normal[0]:.4f}, {optimized_normal[1]:.4f}, {optimized_normal[2]:.4f}]")
            print(f"Optimized plane point (Å): [{optimized_point[0]:.4f}, {optimized_point[1]:.4f}, {optimized_point[2]:.4f}]")
        else:
            print("\nNo neighbors found, using initial plane")
            optimized_normal = initial_normal
            optimized_point = initial_point
        
        # Step 6: Calculate distance from center of mass to plane
        distance = self.point_to_plane_distance(center_of_mass, optimized_normal, optimized_point)
        
        print(f"\n{'='*60}")
        print(f"DISTANCE FROM CENTER OF MASS TO SURFACE: {distance:.6f} Å")
        print(f"{'='*60}\n")
        
        return distance
    
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
        description='Displace atoms in a VASP POSCAR file or calculate surface distance',
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
  
  %(prog)s POSCAR --atom Fe --calculate-surface-distance
      Calculate distance from center of mass of Fe atoms to closest surface
  
  %(prog)s POSCAR --atom 1,2,3,4,5 --calculate-surface-distance --neighbor-threshold 6.0
      Calculate surface distance for atoms 1-5 with 6.0 Å neighbor threshold
        """
    )
    
    parser.add_argument('input_file', help='Input POSCAR file')
    parser.add_argument('--atom', required=True, 
                       help='Comma-separated atom selectors: atom IDs (1-based integers) and/or element symbols (e.g., "1,3,5" or "Fe" or "2,4,6,Fe")')
    parser.add_argument('--direction', choices=['x', 'y', 'z', 'a', 'b', 'c', 'X', 'Y', 'Z', 'A', 'B', 'C'],
                       help='Displacement direction: x/y/z (cartesian) or a/b/c (crystallographic)')
    parser.add_argument('--amount', type=float,
                       help='Displacement amount in Angstrom')
    parser.add_argument('--output', default='POSCAR_disp.vasp',
                       help='Output filename (default: POSCAR_disp.vasp)')
    parser.add_argument('--calculate-surface-distance', action='store_true',
                       help='Calculate distance from center of mass of selected atoms to closest surface')
    parser.add_argument('--neighbor-threshold', type=float, default=5.0,
                       help='Distance threshold for finding neighbors (default: 5.0 Angstrom)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.calculate_surface_distance:
        # Surface distance calculation mode
        if args.direction or args.amount:
            print("Warning: --direction and --amount are ignored in surface distance calculation mode", file=sys.stderr)
    else:
        # Displacement mode
        if not args.direction or args.amount is None:
            parser.error("--direction and --amount are required for displacement mode")
    
    try:
        # Read POSCAR
        print(f"Reading {args.input_file}...")
        poscar = POSCAR(args.input_file)
        print(f"Found {poscar.total_atoms} atoms: {dict(zip(poscar.elements, poscar.element_counts))}")
        
        if args.calculate_surface_distance:
            # Calculate surface distance
            poscar.calculate_surface_distance(args.atom, args.neighbor_threshold)
        else:
            # Displace atoms
            poscar.displace_atoms(args.atom, args.direction, args.amount)
            
            # Write output
            poscar.write_poscar(args.output)
            print(f"Written displaced structure to {args.output}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()