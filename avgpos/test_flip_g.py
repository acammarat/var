#!/usr/bin/env python3
"""
Test script to validate the --flip-g option functionality.
Tests that the g values are flipped when --flip-g is used.
"""

import sys
import os
import subprocess
import tempfile
import numpy as np


def test_default_g():
    """Test that default behavior uses g = average_position - distance."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    try:
        # Run command without --flip-g
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False, None
        
        # Read data file
        data = np.loadtxt(output_file)
        g_values = data[:, 2]  # Third column is g
        
        # Read and verify header comment
        with open(output_file, 'r') as f:
            lines = f.readlines()
            header = ''.join(lines[:10])
        
        if 'average_position - distance_from_plane' not in header:
            print(f"FAILED: Expected header not found")
            print(f"Header: {header}")
            return False, None
        
        print(f"PASSED: Default behavior (g = average_position - distance_from_plane)")
        return True, g_values
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)


def test_flipped_g():
    """Test that --flip-g flips the sign of g values."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    try:
        # Run command with --flip-g
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--flip-g'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False, None
        
        # Read data file
        data = np.loadtxt(output_file)
        g_values = data[:, 2]  # Third column is g
        
        # Read and verify header comment
        with open(output_file, 'r') as f:
            lines = f.readlines()
            header = ''.join(lines[:10])
        
        if 'distance_from_plane - average_position' not in header:
            print(f"FAILED: Expected header not found in flipped mode")
            print(f"Header: {header}")
            return False, None
        
        # Verify output mentions [flipped]
        if '[flipped]' not in result.stdout:
            print(f"FAILED: Output should mention [flipped]")
            print(f"STDOUT: {result.stdout}")
            return False, None
        
        print(f"PASSED: --flip-g (g = distance_from_plane - average_position)")
        return True, g_values
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)


def test_g_values_are_opposite():
    """Test that g values with --flip-g are the negative of default g values."""
    
    # Get default g values
    success1, g_default = test_default_g()
    if not success1:
        return False
    
    # Get flipped g values
    success2, g_flipped = test_flipped_g()
    if not success2:
        return False
    
    # Check that they are opposite
    if g_default is None or g_flipped is None:
        print(f"FAILED: Could not get g values for comparison")
        return False
    
    # They should be negatives of each other
    if np.allclose(g_flipped, -g_default):
        print(f"PASSED: Flipped g values are opposite of default g values")
        return True
    else:
        print(f"FAILED: Flipped g values are not opposite of default g values")
        print(f"Default g: {g_default}")
        print(f"Flipped g: {g_flipped}")
        print(f"Expected:  {-g_default}")
        return False


def test_flip_g_with_labels():
    """Test that --flip-g works correctly with --labels."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    try:
        # Run command with --flip-g and --labels
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--flip-g',
            '--labels'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Read and verify data file has labels
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        # Check data lines have 4 columns
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
        if len(data_lines) == 0:
            print(f"FAILED: No data lines found")
            return False
        
        # Verify each data line has 4 parts (e, f, g, label)
        for line in data_lines:
            parts = line.split()
            if len(parts) != 4:
                print(f"FAILED: Expected 4 columns, got {len(parts)}")
                return False
        
        # Verify header mentions flipped
        header = ''.join(lines[:10])
        if 'distance_from_plane - average_position' not in header:
            print(f"FAILED: Expected flipped header not found with labels")
            return False
        
        print(f"PASSED: --flip-g with --labels")
        return True
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)


def main():
    """Run all tests."""
    print("Testing --flip-g option functionality")
    print("=" * 60)
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    tests = [
        ("g values are opposite", test_g_values_are_opposite),
        ("--flip-g with --labels", test_flip_g_with_labels)
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        print(f"\nTest: {test_name}")
        if test_func():
            passed += 1
        else:
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
