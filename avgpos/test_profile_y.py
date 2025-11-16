#!/usr/bin/env python3
"""
Test script to validate the --profile-y option functionality.
Tests that profile extraction along x at a specified y value works correctly.
"""

import sys
import os
import subprocess
import tempfile
import numpy as np


def test_profile_extraction():
    """Test that profile extraction creates a valid output file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    base_name = os.path.splitext(output_file)[0]
    profile_file = f"{base_name}_profile.dat"
    
    try:
        # Run command with --profile-y
        y_value = 0.5
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            f'--profile-y={y_value}'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that profile file was created
        if not os.path.exists(profile_file):
            print(f"FAILED: Profile file {profile_file} was not created")
            return False
        
        # Read profile data file
        with open(profile_file, 'r') as f:
            lines = f.readlines()
        
        # Check header
        if not lines[0].startswith('# x g'):
            print(f"FAILED: Expected '# x g' header")
            return False
        
        # Check that y value is mentioned in header
        header = ''.join(lines[:10])
        if f'y coordinate requested: {y_value:.6f}' not in header:
            print(f"FAILED: Expected y coordinate requested in header")
            return False
        if 'y coordinate used:' not in header:
            print(f"FAILED: Expected y coordinate used in header")
            return False
        
        # Load data (skip comment lines)
        data_lines = [l for l in lines if not l.startswith('#') and l.strip()]
        if len(data_lines) == 0:
            print(f"FAILED: No data lines found in profile file")
            return False
        
        # Parse data
        data = []
        for line in data_lines:
            parts = line.split()
            if len(parts) != 2:
                print(f"FAILED: Expected 2 columns (x, g), got {len(parts)}")
                return False
            data.append([float(parts[0]), float(parts[1])])
        
        data = np.array(data)
        
        # Check that we have at least one data point
        if len(data) < 1:
            print(f"FAILED: No data points ({len(data)})")
            return False
        
        # Check that x values are sorted
        if len(data) > 1:
            x_values = data[:, 0]
            if not np.all(x_values[1:] >= x_values[:-1]):
                print(f"FAILED: x values are not sorted")
                return False
        
        # Check that output mentions profile extraction
        if 'Profile data written to:' not in result.stdout:
            print(f"FAILED: Output should mention profile data file")
            return False
        
        if f'y = {y_value:.6f}' not in result.stdout:
            print(f"FAILED: Output should mention y coordinate")
            return False
        
        print(f"PASSED: Profile extraction (y={y_value})")
        return True
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(profile_file):
            os.remove(profile_file)
        # Also clean up the generated plot script
        script_file = f"{base_name}_plot.py"
        if os.path.exists(script_file):
            os.remove(script_file)


def test_profile_with_different_y_values():
    """Test that different y values that map to different actual y values produce different profiles."""
    profiles = []
    actual_y_used = []
    # Use y values that will map to different actual y values
    # Based on the example POSCAR, actual y values are 0.0 and ~0.959367
    y_values = [0.0, 0.959367]
    
    for y_val in y_values:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
            output_file = f.name
        
        base_name = os.path.splitext(output_file)[0]
        profile_file = f"{base_name}_profile.dat"
        
        try:
            cmd = [
                sys.executable,
                'avgpos.py',
                'example/POSCAR',
                '-s', 'Se',
                '-d', 'z',
                '-o', output_file,
                '--plot',
                f'--profile-y={y_val}'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"FAILED: Command failed for y={y_val}")
                return False
            
            # Read profile data
            data = np.loadtxt(profile_file)
            g_values = data[:, 1]
            profiles.append(g_values)
            
            # Read the actual y value used from the header
            with open(profile_file, 'r') as f:
                for line in f:
                    if 'y coordinate used:' in line:
                        actual_y = float(line.split(':')[1].strip())
                        actual_y_used.append(actual_y)
                        break
        
        finally:
            # Clean up
            if os.path.exists(output_file):
                os.remove(output_file)
            if os.path.exists(profile_file):
                os.remove(profile_file)
            script_file = f"{base_name}_plot.py"
            if os.path.exists(script_file):
                os.remove(script_file)
    
    # Check that profiles are different since we used different actual y values
    if not np.allclose(profiles[0], profiles[1]):
        print(f"PASSED: Different y values produce different profiles")
        return True
    else:
        print(f"FAILED: Profiles should be different but are the same")
        return False


def test_profile_without_plot():
    """Test that --profile-y without --plot shows warning."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    try:
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--profile-y=0.5'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed")
            return False
        
        # Check that warning is shown
        if '--profile-y flag requires both -o/--output and --plot' not in result.stdout:
            print(f"FAILED: Expected warning about --profile-y requiring --plot")
            return False
        
        # Check that profile file was NOT created
        base_name = os.path.splitext(output_file)[0]
        profile_file = f"{base_name}_profile.dat"
        if os.path.exists(profile_file):
            print(f"FAILED: Profile file should not be created without --plot")
            os.remove(profile_file)
            return False
        
        print(f"PASSED: Warning shown when --plot is not specified")
        return True
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)


def test_profile_without_output():
    """Test that --profile-y without -o shows warning."""
    cmd = [
        sys.executable,
        'avgpos.py',
        'example/POSCAR',
        '-s', 'Se',
        '-d', 'z',
        '--profile-y=0.5'
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"FAILED: Command failed")
        return False
    
    # Check that warning is shown
    if '--profile-y flag requires both -o/--output and --plot' not in result.stdout:
        print(f"FAILED: Expected warning about --profile-y requiring -o")
        return False
    
    print(f"PASSED: Warning shown when -o is not specified")
    return True


def main():
    """Run all tests."""
    print("Testing --profile-y option functionality")
    print("=" * 60)
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    tests = [
        ("Profile extraction", test_profile_extraction),
        ("Different y values", test_profile_with_different_y_values),
        ("Without --plot flag", test_profile_without_plot),
        ("Without -o flag", test_profile_without_output)
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
