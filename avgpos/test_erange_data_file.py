#!/usr/bin/env python3
"""
Test script to validate that --erange creates a filtered data file.
Tests that the generated _erange.dat file contains only points within the specified range.
"""

import sys
import os
import subprocess
import tempfile
import numpy as np


def test_erange_data_file_created():
    """Test that _erange.dat file is created when --erange is specified."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    erange_file = output_file.replace('.dat', '_erange.dat')
    
    try:
        # Run command with --erange
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--erange=0,2,0,2'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that erange data file was created
        if not os.path.exists(erange_file):
            print(f"FAILED: Erange data file {erange_file} was not created")
            return False
        
        # Check that output mentions the file
        if 'Filtered data file written to:' not in result.stdout:
            print(f"FAILED: Filtered data file message not in output")
            print(f"STDOUT: {result.stdout}")
            return False
        
        print("PASSED: Erange data file created")
        return True
    
    finally:
        # Clean up
        for fname in [output_file, plot_script, erange_file]:
            if os.path.exists(fname):
                os.remove(fname)


def test_erange_data_file_filtering():
    """Test that the _erange.dat file contains only points within the specified range."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    erange_file = output_file.replace('.dat', '_erange.dat')
    
    try:
        # Run command with restrictive --erange
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--erange=0,1,0,1'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Read the original data file
        original_data = np.loadtxt(output_file)
        
        # Read the filtered data file
        filtered_data = np.loadtxt(erange_file)
        
        # Check that filtered data is a subset of original data
        if len(filtered_data) > len(original_data):
            print(f"FAILED: Filtered data has more points than original data")
            return False
        
        # Check that all filtered points have e and f within [0, 1]
        e_coords = filtered_data[:, 0]
        f_coords = filtered_data[:, 1]
        
        if not np.all((e_coords >= 0) & (e_coords <= 1)):
            print(f"FAILED: Some e coordinates are outside [0, 1]")
            return False
        
        if not np.all((f_coords >= 0) & (f_coords <= 1)):
            print(f"FAILED: Some f coordinates are outside [0, 1]")
            return False
        
        # Check that the output mentions the correct number of points
        if '(2 points within erange)' not in result.stdout:
            print(f"FAILED: Expected 2 points in filtered data")
            print(f"STDOUT: {result.stdout}")
            return False
        
        print("PASSED: Erange data file contains correctly filtered points")
        return True
    
    finally:
        # Clean up
        for fname in [output_file, plot_script, erange_file]:
            if os.path.exists(fname):
                os.remove(fname)


def test_erange_data_file_with_labels():
    """Test that the _erange.dat file preserves labels when used with --labels."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    erange_file = output_file.replace('.dat', '_erange.dat')
    
    try:
        # Run command with --labels and --erange
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--labels',
            '--erange=0,1,0,1'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Read the filtered data file
        with open(erange_file, 'r') as f:
            lines = f.readlines()
        
        # Check that the file has labels
        data_lines = [l for l in lines if not l.startswith('#')]
        if len(data_lines) == 0:
            print(f"FAILED: No data lines in filtered file")
            return False
        
        # Check that each data line has 4 columns (e, f, g, label)
        for line in data_lines:
            parts = line.split()
            if len(parts) != 4:
                print(f"FAILED: Expected 4 columns in data line, got {len(parts)}")
                return False
            # Check that the 4th column is a label (starts with 'Se')
            if not parts[3].startswith('Se'):
                print(f"FAILED: Expected label to start with 'Se', got {parts[3]}")
                return False
        
        print("PASSED: Erange data file preserves labels")
        return True
    
    finally:
        # Clean up
        for fname in [output_file, plot_script, erange_file]:
            if os.path.exists(fname):
                os.remove(fname)


def test_no_erange_no_file():
    """Test that no _erange.dat file is created when --erange is not specified."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    erange_file = output_file.replace('.dat', '_erange.dat')
    
    try:
        # Run command without --erange
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that erange data file was NOT created
        if os.path.exists(erange_file):
            print(f"FAILED: Erange data file should not have been created")
            os.remove(erange_file)
            return False
        
        # Check that output does NOT mention filtered file
        if 'Filtered data file written to:' in result.stdout:
            print(f"FAILED: Filtered data file should not be mentioned in output")
            return False
        
        print("PASSED: No erange data file created when --erange is not specified")
        return True
    
    finally:
        # Clean up
        for fname in [output_file, plot_script]:
            if os.path.exists(fname):
                os.remove(fname)


def test_erange_header():
    """Test that the _erange.dat file has proper header information."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    erange_file = output_file.replace('.dat', '_erange.dat')
    
    try:
        # Run command with --erange
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--erange=-5,15,-10,10'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Read the filtered data file
        with open(erange_file, 'r') as f:
            content = f.read()
        
        # Check that the header mentions the erange
        if 'Data filtered for erange: e=[-5.0, 15.0], f=[-10.0, 10.0]' not in content:
            print(f"FAILED: Header does not mention correct erange")
            print(f"Content: {content[:500]}")
            return False
        
        print("PASSED: Erange data file has proper header")
        return True
    
    finally:
        # Clean up
        for fname in [output_file, plot_script, erange_file]:
            if os.path.exists(fname):
                os.remove(fname)


def main():
    """Run all tests."""
    print("Testing --erange data file output functionality")
    print("=" * 60)
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    tests = [
        ("Erange data file created", test_erange_data_file_created),
        ("Erange data file filtering", test_erange_data_file_filtering),
        ("Erange data file with labels", test_erange_data_file_with_labels),
        ("No erange file without --erange", test_no_erange_no_file),
        ("Erange data file header", test_erange_header)
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
