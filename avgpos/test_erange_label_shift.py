#!/usr/bin/env python3
"""
Test script to validate that --erange shifts axis labels to start from 0.
This is the new behavior implemented as per the feature request.
"""

import sys
import os
import subprocess
import tempfile
import re


def test_erange_shifts_labels_to_zero():
    """Test that erange shifts axis labels to start from 0."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with --erange=10,20,11,25
        # This should create a plot where:
        # - Data range is e=[10,20], f=[11,25]
        # - Displayed labels are e=[0,10], f=[0,14]
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--erange=10,20,11,25'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that plot script was generated
        if not os.path.exists(plot_script):
            print(f"FAILED: Plot script {plot_script} was not generated")
            return False
        
        # Read plot script and verify it has the shift variables
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for the shift calculation
        if 'e_display_shift = e_min' not in script_content:
            print("FAILED: e_display_shift calculation not found")
            return False
        
        if 'f_display_shift = f_min' not in script_content:
            print("FAILED: f_display_shift calculation not found")
            return False
        
        # Check that display coordinates are used
        if 'e_display = e - e_display_shift' not in script_content:
            print("FAILED: e_display calculation not found")
            return False
        
        if 'f_display = f - f_display_shift' not in script_content:
            print("FAILED: f_display calculation not found")
            return False
        
        # Check that display coordinates are used in plotting
        if 'e_mesh_display' not in script_content or 'f_mesh_display' not in script_content:
            print("FAILED: Display mesh coordinates not found")
            return False
        
        if 'pcolormesh(e_mesh_display, f_mesh_display' not in script_content:
            print("FAILED: pcolormesh not using display coordinates")
            return False
        
        if 'scatter(e_display, f_display' not in script_content:
            print("FAILED: scatter not using display coordinates")
            return False
        
        print("PASSED: erange shifts axis labels to start from 0")
        return True
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_no_shift_without_erange():
    """Test that without erange, no shift is applied (backward compatibility)."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
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
        
        # Check that plot script was generated
        if not os.path.exists(plot_script):
            print(f"FAILED: Plot script {plot_script} was not generated")
            return False
        
        # Read plot script and verify shift is 0
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check that shift is set to 0 when erange is None
        if 'e_display_shift = 0' not in script_content:
            print("FAILED: e_display_shift should be 0 without erange")
            return False
        
        if 'f_display_shift = 0' not in script_content:
            print("FAILED: f_display_shift should be 0 without erange")
            return False
        
        print("PASSED: No shift applied without erange (backward compatible)")
        return True
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_label_coordinates_shifted():
    """Test that labels are positioned using shifted coordinates."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with --erange and --labels
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--labels',
            '--erange=10,20,11,25'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that plot script was generated
        if not os.path.exists(plot_script):
            print(f"FAILED: Plot script {plot_script} was not generated")
            return False
        
        # Read plot script and verify labels use display coordinates
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check that annotate uses e_display and f_display
        if 'annotate(labels[i], (e_display[i], f_display[i])' not in script_content:
            print("FAILED: annotate not using display coordinates")
            return False
        
        print("PASSED: Labels positioned using shifted coordinates")
        return True
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def main():
    """Run all tests."""
    print("Testing erange label shifting functionality")
    print("=" * 60)
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    tests = [
        ("erange shifts labels to start from 0", test_erange_shifts_labels_to_zero),
        ("No shift without erange (backward compatible)", test_no_shift_without_erange),
        ("Labels positioned with shifted coordinates", test_label_coordinates_shifted)
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
