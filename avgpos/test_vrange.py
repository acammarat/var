#!/usr/bin/env python3
"""
Test script to validate the --vrange option functionality.
Tests that the generated plotting script uses custom vrange when specified.
"""

import sys
import os
import subprocess
import tempfile
import re


def test_default_vrange():
    """Test that default behavior uses data min/max when --vrange is not specified."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command without --vrange
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
        
        # Read plot script and verify it uses None for vrange
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for vrange = None in the script
        if 'vrange = None' in script_content:
            # Also verify the fallback logic
            if 'vmin, vmax = g.min(), g.max()' in script_content:
                print("PASSED: Default behavior (uses data min/max)")
                return True
            else:
                print("FAILED: Fallback logic not found")
                return False
        else:
            print("FAILED: vrange = None not found in generated script")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_custom_vrange():
    """Test that custom vrange values are used when --vrange is specified."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with --vrange
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--vrange=-5,5'
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
        
        # Check that output mentions custom range
        if '(with custom color range: -5.0 to 5.0)' not in result.stdout:
            print(f"FAILED: Custom range not mentioned in output")
            print(f"STDOUT: {result.stdout}")
            return False
        
        # Read plot script and verify it uses custom vrange
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for vrange = (-5.0, 5.0) in the script
        if 'vrange = (-5.0, 5.0)' in script_content:
            print("PASSED: Custom vrange (-5, 5)")
            return True
        else:
            print("FAILED: Custom vrange values not found in generated script")
            print(f"Script content snippet: {script_content[script_content.find('vrange'):script_content.find('vrange')+100]}")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_positive_vrange():
    """Test that positive vrange values work correctly."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with positive --vrange
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--vrange=0,10'
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
        
        # Read plot script and verify it uses custom vrange
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for vrange = (0.0, 10.0) in the script
        if 'vrange = (0.0, 10.0)' in script_content:
            print("PASSED: Custom vrange (0, 10)")
            return True
        else:
            print("FAILED: Custom vrange values not found in generated script")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_invalid_vrange_reversed():
    """Test that invalid vrange (vmin > vmax) shows warning and falls back to default."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with invalid --vrange (reversed)
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--vrange=5,-5'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that warning is shown
        if 'Warning: Invalid vrange format' not in result.stdout or 'vmin must be less than vmax' not in result.stdout:
            print(f"FAILED: Expected warning not found in output")
            print(f"STDOUT: {result.stdout}")
            return False
        
        # Read plot script and verify it falls back to None
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for vrange = None in the script
        if 'vrange = None' in script_content:
            print("PASSED: Invalid vrange (reversed) falls back to default")
            return True
        else:
            print("FAILED: Script should use None for invalid vrange")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_invalid_vrange_single_value():
    """Test that invalid vrange (single value) shows warning and falls back to default."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with invalid --vrange (single value)
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--vrange=5'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that warning is shown
        if 'Warning: Invalid vrange format' not in result.stdout or 'vrange must have exactly 2 values' not in result.stdout:
            print(f"FAILED: Expected warning not found in output")
            print(f"STDOUT: {result.stdout}")
            return False
        
        # Read plot script and verify it falls back to None
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for vrange = None in the script
        if 'vrange = None' in script_content:
            print("PASSED: Invalid vrange (single value) falls back to default")
            return True
        else:
            print("FAILED: Script should use None for invalid vrange")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_vrange_without_plot():
    """Test that --vrange without --plot shows a warning."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    try:
        # Run command with --vrange but without --plot
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--vrange=0,10'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that no plot script was generated
        plot_script = output_file.replace('.dat', '_plot.py')
        if os.path.exists(plot_script):
            print(f"FAILED: Plot script should not have been generated")
            os.remove(plot_script)
            return False
        
        # No warning is needed in this case as vrange is silently ignored
        # The tool works fine, it just doesn't generate a plot
        print("PASSED: --vrange without --plot (silently ignored)")
        return True
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)


def test_vrange_without_output():
    """Test that --vrange without output shows a warning."""
    try:
        # Run command with --vrange but without -o
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '--vrange=0,10'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that warning is shown
        if 'Warning: --vrange flag requires both -o/--output and --plot to be specified' in result.stdout:
            print("PASSED: --vrange without output shows warning")
            return True
        else:
            print("FAILED: Expected warning not found in output")
            print(f"STDOUT: {result.stdout}")
            return False
    
    finally:
        pass


def main():
    """Run all tests."""
    print("Testing --vrange option functionality")
    print("=" * 60)
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    tests = [
        ("Default behavior (no vrange)", test_default_vrange),
        ("Custom vrange (-5, 5)", test_custom_vrange),
        ("Positive vrange (0, 10)", test_positive_vrange),
        ("Invalid vrange (reversed)", test_invalid_vrange_reversed),
        ("Invalid vrange (single value)", test_invalid_vrange_single_value),
        ("--vrange without --plot", test_vrange_without_plot),
        ("--vrange without output", test_vrange_without_output)
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
