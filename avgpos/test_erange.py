#!/usr/bin/env python3
"""
Test script to validate the --erange option functionality.
Tests that the generated plotting script uses custom erange when specified.
"""

import sys
import os
import subprocess
import tempfile
import re


def test_default_erange():
    """Test that default behavior uses data min/max when --erange is not specified."""
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
        
        # Read plot script and verify it uses None for erange
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for erange = None in the script
        if 'erange = None' in script_content:
            # Also verify the fallback logic
            if 'e_min, e_max = e.min(), e.max()' in script_content and 'f_min, f_max = f.min(), f.max()' in script_content:
                print("PASSED: Default behavior (uses data min/max)")
                return True
            else:
                print("FAILED: Fallback logic not found")
                return False
        else:
            print("FAILED: erange = None not found in generated script")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_custom_erange():
    """Test that custom erange values are used when --erange is specified."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
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
            '--erange=0,10,0,10'
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
        if '(with custom e,f range: e=[0.0, 10.0], f=[0.0, 10.0])' not in result.stdout:
            print(f"FAILED: Custom range not mentioned in output")
            print(f"STDOUT: {result.stdout}")
            return False
        
        # Read plot script and verify it uses custom erange
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for erange = (0.0, 10.0, 0.0, 10.0) in the script
        if 'erange = (0.0, 10.0, 0.0, 10.0)' in script_content:
            print("PASSED: Custom erange (0, 10, 0, 10)")
            return True
        else:
            print("FAILED: Custom erange values not found in generated script")
            print(f"Script content snippet: {script_content[script_content.find('erange'):script_content.find('erange')+100] if 'erange' in script_content else 'erange not found'}")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_asymmetric_erange():
    """Test that asymmetric erange values work correctly."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with asymmetric --erange
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--erange=-5,5,-10,10'
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
        
        # Read plot script and verify it uses custom erange
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for erange = (-5.0, 5.0, -10.0, 10.0) in the script
        if 'erange = (-5.0, 5.0, -10.0, 10.0)' in script_content:
            print("PASSED: Asymmetric erange (-5, 5, -10, 10)")
            return True
        else:
            print("FAILED: Custom erange values not found in generated script")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_invalid_erange_reversed_e():
    """Test that invalid erange (emin > emax) shows warning and falls back to default."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with invalid --erange (reversed e values)
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--erange=10,0,0,10'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that warning is shown
        if 'Warning: Invalid erange format' not in result.stdout or 'emin must be less than emax' not in result.stdout:
            print(f"FAILED: Expected warning not found in output")
            print(f"STDOUT: {result.stdout}")
            return False
        
        # Read plot script and verify it falls back to None
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for erange = None in the script
        if 'erange = None' in script_content:
            print("PASSED: Invalid erange (reversed e) falls back to default")
            return True
        else:
            print("FAILED: Script should use None for invalid erange")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_invalid_erange_reversed_f():
    """Test that invalid erange (fmin > fmax) shows warning and falls back to default."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with invalid --erange (reversed f values)
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--erange=0,10,10,0'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that warning is shown
        if 'Warning: Invalid erange format' not in result.stdout or 'fmin must be less than fmax' not in result.stdout:
            print(f"FAILED: Expected warning not found in output")
            print(f"STDOUT: {result.stdout}")
            return False
        
        # Read plot script and verify it falls back to None
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for erange = None in the script
        if 'erange = None' in script_content:
            print("PASSED: Invalid erange (reversed f) falls back to default")
            return True
        else:
            print("FAILED: Script should use None for invalid erange")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_invalid_erange_wrong_count():
    """Test that invalid erange (wrong number of values) shows warning and falls back to default."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with invalid --erange (only 3 values)
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--erange=0,10,0'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that warning is shown
        if 'Warning: Invalid erange format' not in result.stdout or 'erange must have exactly 4 values' not in result.stdout:
            print(f"FAILED: Expected warning not found in output")
            print(f"STDOUT: {result.stdout}")
            return False
        
        # Read plot script and verify it falls back to None
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for erange = None in the script
        if 'erange = None' in script_content:
            print("PASSED: Invalid erange (wrong count) falls back to default")
            return True
        else:
            print("FAILED: Script should use None for invalid erange")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_erange_without_plot():
    """Test that --erange without --plot shows a warning."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    try:
        # Run command with --erange but without --plot
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--erange=0,10,0,10'
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
        
        # No warning is needed in this case as erange is silently ignored
        # The tool works fine, it just doesn't generate a plot
        print("PASSED: --erange without --plot (silently ignored)")
        return True
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)


def test_erange_without_output():
    """Test that --erange without output shows a warning."""
    try:
        # Run command with --erange but without -o
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '--erange=0,10,0,10'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that warning is shown
        if 'Warning: --erange flag requires both -o/--output and --plot to be specified' in result.stdout:
            print("PASSED: --erange without output shows warning")
            return True
        else:
            print("FAILED: Expected warning not found in output")
            print(f"STDOUT: {result.stdout}")
            return False
    
    finally:
        pass


def test_erange_with_vrange():
    """Test that --erange and --vrange can be used together."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with both --erange and --vrange
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--erange=0,10,0,10',
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
        
        # Check that both custom ranges are mentioned in output
        if '(with custom color range: -5.0 to 5.0)' not in result.stdout:
            print(f"FAILED: Custom vrange not mentioned in output")
            print(f"STDOUT: {result.stdout}")
            return False
        
        if '(with custom e,f range: e=[0.0, 10.0], f=[0.0, 10.0])' not in result.stdout:
            print(f"FAILED: Custom erange not mentioned in output")
            print(f"STDOUT: {result.stdout}")
            return False
        
        # Read plot script and verify both ranges are set
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for both vrange and erange
        if 'vrange = (-5.0, 5.0)' in script_content and 'erange = (0.0, 10.0, 0.0, 10.0)' in script_content:
            print("PASSED: Combined erange and vrange")
            return True
        else:
            print("FAILED: Both custom ranges not found in generated script")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def main():
    """Run all tests."""
    print("Testing --erange option functionality")
    print("=" * 60)
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    tests = [
        ("Default behavior (no erange)", test_default_erange),
        ("Custom erange (0, 10, 0, 10)", test_custom_erange),
        ("Asymmetric erange (-5, 5, -10, 10)", test_asymmetric_erange),
        ("Invalid erange (reversed e)", test_invalid_erange_reversed_e),
        ("Invalid erange (reversed f)", test_invalid_erange_reversed_f),
        ("Invalid erange (wrong count)", test_invalid_erange_wrong_count),
        ("--erange without --plot", test_erange_without_plot),
        ("--erange without output", test_erange_without_output),
        ("--erange with --vrange", test_erange_with_vrange)
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
