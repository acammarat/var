#!/usr/bin/env python3
"""
Test script to validate the --label-at-projection option functionality.
Tests that labels are positioned at projection coordinates and circles are hidden.
"""

import sys
import os
import subprocess
import tempfile
import re


def test_default_behavior():
    """Test that default behavior has circles and labels with offset."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command without --label-at-projection
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--labels',
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
        
        # Read plot script and verify it contains circles and offset labels
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check for scatter (circles)
        if 'ax.scatter(' not in script_content:
            print("FAILED: scatter not found in generated script (should be present by default)")
            return False
        
        # Check for offset labels (xytext parameter)
        if 'xytext=(5, 5)' not in script_content:
            print("FAILED: xytext offset not found in generated script (should be present by default)")
            return False
        
        print("PASSED: Default behavior (circles + offset labels)")
        return True
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_label_at_projection():
    """Test that --label-at-projection positions labels at projection coordinates and hides circles."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with --label-at-projection
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--labels',
            '--plot',
            '--label-at-projection'
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
        
        # Read plot script and verify it does NOT contain circles
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check that scatter is NOT in the script (circles should be hidden)
        if 'ax.scatter(' in script_content:
            print("FAILED: scatter found in generated script (should not be present with --label-at-projection)")
            return False
        
        # Check that labels are centered (ha='center', va='center')
        if "ha='center'" not in script_content or "va='center'" not in script_content:
            print("FAILED: Labels are not centered (should have ha='center' and va='center')")
            return False
        
        # Check that labels do NOT have xytext offset
        if 'xytext=' in script_content:
            print("FAILED: xytext found in generated script (labels should not be offset)")
            return False
        
        # Verify message in stdout
        if "at projection coordinates (no circles)" not in result.stdout:
            print("FAILED: Expected output message about labels at projection coordinates not found")
            return False
        
        print("PASSED: Labels at projection coordinates (no circles)")
        return True
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_label_at_projection_with_no_box():
    """Test that --label-at-projection works with --label-no-box."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with --label-at-projection and --label-no-box
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--labels',
            '--plot',
            '--label-at-projection',
            '--label-no-box'
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
        
        # Read plot script and verify it does NOT contain circles or bbox
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Check that scatter is NOT in the script
        if 'ax.scatter(' in script_content:
            print("FAILED: scatter found in generated script (should not be present)")
            return False
        
        # Check that bbox is NOT in the script
        if 'bbox=dict(' in script_content:
            print("FAILED: bbox found in generated script (should not be present with --label-no-box)")
            return False
        
        # Check that labels are centered
        if "ha='center'" not in script_content or "va='center'" not in script_content:
            print("FAILED: Labels are not centered")
            return False
        
        print("PASSED: Labels at projection coordinates without box")
        return True
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def test_label_at_projection_without_labels():
    """Test that --label-at-projection has no effect when --labels is not used."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    plot_script = output_file.replace('.dat', '_plot.py')
    
    try:
        # Run command with --label-at-projection but WITHOUT --labels
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--plot',
            '--label-at-projection'
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
        
        # Read plot script and verify it contains circles (default behavior)
        with open(plot_script, 'r') as f:
            script_content = f.read()
        
        # Since --labels was not used, circles should be present (default)
        if 'ax.scatter(' not in script_content:
            print("FAILED: scatter not found (should be present when --labels is not used)")
            return False
        
        # Labels should not be added at all
        if 'ax.annotate(' in script_content:
            print("FAILED: Labels found in script even though --labels was not used")
            return False
        
        print("PASSED: --label-at-projection has no effect without --labels")
        return True
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)
        if os.path.exists(plot_script):
            os.remove(plot_script)


def main():
    """Run all tests."""
    print("Testing --label-at-projection option functionality")
    print("=" * 60)
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    tests = [
        ("Default behavior", test_default_behavior),
        ("With --label-at-projection flag", test_label_at_projection),
        ("--label-at-projection with --label-no-box", test_label_at_projection_with_no_box),
        ("--label-at-projection without --labels", test_label_at_projection_without_labels)
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
