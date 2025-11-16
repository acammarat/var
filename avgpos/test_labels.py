#!/usr/bin/env python3
"""
Simple test script to validate the --labels option functionality.
Tests different label formats: type, id, and both.
"""

import sys
import os
import subprocess
import tempfile

def run_test(label_format, expected_labels):
    """Run avgpos with a specific label format and verify output."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as f:
        output_file = f.name
    
    try:
        # Construct command
        cmd = [
            sys.executable, 
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', output_file,
            '--labels'
        ]
        
        if label_format:
            cmd.append(label_format)
        
        # Run command
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Read and verify output
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        # Extract labels from data lines (skip comments)
        actual_labels = []
        for line in lines:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                if len(parts) == 4:  # e, f, g, label
                    actual_labels.append(parts[3])
        
        # Verify labels
        if actual_labels == expected_labels:
            print(f"PASSED: --labels {label_format if label_format else '(no arg)'}")
            return True
        else:
            print(f"FAILED: --labels {label_format if label_format else '(no arg)'}")
            print(f"  Expected: {expected_labels}")
            print(f"  Got:      {actual_labels}")
            return False
    
    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)

def main():
    """Run all tests."""
    print("Testing --labels option functionality")
    print("=" * 60)
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    tests = [
        ('type', ['Se', 'Se', 'Se', 'Se']),
        ('id', ['2', '3', '4', '5']),
        ('both', ['Se2', 'Se3', 'Se4', 'Se5']),
        (None, ['Se2', 'Se3', 'Se4', 'Se5'])  # Default to 'both'
    ]
    
    passed = 0
    failed = 0
    
    for label_format, expected in tests:
        if run_test(label_format, expected):
            passed += 1
        else:
            failed += 1
    
    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    
    return 0 if failed == 0 else 1

if __name__ == '__main__':
    sys.exit(main())
