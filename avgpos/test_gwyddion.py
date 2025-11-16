#!/usr/bin/env python3
"""
Test script to validate the --gwyddion option functionality.
Tests the generation of Gwyddion Simple Field (GSF) format files.
"""

import sys
import os
import subprocess
import tempfile
import struct

def read_gsf_header(filename):
    """Read and parse the GSF header."""
    with open(filename, 'rb') as f:
        # Read lines until we hit the 4 null bytes
        header_lines = []
        header_bytes = b''
        
        while True:
            byte = f.read(1)
            if not byte:
                break
            header_bytes += byte
            
            # Check for header terminator (4 null bytes)
            if header_bytes.endswith(b'\x00\x00\x00\x00'):
                # Remove the terminator
                header_bytes = header_bytes[:-4]
                break
        
        # Parse header
        header_text = header_bytes.decode('utf-8')
        header_dict = {}
        for line in header_text.strip().split('\n'):
            if '=' in line:
                key, value = line.split('=', 1)
                header_dict[key.strip()] = value.strip()
        
        # Read binary data
        data_bytes = f.read()
        n_floats = len(data_bytes) // 4
        data = struct.unpack(f'{n_floats}f', data_bytes)
        
        return header_dict, data


def test_gwyddion_basic():
    """Test basic Gwyddion GSF file generation."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as dat_f:
        dat_file = dat_f.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gsf', delete=False) as gsf_f:
        gsf_file = gsf_f.name
    
    try:
        # Run avgpos with --gwyddion option
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', dat_file,
            '--gwyddion', gsf_file
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that GSF file was created
        if not os.path.exists(gsf_file):
            print(f"FAILED: GSF file was not created: {gsf_file}")
            return False
        
        # Read and validate GSF header
        try:
            header, data = read_gsf_header(gsf_file)
        except Exception as e:
            print(f"FAILED: Error reading GSF file: {e}")
            return False
        
        # Validate required header fields
        required_fields = ['XRes', 'YRes', 'XReal', 'YReal']
        for field in required_fields:
            if field not in header:
                print(f"FAILED: Missing required header field: {field}")
                return False
        
        # Validate data size matches header
        xres = int(header['XRes'])
        yres = int(header['YRes'])
        expected_size = xres * yres
        
        if len(data) != expected_size:
            print(f"FAILED: Data size mismatch. Expected {expected_size}, got {len(data)}")
            return False
        
        # Check that file contains valid float data
        if not all(isinstance(x, float) for x in data[:10]):  # Check first 10 values
            print(f"FAILED: Data does not contain valid floats")
            return False
        
        # Verify success message in output
        if 'Gwyddion GSF file written to' not in result.stdout:
            print(f"FAILED: Expected success message not found in output")
            return False
        
        print(f"PASSED: Basic Gwyddion GSF file generation")
        print(f"  XRes: {xres}, YRes: {yres}")
        print(f"  Data points: {len(data)}")
        print(f"  File size: {os.path.getsize(gsf_file)} bytes")
        return True
    
    finally:
        # Clean up
        for f in [dat_file, gsf_file]:
            if os.path.exists(f):
                os.remove(f)


def test_gwyddion_without_output():
    """Test that --gwyddion requires -o option."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gsf', delete=False) as gsf_f:
        gsf_file = gsf_f.name
    
    try:
        # Run avgpos with --gwyddion but without -o
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '--gwyddion', gsf_file
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Should succeed but show warning
        if result.returncode != 0:
            print(f"FAILED: Command failed unexpectedly")
            return False
        
        # Check for warning message
        if 'Warning' not in result.stdout or 'gwyddion' not in result.stdout.lower():
            print(f"FAILED: Expected warning about missing -o option")
            return False
        
        print(f"PASSED: Warning shown when -o is not specified")
        return True
    
    finally:
        # Clean up
        if os.path.exists(gsf_file):
            os.remove(gsf_file)


def test_gwyddion_with_labels():
    """Test Gwyddion export with labels (labels should not affect GSF output)."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as dat_f:
        dat_file = dat_f.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='.gsf', delete=False) as gsf_f:
        gsf_file = gsf_f.name
    
    try:
        # Run avgpos with --gwyddion and --labels
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', dat_file,
            '--gwyddion', gsf_file,
            '--labels', 'both'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that GSF file was created
        if not os.path.exists(gsf_file):
            print(f"FAILED: GSF file was not created: {gsf_file}")
            return False
        
        # Read GSF file to ensure it's valid
        try:
            header, data = read_gsf_header(gsf_file)
            xres = int(header['XRes'])
            yres = int(header['YRes'])
        except Exception as e:
            print(f"FAILED: Error reading GSF file: {e}")
            return False
        
        print(f"PASSED: Gwyddion export works with --labels option")
        print(f"  Grid size: {xres}x{yres}")
        return True
    
    finally:
        # Clean up
        for f in [dat_file, gsf_file]:
            if os.path.exists(f):
                os.remove(f)


def main():
    """Run all tests."""
    print("Testing --gwyddion option functionality")
    print("=" * 60)
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    tests = [
        ('Basic GSF generation', test_gwyddion_basic),
        ('Warning without -o', test_gwyddion_without_output),
        ('GSF with labels', test_gwyddion_with_labels),
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            if test_func():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"FAILED: {test_name} - Exception: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    
    return 0 if failed == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
