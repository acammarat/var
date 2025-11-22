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

def read_ascii_matrix(filename):
    """Read and parse the ASCII matrix file."""
    import numpy as np
    
    header_dict = {}
    data_lines = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                # Parse header
                if 'Grid resolution:' in line:
                    # Extract XRes x YRes
                    parts = line.split(':')[1].strip().split('x')
                    header_dict['XRes'] = parts[0].strip()
                    header_dict['YRes'] = parts[1].strip()
                elif 'X range' in line:
                    # Extract range info
                    parts = line.split('(width:')[1].split('Å')[0].strip()
                    header_dict['XReal'] = parts
                elif 'Y range' in line:
                    # Extract range info
                    parts = line.split('(width:')[1].split('Å')[0].strip()
                    header_dict['YReal'] = parts
            elif line:  # Non-empty, non-comment line
                # Parse data line
                data_lines.append(line)
    
    # Parse all data values
    data = []
    for line in data_lines:
        values = [float(v) for v in line.split()]
        data.extend(values)
    
    return header_dict, data


def test_gwyddion_basic():
    """Test basic Gwyddion ASCII matrix file generation."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as dat_f:
        dat_file = dat_f.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as gsf_f:
        gsf_file = gsf_f.name
    
    # Get plot script filename
    plot_script = dat_file.replace('.dat', '_plot.py')
    
    try:
        # Run avgpos with --gwyddion and --plot options
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', dat_file,
            '--plot',
            '--gwyddion', gsf_file
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Verify message about ASCII matrix file
        if 'ASCII matrix file will be written to' not in result.stdout:
            print(f"FAILED: Expected message about ASCII matrix not found in output")
            return False
        
        # Run the generated plot script to create the ASCII file
        if not os.path.exists(plot_script):
            print(f"FAILED: Plot script was not created: {plot_script}")
            return False
        
        plot_result = subprocess.run([sys.executable, plot_script], capture_output=True, text=True)
        
        if plot_result.returncode != 0:
            print(f"FAILED: Plot script failed with return code {plot_result.returncode}")
            print(f"STDERR: {plot_result.stderr}")
            return False
        
        # Check that ASCII file was created by the plot script
        if not os.path.exists(gsf_file):
            print(f"FAILED: ASCII matrix file was not created: {gsf_file}")
            return False
        
        # Read and validate ASCII matrix file
        try:
            header, data = read_ascii_matrix(gsf_file)
        except Exception as e:
            print(f"FAILED: Error reading ASCII matrix file: {e}")
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
        
        print(f"PASSED: Basic Gwyddion ASCII matrix file generation")
        print(f"  XRes: {xres}, YRes: {yres}")
        print(f"  Data points: {len(data)}")
        print(f"  File size: {os.path.getsize(gsf_file)} bytes")
        return True
    
    finally:
        # Clean up
        for f in [dat_file, gsf_file, plot_script]:
            if os.path.exists(f):
                os.remove(f)
        # Also clean up the PNG file
        png_file = dat_file.replace('.dat', '_heatmap.png')
        if os.path.exists(png_file):
            os.remove(png_file)


def test_gwyddion_without_output():
    """Test that --gwyddion requires both -o and --plot options."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as gsf_f:
        gsf_file = gsf_f.name
    
    try:
        # Run avgpos with --gwyddion but without -o or --plot
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
            print(f"FAILED: Expected warning about missing -o and --plot options")
            return False
        
        print(f"PASSED: Warning shown when -o is not specified")
        return True
    
    finally:
        # Clean up
        if os.path.exists(gsf_file):
            os.remove(gsf_file)


def test_gwyddion_with_labels():
    """Test Gwyddion export with labels (labels should not affect ASCII matrix output)."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as dat_f:
        dat_file = dat_f.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as gsf_f:
        gsf_file = gsf_f.name
    
    # Get plot script filename
    plot_script = dat_file.replace('.dat', '_plot.py')
    
    try:
        # Run avgpos with --gwyddion, --plot and --labels
        cmd = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', dat_file,
            '--plot',
            '--gwyddion', gsf_file,
            '--labels', 'both'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"FAILED: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Run the generated plot script to create the ASCII file
        if not os.path.exists(plot_script):
            print(f"FAILED: Plot script was not created: {plot_script}")
            return False
        
        plot_result = subprocess.run([sys.executable, plot_script], capture_output=True, text=True)
        
        if plot_result.returncode != 0:
            print(f"FAILED: Plot script failed with return code {plot_result.returncode}")
            print(f"STDERR: {plot_result.stderr}")
            return False
        
        # Check that ASCII file was created by the plot script
        if not os.path.exists(gsf_file):
            print(f"FAILED: ASCII matrix file was not created: {gsf_file}")
            return False
        
        # Read ASCII matrix file to ensure it's valid
        try:
            header, data = read_ascii_matrix(gsf_file)
            xres = int(header['XRes'])
            yres = int(header['YRes'])
        except Exception as e:
            print(f"FAILED: Error reading ASCII matrix file: {e}")
            return False
        
        print(f"PASSED: Gwyddion export works with --labels option")
        print(f"  Grid size: {xres}x{yres}")
        return True
    
    finally:
        # Clean up
        for f in [dat_file, gsf_file, plot_script]:
            if os.path.exists(f):
                os.remove(f)
        # Also clean up the PNG file
        png_file = dat_file.replace('.dat', '_heatmap.png')
        if os.path.exists(png_file):
            os.remove(png_file)


def test_gwyddion_not_affected_by_erange():
    """Test that Gwyddion export always uses full data range, not affected by --erange."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.dat', delete=False) as dat_f:
        dat_file = dat_f.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='_full.txt', delete=False) as gsf_full_f:
        gsf_full_file = gsf_full_f.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='_erange.txt', delete=False) as gsf_erange_f:
        gsf_erange_file = gsf_erange_f.name
    
    # Get plot script filenames
    plot_script_full = dat_file.replace('.dat', '_full_plot.py')
    plot_script_erange = dat_file.replace('.dat', '_erange_plot.py')
    
    try:
        # First, generate Gwyddion file WITHOUT erange
        cmd_full = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', dat_file.replace('.dat', '_full.dat'),
            '--plot',
            '--gwyddion', gsf_full_file
        ]
        
        result_full = subprocess.run(cmd_full, capture_output=True, text=True)
        if result_full.returncode != 0:
            print(f"FAILED: Full command failed")
            return False
        
        if not os.path.exists(plot_script_full):
            print(f"FAILED: Full plot script not created")
            return False
        
        plot_result_full = subprocess.run([sys.executable, plot_script_full], capture_output=True, text=True)
        if plot_result_full.returncode != 0:
            print(f"FAILED: Full plot script failed")
            return False
        
        # Now generate Gwyddion file WITH erange (should be the same)
        cmd_erange = [
            sys.executable,
            'avgpos.py',
            'example/POSCAR',
            '-s', 'Se',
            '-d', 'z',
            '-o', dat_file.replace('.dat', '_erange.dat'),
            '--plot',
            '--gwyddion', gsf_erange_file,
            '--erange=0,1,0,1'  # Restrictive erange
        ]
        
        result_erange = subprocess.run(cmd_erange, capture_output=True, text=True)
        if result_erange.returncode != 0:
            print(f"FAILED: Erange command failed")
            return False
        
        plot_script_erange = dat_file.replace('.dat', '_erange_plot.py')
        if not os.path.exists(plot_script_erange):
            print(f"FAILED: Erange plot script not created")
            return False
        
        plot_result_erange = subprocess.run([sys.executable, plot_script_erange], capture_output=True, text=True)
        if plot_result_erange.returncode != 0:
            print(f"FAILED: Erange plot script failed")
            return False
        
        # Read both ASCII files
        try:
            header_full, data_full = read_ascii_matrix(gsf_full_file)
            header_erange, data_erange = read_ascii_matrix(gsf_erange_file)
        except Exception as e:
            print(f"FAILED: Error reading ASCII matrix files: {e}")
            return False
        
        # Verify both have the same dimensions (full data)
        if header_full['XRes'] != header_erange['XRes'] or header_full['YRes'] != header_erange['YRes']:
            print(f"FAILED: Grid dimensions differ")
            print(f"  Full: {header_full['XRes']}x{header_full['YRes']}")
            print(f"  Erange: {header_erange['XRes']}x{header_erange['YRes']}")
            return False
        
        # Verify both have the same data range (full data)
        if header_full['XReal'] != header_erange['XReal'] or header_full['YReal'] != header_erange['YReal']:
            print(f"FAILED: Data ranges differ")
            print(f"  Full: X={header_full['XReal']}, Y={header_full['YReal']}")
            print(f"  Erange: X={header_erange['XReal']}, Y={header_erange['YReal']}")
            return False
        
        # Verify data is identical
        if len(data_full) != len(data_erange):
            print(f"FAILED: Data size differs")
            return False
        
        # Check that data values are the same (within floating point tolerance)
        import numpy as np
        if not np.allclose(data_full, data_erange, rtol=1e-9):
            print(f"FAILED: Data values differ")
            return False
        
        print(f"PASSED: Gwyddion export not affected by --erange")
        print(f"  Both files have identical: {header_full['XRes']}x{header_full['YRes']} grid, {len(data_full)} points")
        return True
    
    finally:
        # Clean up
        cleanup_files = [
            dat_file, gsf_full_file, gsf_erange_file, 
            plot_script_full, plot_script_erange,
            dat_file.replace('.dat', '_full.dat'),
            dat_file.replace('.dat', '_erange.dat'),
            dat_file.replace('.dat', '_full_heatmap.png'),
            dat_file.replace('.dat', '_erange_heatmap.png'),
        ]
        for f in cleanup_files:
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
        ('GSF not affected by erange', test_gwyddion_not_affected_by_erange),
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
