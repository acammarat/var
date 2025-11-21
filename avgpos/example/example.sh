#!/bin/bash
# Example usage of avgpos tool

echo "====== Example 1: Average position of Se atoms along z-axis ======"
python3 ../avgpos.py POSCAR -s Se -d z

echo ""
echo "====== Example 2: Average position of atoms 2,3,4 along c lattice vector ======"
python3 ../avgpos.py POSCAR -i 2,3,4 -d c

echo ""
echo "====== Example 3: Average position of W and Mo atoms along [1,1,0] direction ======"
python3 ../avgpos.py POSCAR -s W,Mo -d "[1,1,0]"

echo ""
echo "====== Example 4: Average position of all Se atoms along x-axis ======"
python3 ../avgpos.py POSCAR -s Se -d x

echo ""
echo "====== Example 5: Advanced visualization with labels, replication, erange, and Gwyddion export ======"
python3 ../avgpos.py MPOSCAR_135m.vasp -i 5,29,1,21,17,41 -d b --plot --labels type --label-no-box --replicate 6,6 --vrange=-0.3,0.3 --label-at-projection -o mod135m_s1.dat --gwyddion mod135m_s1_matrix.dat --erange=10,30,5,10
echo ""
echo "To generate the visualization, run: python3 mod135m_s1_plot.py"
echo "This will create:"
echo "  - mod135m_s1_heatmap.png (visualization)"
echo "  - mod135m_s1_matrix.dat (Gwyddion ASCII matrix)"
