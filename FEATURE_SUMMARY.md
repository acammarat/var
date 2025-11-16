# Label at Projection Feature

## Overview
This feature adds the `--label-at-projection` option to position atom labels at the exact projection coordinates on the plane, replacing the circles.

## Visual Comparison

### Before (Default behavior)
- Labels are offset from circles by (5, 5) points
- Circles are visible at projection coordinates
- Command: `avgpos.py POSCAR -s Se -d z -o out.dat --plot --labels`
- See: `comparison_default_labels.png`

### After (With --label-at-projection)
- Labels are centered at exact projection coordinates
- Circles are hidden (labels replace them)
- Command: `avgpos.py POSCAR -s Se -d z -o out.dat --plot --labels --label-at-projection`
- See: `comparison_label_at_projection.png`

## Code Changes

### Label Positioning
**Before:**
```python
ax.annotate(labels[i], (e[i], f[i]), 
            xytext=(5, 5), textcoords='offset points',
            fontsize=10, fontweight='bold', color='black',
            bbox=dict(...))
```

**After (with --label-at-projection):**
```python
ax.annotate(labels[i], (e[i], f[i]), 
            ha='center', va='center',
            fontsize=10, fontweight='bold', color='black',
            bbox=dict(...))
```

### Circle Visibility
- Circles are hidden when both `--label-at-projection` AND `--labels` are used
- Otherwise circles follow existing behavior

## Testing
All tests pass (20/20):
- ✅ test_labels.py: 4/4 passed
- ✅ test_label_no_box.py: 3/3 passed
- ✅ test_label_at_projection.py: 4/4 passed (NEW)
- ✅ test_flip_g.py: 2/2 passed
- ✅ test_vrange.py: 7/7 passed

## Security
- ✅ CodeQL scan: 0 vulnerabilities found

## Use Cases
This feature is useful when:
- Atom identities are more important than circles
- You want a cleaner, less cluttered plot
- Labels should be the primary visual element at each projection point
