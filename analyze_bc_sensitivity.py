#!/usr/bin/env python3
"""
Boundary Condition Sensitivity Analysis
Compares uniform vs parabolic inlet velocity profiles.
"""

import numpy as np
from pathlib import Path

Z_EXPANSION = 0.122685
RHO = 1056.0

def read_openfoam_sample(case_dir, set_name):
    """Read sample data from postProcessing."""
    post = case_dir / 'postProcessing' / 'sampleDict'
    if not post.exists():
        return None, None, None
    
    time_dirs = sorted([d for d in post.iterdir() if d.is_dir()],
                      key=lambda x: float(x.name) if x.name.replace('.','').isdigit() else 0)
    if not time_dirs:
        return None, None, None
    
    sample_file = time_dirs[-1] / f"{set_name}_p_U.xy"
    if not sample_file.exists():
        return None, None, None
    
    data = np.loadtxt(sample_file)
    return data[:, 0], data[:, 1], data[:, 2:5]


def get_centerline_velocity_at_z(case_dir, z_exp):
    """Get centerline axial velocity at specified z (experimental coords)."""
    pos, p, U = read_openfoam_sample(case_dir, 'centerline')
    if pos is None:
        return None
    z_sim = pos - Z_EXPANSION
    idx = np.argmin(np.abs(z_sim - z_exp))
    return U[idx, 2]


def get_pressure_drop(case_dir):
    """Get pressure drop from inlet to outlet."""
    pos, p_kin, U = read_openfoam_sample(case_dir, 'centerline')
    if pos is None:
        return None
    p = p_kin * RHO
    inlet_idx = np.argmin(np.abs(pos - 0.0))
    outlet_idx = np.argmin(np.abs(pos - 0.242685))
    return p[inlet_idx] - p[outlet_idx]


def get_jet_width(case_dir, set_name):
    """Get jet width at specified location."""
    pos, p, U = read_openfoam_sample(case_dir, set_name)
    if pos is None:
        return None
    
    center_idx = np.argmin(np.abs(pos))
    u_center = U[center_idx, 2]
    if u_center <= 0:
        return None
    
    u_half = u_center / 2
    r_pos = pos[pos >= 0]
    u_pos = U[pos >= 0, 2]
    
    for i in range(len(u_pos) - 1):
        if u_pos[i] >= u_half and u_pos[i+1] < u_half:
            t = (u_half - u_pos[i]) / (u_pos[i+1] - u_pos[i])
            r_half = r_pos[i] + t * (r_pos[i+1] - r_pos[i])
            return 2 * r_half * 1000
    return None


def main():
    print("=" * 70)
    print("BC Sensitivity Analysis - Uniform vs Parabolic Inlet Profile")
    print("=" * 70)
    
    bc_dir = Path(__file__).parent / 'BC_sensitivity'
    uniform_dir = bc_dir / 'uniform_inlet'
    parabolic_dir = bc_dir / 'parabolic_inlet'
    
    # Quantities to compare
    results = []
    
    # Centerline velocity at key locations
    for z_exp, name in [(-0.020, 'Throat (z=-20mm)'),
                        (0.0, 'Expansion (z=0)'),
                        (0.032, 'Downstream (z=+32mm)'),
                        (0.080, 'Far field (z=+80mm)')]:
        u_uniform = get_centerline_velocity_at_z(uniform_dir, z_exp)
        u_parabolic = get_centerline_velocity_at_z(parabolic_dir, z_exp)
        if u_uniform and u_parabolic:
            diff = abs(u_uniform - u_parabolic) / u_parabolic * 100
            results.append((f'U_z at {name}', u_uniform, u_parabolic, diff, 'm/s'))
    
    # Pressure drop
    dp_uniform = get_pressure_drop(uniform_dir)
    dp_parabolic = get_pressure_drop(parabolic_dir)
    if dp_uniform and dp_parabolic:
        diff = abs(dp_uniform - dp_parabolic) / dp_parabolic * 100
        results.append(('Pressure Drop', dp_uniform, dp_parabolic, diff, 'Pa'))
    
    # Jet width
    for set_name, name in [('radial_z_000', 'z=0'), 
                           ('radial_z_plus032', 'z=+32mm')]:
        jw_uniform = get_jet_width(uniform_dir, set_name)
        jw_parabolic = get_jet_width(parabolic_dir, set_name)
        if jw_uniform and jw_parabolic:
            diff = abs(jw_uniform - jw_parabolic) / jw_parabolic * 100
            results.append((f'Jet Width at {name}', jw_uniform, jw_parabolic, diff, 'mm'))
    
    # Print results
    print(f"\n{'Quantity':<30} {'Uniform':>12} {'Parabolic':>12} {'Diff %':>10}")
    print("-" * 70)
    
    for name, val_u, val_p, diff, unit in results:
        print(f"{name:<30} {val_u:>11.4f} {val_p:>12.4f} {diff:>9.2f}%")
    
    print("\n" + "=" * 70)
    print("Interpretation:")
    print("=" * 70)
    
    avg_diff = np.mean([r[3] for r in results])
    max_diff = max([r[3] for r in results])
    
    print(f"\nAverage sensitivity: {avg_diff:.2f}%")
    print(f"Maximum sensitivity: {max_diff:.2f}%")
    
    if max_diff < 2:
        print("\n✅ LOW SENSITIVITY: Inlet profile has minimal impact on results.")
        print("   Either BC is acceptable for this geometry.")
    elif max_diff < 5:
        print("\n⚠️  MODERATE SENSITIVITY: Inlet profile affects some quantities.")
        print("   Parabolic (fully-developed) is more physically accurate.")
    else:
        print("\n❌ HIGH SENSITIVITY: Results are strongly dependent on inlet BC.")
        print("   Must use correct (parabolic) profile for validation.")


if __name__ == "__main__":
    main()
