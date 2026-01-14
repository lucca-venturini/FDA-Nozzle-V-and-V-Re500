#!/usr/bin/env python3
"""
FDA Nozzle Benchmark V&V - Multi-Case Comparison
Plots all simulation cases on the same graphs for direct comparison.

Run from FDA_V&V directory:
    python compare_all_cases.py
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import re

# Set non-interactive backend
import matplotlib
matplotlib.use('Agg')

# Script location is FDA_V&V root
ROOT_DIR = Path(__file__).parent

# Constants
Z_EXPANSION = 0.122685  # Expansion plane location in simulation coordinates (m)
RHO = 1056.0            # Fluid density (kg/m3)

# Case configurations with colors
CASES = [
    {
        'name': 'Run1 (simpleFoam)',
        'path': ROOT_DIR / 'Re500_run1',
        'color': 'blue',
        'marker': 's',
        'linestyle': '-'
    },
    {
        'name': 'Run2 (simpleFoam)',
        'path': ROOT_DIR / 'Re500_run2',
        'color': 'green',
        'marker': '^',
        'linestyle': '--'
    },
    {
        'name': 'Transient (pimpleFoam)',
        'path': ROOT_DIR / 'Re500_transient',
        'color': 'red',
        'marker': 'o',
        'linestyle': ':'
    },
]


def parse_experimental_file(filepath):
    """Parse FDA experimental data file."""
    data_sections = {}
    current_key = None
    current_data = []
    expected_rows = 0
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('plot-') or line.startswith('geometry-') or line.startswith('fluid-') or line.startswith('dataset-'):
                if current_key and current_data:
                    data_sections[current_key] = np.array(current_data)
                current_key = line
                current_data = []
                expected_rows = 0
            elif current_key:
                parts = line.split()
                if len(parts) == 1 and parts[0].isdigit():
                    expected_rows = int(parts[0])
                elif len(parts) >= 2:
                    try:
                        row = [float(x) for x in parts[:2]]
                        current_data.append(row)
                    except ValueError:
                        pass
        
        if current_key and current_data:
            data_sections[current_key] = np.array(current_data)
    
    return data_sections


def read_openfoam_sample(postprocess_dir, set_name):
    """Read OpenFOAM sample data from postProcessing directory."""
    time_dirs = [d for d in postprocess_dir.iterdir() if d.is_dir()]
    if not time_dirs:
        return None, None, None
    
    time_dirs.sort(key=lambda x: float(x.name) if x.name.replace('.', '').isdigit() else 0)
    latest = time_dirs[-1]
    
    sample_file = latest / f"{set_name}_p_U.xy"
    if not sample_file.exists():
        return None, None, None
    
    data = np.loadtxt(sample_file)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    
    return data[:, 0], data[:, 1], data[:, 2:5]


def get_latest_time_dir(sim_data_dir):
    """Get the latest time directory from simulation_data."""
    if not sim_data_dir.exists():
        return None
    time_dirs = [d for d in sim_data_dir.iterdir() if d.is_dir()]
    if not time_dirs:
        return None
    time_dirs.sort(key=lambda x: float(x.name) if x.name.replace('.', '').isdigit() else 0)
    return time_dirs[-1]


def load_case_data(case):
    """Load simulation data directory for a case."""
    sim_data_dir = case['path'] / 'simulation_data'
    if not sim_data_dir.exists():
        return None
    return sim_data_dir


def plot_centerline_comparison(ax, cases_data, exp_data):
    """Plot centerline velocity comparison for all cases."""
    exp_centerline = None
    for key in exp_data:
        if 'z-distribution-axial-velocity' in key:
            exp_centerline = exp_data[key]
            break
    
    if exp_centerline is not None:
        ax.plot(exp_centerline[:, 0] * 1000, exp_centerline[:, 1], 
                'ko', markersize=6, label='Experiment (PIV)', zorder=10)
    
    for case, sim_data_dir in cases_data:
        if sim_data_dir is None:
            continue
        sim_pos, sim_p, sim_U = read_openfoam_sample(sim_data_dir, 'centerline')
        if sim_pos is not None and sim_U is not None:
            z_exp = (sim_pos - Z_EXPANSION) * 1000
            ax.plot(z_exp, sim_U[:, 2], 
                    color=case['color'], linestyle=case['linestyle'],
                    linewidth=2, label=case['name'])
    
    ax.set_xlabel('Axial position from expansion (mm)')
    ax.set_ylabel('Centerline Axial Velocity (m/s)')
    ax.set_title('Centerline Axial Velocity Comparison')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)


def plot_radial_profile_comparison(ax, cases_data, exp_data, z_exp_val, set_name, title):
    """Plot radial velocity profile comparison for all cases at a specific z."""
    exp_found = None
    for key in exp_data:
        if 'profile-axial-velocity-at-z' in key:
            match = re.search(r'at-z\s+(-?[\d.]+)', key)
            if match:
                z_in_key = float(match.group(1))
                if abs(z_in_key - z_exp_val) < 0.001:
                    exp_found = exp_data[key]
                    break
    
    if exp_found is not None:
        mask = exp_found[:, 1] != 0
        ax.plot(exp_found[mask, 0] * 1000, exp_found[mask, 1], 
                'ko', markersize=4, label='Experiment', zorder=10)
    
    for case, sim_data_dir in cases_data:
        if sim_data_dir is None:
            continue
        sim_pos, sim_p, sim_U = read_openfoam_sample(sim_data_dir, set_name)
        if sim_pos is not None and sim_U is not None:
            ax.plot(sim_pos * 1000, sim_U[:, 2], 
                    color=case['color'], linestyle=case['linestyle'],
                    linewidth=1.5, label=case['name'])
    
    ax.set_xlabel('r (mm)')
    ax.set_ylabel('Uz (m/s)')
    ax.set_title(title, fontsize=9)
    ax.grid(True, alpha=0.3)


def plot_pressure_comparison(ax, cases_data, exp_data):
    """Plot centerline pressure comparison for all cases."""
    exp_pressure = None
    for key in exp_data:
        if 'z-distribution-pressure' in key:
            exp_pressure = exp_data[key]
            break
    
    if exp_pressure is not None:
        ax.plot(exp_pressure[:, 0] * 1000, exp_pressure[:, 1], 
                'ko', markersize=6, label='Experiment', zorder=10)
    
    for case, sim_data_dir in cases_data:
        if sim_data_dir is None:
            continue
        sim_pos, sim_p_kin, sim_U = read_openfoam_sample(sim_data_dir, 'centerline')
        if sim_pos is not None and sim_p_kin is not None:
            z_exp = (sim_pos - Z_EXPANSION) * 1000
            sim_p = sim_p_kin * RHO
            
            if exp_pressure is not None:
                exp_at_zero = np.interp(0, exp_pressure[:, 0] * 1000, exp_pressure[:, 1])
                sim_at_zero = np.interp(0, z_exp, sim_p)
                sim_p_shifted = sim_p + (exp_at_zero - sim_at_zero)
            else:
                sim_p_shifted = sim_p
            
            ax.plot(z_exp, sim_p_shifted, 
                    color=case['color'], linestyle=case['linestyle'],
                    linewidth=2, label=case['name'])
    
    ax.set_xlabel('Axial position from expansion (mm)')
    ax.set_ylabel('Pressure (Pa)')
    ax.set_title('Centerline Pressure Comparison')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)


def calculate_jet_width(r, u_axial):
    """Calculate jet full width (diameter) at half maximum."""
    if r is None or u_axial is None or len(r) < 3:
        return None
    
    center_idx = np.argmin(np.abs(r))
    u_center = u_axial[center_idx]
    
    if u_center <= 0:
        return None
    
    u_half = u_center / 2
    positive_r_mask = r >= 0
    r_pos = r[positive_r_mask]
    u_pos = u_axial[positive_r_mask]
    
    for i in range(len(u_pos) - 1):
        if u_pos[i] >= u_half and u_pos[i+1] < u_half:
            t = (u_half - u_pos[i]) / (u_pos[i+1] - u_pos[i])
            r_half = r_pos[i] + t * (r_pos[i+1] - r_pos[i])
            return 2 * r_half
    
    return None


def plot_jet_width_comparison(ax, cases_data, exp_data):
    """Plot jet width comparison for all cases."""
    exp_jw = None
    for key in exp_data:
        if 'jet-width' in key:
            exp_jw = exp_data[key]
            break
    
    if exp_jw is not None:
        ax.plot(exp_jw[:, 0] * 1000, exp_jw[:, 1] * 1000, 
                'ko', markersize=6, label='Experiment', zorder=10)
    
    z_locations = [
        (-0.088, 'radial_z_minus088'),
        (-0.064, 'radial_z_minus064'),
        (-0.048, 'radial_z_minus048'),
        (-0.042, 'radial_z_minus042'),
        (-0.020, 'radial_z_minus020'),
        (-0.008, 'radial_z_minus008'),
        (0.000, 'radial_z_000'),
        (0.008, 'radial_z_plus008'),
        (0.016, 'radial_z_plus016'),
        (0.024, 'radial_z_plus024'),
        (0.032, 'radial_z_plus032'),
        (0.040, 'radial_z_plus040'),
        (0.048, 'radial_z_plus048'),
        (0.060, 'radial_z_plus060'),
        (0.080, 'radial_z_plus080'),
    ]
    
    for case, sim_data_dir in cases_data:
        if sim_data_dir is None:
            continue
        
        sim_z = []
        sim_jw = []
        
        for z_exp, set_name in z_locations:
            sim_pos, sim_p, sim_U = read_openfoam_sample(sim_data_dir, set_name)
            if sim_pos is not None and sim_U is not None:
                jw = calculate_jet_width(sim_pos, sim_U[:, 2])
                if jw is not None:
                    sim_z.append(z_exp * 1000)
                    sim_jw.append(jw * 1000)
        
        if sim_z:
            ax.plot(sim_z, sim_jw, 
                    color=case['color'], linestyle=case['linestyle'],
                    marker=case['marker'], markersize=5,
                    linewidth=2, label=case['name'])
    
    ax.set_xlabel('Axial position from expansion (mm)')
    ax.set_ylabel('Jet Width (mm)')
    ax.set_title('Jet Width Comparison')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)


def main():
    print("=" * 60)
    print("FDA Nozzle Benchmark - Multi-Case Comparison")
    print("=" * 60)
    
    # Load experimental data
    exp_file = CASES[0]['path'] / 'experimental_data' / 'PIV_Sudden_Expansion_500_243.txt'
    if not exp_file.exists():
        print(f"Error: Experimental file not found: {exp_file}")
        return
    
    exp_data = parse_experimental_file(exp_file)
    print(f"Loaded experimental data with {len(exp_data)} sections")
    
    # Load all case data
    cases_data = []
    for case in CASES:
        sim_data_dir = load_case_data(case)
        if sim_data_dir:
            time_dir = get_latest_time_dir(sim_data_dir)
            if time_dir:
                print(f"  ✓ {case['name']}: time = {time_dir.name}")
            else:
                print(f"  ✗ {case['name']}: No time directories found")
                sim_data_dir = None
        else:
            print(f"  ✗ {case['name']}: simulation_data not found")
        cases_data.append((case, sim_data_dir))
    
    # Create comparison plots
    fig = plt.figure(figsize=(16, 12))
    fig.suptitle('FDA Nozzle Benchmark V&V - Multi-Case Comparison\nRe = 500 (Laminar), Sudden Expansion', 
                 fontsize=14, fontweight='bold')
    
    ax1 = fig.add_subplot(2, 2, 1)
    plot_centerline_comparison(ax1, cases_data, exp_data)
    
    ax2 = fig.add_subplot(2, 2, 2)
    plot_pressure_comparison(ax2, cases_data, exp_data)
    
    ax3 = fig.add_subplot(2, 2, 3)
    plot_jet_width_comparison(ax3, cases_data, exp_data)
    
    ax4 = fig.add_subplot(2, 2, 4)
    plot_radial_profile_comparison(ax4, cases_data, exp_data, 0.0, 'radial_z_000', 
                                   'Radial Profile at z = 0 (Expansion)')
    ax4.legend(fontsize=8)
    
    plt.tight_layout()
    
    # Save to FDA_V&V directory
    output_file = ROOT_DIR / 'multi_case_comparison.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved comparison plot to: {output_file}")
    
    plt.close()
    print("Done!")


if __name__ == "__main__":
    main()
