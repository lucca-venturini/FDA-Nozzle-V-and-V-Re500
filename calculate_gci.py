#!/usr/bin/env python3
"""
GCI (Grid Convergence Index) Analysis
Per Roache/ASME V&V 20-2009 methodology

Calculates discretization uncertainty for FDA Nozzle Benchmark.
Generates convergence plots for documentation.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Constants
Z_EXPANSION = 0.122685  # Expansion plane in simulation coordinates
RHO = 1056.0
PLOTS_DIR = Path(__file__).parent / 'GCI' / 'plots'

def read_openfoam_sample(case_dir, set_name):
    """Read sample data from postProcessing or simulation_data."""
    # Try simulation_data first
    sim_data = case_dir / 'simulation_data'
    if sim_data.exists():
        time_dirs = sorted([d for d in sim_data.iterdir() if d.is_dir()], 
                          key=lambda x: float(x.name) if x.name.replace('.','').isdigit() else 0)
        if time_dirs:
            sample_file = time_dirs[-1] / f"{set_name}_p_U.xy"
            if sample_file.exists():
                data = np.loadtxt(sample_file)
                return data[:, 0], data[:, 1], data[:, 2:5]
    
    # Try postProcessing
    post = case_dir / 'postProcessing' / 'sampleDict'
    if post.exists():
        time_dirs = sorted([d for d in post.iterdir() if d.is_dir()],
                          key=lambda x: float(x.name) if x.name.replace('.','').isdigit() else 0)
        if time_dirs:
            sample_file = time_dirs[-1] / f"{set_name}_p_U.xy"
            if sample_file.exists():
                data = np.loadtxt(sample_file)
                return data[:, 0], data[:, 1], data[:, 2:5]
    
    return None, None, None


def get_centerline_velocity_at_z0(case_dir):
    """Get centerline axial velocity at z=0 (expansion plane)."""
    pos, p, U = read_openfoam_sample(case_dir, 'centerline')
    if pos is None:
        return None
    
    # Find z=0 in experimental coordinates
    z_exp = pos - Z_EXPANSION
    idx = np.argmin(np.abs(z_exp))
    return U[idx, 2]  # Uz at expansion


def get_pressure_drop(case_dir):
    """Get pressure drop from inlet to outlet centerline."""
    pos, p_kin, U = read_openfoam_sample(case_dir, 'centerline')
    if pos is None:
        return None
    
    p = p_kin * RHO  # Convert to static pressure
    
    # Pressure at inlet (z=0 in sim) and outlet (z=0.242685)
    inlet_idx = np.argmin(np.abs(pos - 0.0))
    outlet_idx = np.argmin(np.abs(pos - 0.242685))
    
    return p[inlet_idx] - p[outlet_idx]


def get_jet_width_at_z32(case_dir):
    """Get jet width at z=+32mm from expansion."""
    pos, p, U = read_openfoam_sample(case_dir, 'radial_z_plus032')
    if pos is None:
        return None
    
    # Calculate jet half-width
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
            return 2 * r_half * 1000  # Full width in mm
    
    return None


def calculate_gci(f1, f2, f3, r21, r32):
    """
    Calculate GCI using Roache methodology.
    
    f1 = fine mesh result
    f2 = medium mesh result
    f3 = coarse mesh result
    r21 = h2/h1 (refinement ratio medium/fine)
    r32 = h3/h2 (refinement ratio coarse/medium)
    
    Returns: (p, f_exact, GCI_fine, GCI_medium)
    """
    # Calculate apparent order of convergence
    epsilon_32 = f3 - f2
    epsilon_21 = f2 - f1
    
    if epsilon_21 == 0 or epsilon_32 == 0:
        return None, None, None, None
    
    # Check for oscillatory convergence
    s = np.sign(epsilon_32 / epsilon_21)
    
    if s < 0:
        print("  Warning: Oscillatory convergence detected")
    
    # Iterative solution for p (order of convergence)
    # Using simplified formula for constant r
    r = r21  # Assuming r21 = r32 = r
    
    # p = ln(|epsilon_32/epsilon_21|) / ln(r)
    if abs(epsilon_32/epsilon_21) <= 0:
        return None, None, None, None
        
    p = abs(np.log(abs(epsilon_32/epsilon_21))) / np.log(r)
    
    # Limit p to reasonable range (0.5 to 5)
    p = max(0.5, min(5.0, p))
    
    # Richardson extrapolation
    f_exact = f1 + (f1 - f2) / (r**p - 1)
    
    # GCI calculation (Fs = 1.25 for 3 grids)
    Fs = 1.25
    
    # Relative errors
    e_21 = abs((f1 - f2) / f1) if f1 != 0 else 0
    e_32 = abs((f2 - f3) / f2) if f2 != 0 else 0
    
    GCI_fine = Fs * e_21 / (r**p - 1) * 100  # as percentage
    GCI_medium = Fs * e_32 / (r**p - 1) * 100
    
    return p, f_exact, GCI_fine, GCI_medium


def main():
    print("=" * 70)
    print("GCI Analysis - FDA Nozzle Benchmark (Re = 500)")
    print("Methodology: Roache/ASME V&V 20-2009")
    print("=" * 70)
    
    gci_dir = Path(__file__).parent / 'GCI'
    
    coarse_dir = gci_dir / 'coarse'
    medium_dir = gci_dir / 'medium'
    fine_dir = gci_dir / 'fine'
    
    # Refinement ratio (cell size ratio, not count)
    r = 2.0  # h_coarse/h_medium = h_medium/h_fine = 2
    
    print(f"\nRefinement ratio: r = {r}")
    print(f"Cell count ratio: r³ = {r**3:.1f}")
    
    # Extract quantities from each mesh
    quantities = {
        'Centerline Velocity at z=0 (m/s)': {
            'func': get_centerline_velocity_at_z0,
            'unit': 'm/s'
        },
        'Pressure Drop (Pa)': {
            'func': get_pressure_drop,
            'unit': 'Pa'
        },
        'Jet Width at z=+32mm': {
            'func': get_jet_width_at_z32,
            'unit': 'mm'
        }
    }
    
    print("\n" + "-" * 70)
    print("Extracted Values:")
    print("-" * 70)
    print(f"{'Quantity':<35} {'Coarse':>12} {'Medium':>12} {'Fine':>12}")
    print("-" * 70)
    
    results = {}
    
    for name, info in quantities.items():
        f3 = info['func'](coarse_dir)  # Coarse
        f2 = info['func'](medium_dir)  # Medium
        f1 = info['func'](fine_dir)    # Fine
        
        results[name] = {'f1': f1, 'f2': f2, 'f3': f3, 'unit': info['unit']}
        
        f3_str = f"{f3:.6f}" if f3 is not None else "N/A"
        f2_str = f"{f2:.6f}" if f2 is not None else "N/A"
        f1_str = f"{f1:.6f}" if f1 is not None else "N/A"
        
        print(f"{name:<35} {f3_str:>12} {f2_str:>12} {f1_str:>12}")
    
    print("\n" + "=" * 70)
    print("GCI Results:")
    print("=" * 70)
    print(f"{'Quantity':<30} {'p':>6} {'f_exact':>12} {'GCI_fine':>10} {'GCI_med':>10}")
    print("-" * 70)
    
    for name, data in results.items():
        f1, f2, f3 = data['f1'], data['f2'], data['f3']
        
        if f1 is None or f2 is None or f3 is None:
            print(f"{name:<30} {'N/A':>6} {'N/A':>12} {'N/A':>10} {'N/A':>10}")
            continue
        
        p, f_exact, gci_fine, gci_med = calculate_gci(f1, f2, f3, r, r)
        
        if p is None:
            print(f"{name:<30} {'N/A':>6} {'N/A':>12} {'N/A':>10} {'N/A':>10}")
            continue
        
        print(f"{name:<30} {p:>6.2f} {f_exact:>12.6f} {gci_fine:>9.2f}% {gci_med:>9.2f}%")
        
        # Store for summary
        data['p'] = p
        data['f_exact'] = f_exact
        data['gci_fine'] = gci_fine
    
    print("\n" + "=" * 70)
    print("Summary:")
    print("=" * 70)
    print("""
GCI_fine represents the estimated discretization uncertainty band (95% confidence)
for the fine mesh solution.

Interpretation:
- GCI < 2%:  Excellent - mesh is well-converged
- GCI 2-5%:  Good - acceptable for most engineering applications
- GCI 5-10%: Moderate - consider refinement for critical applications
- GCI > 10%: Poor - mesh refinement recommended
""")

    # Generate plots
    generate_gci_plots(results, r)


def generate_gci_plots(results, r):
    """Generate and save GCI convergence plots."""
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Plot 1: Grid convergence for all quantities
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    
    h_ratio = [2.0, 1.0, 0.5]  # Relative cell size (coarse, medium, fine)
    mesh_labels = ['Coarse\n(~250k)', 'Medium\n(~2M)', 'Fine\n(~16M)']
    
    plot_data = [
        ('Centerline Velocity at z=0 (m/s)', 'Centerline Velocity', 'm/s', 'steelblue'),
        ('Pressure Drop (Pa)', 'Pressure Drop', 'Pa', 'darkorange'),
        ('Jet Width at z=+32mm', 'Jet Width', 'mm', 'seagreen')
    ]
    
    for idx, (key, title, unit, color) in enumerate(plot_data):
        ax = axes[idx]
        if key in results:
            data = results[key]
            values = [data['f3'], data['f2'], data['f1']]  # coarse, medium, fine
            
            ax.plot(h_ratio, values, 'o-', color=color, markersize=10, linewidth=2)
            
            # Richardson extrapolated value
            if 'f_exact' in data:
                ax.axhline(y=data['f_exact'], color=color, linestyle='--', alpha=0.5,
                          label=f'Extrapolated: {data["f_exact"]:.4f}')
            
            # Error bars for GCI
            if 'gci_fine' in data:
                gci_err = data['gci_fine'] / 100 * data['f1']
                ax.errorbar([0.5], [data['f1']], yerr=[gci_err], color=color, 
                           capsize=5, capthick=2, fmt='none')
            
            ax.set_xlabel('Relative Cell Size (h/h_fine)', fontsize=11)
            ax.set_ylabel(f'{title} [{unit}]', fontsize=11)
            ax.set_title(title, fontsize=12, fontweight='bold')
            
            # Add mesh labels
            for i, (h, val) in enumerate(zip(h_ratio, values)):
                ax.annotate(mesh_labels[i], (h, val), textcoords="offset points",
                           xytext=(0, 12), ha='center', fontsize=8)
            
            ax.grid(True, alpha=0.3)
            ax.set_xlim(0.3, 2.3)
            if 'f_exact' in data:
                ax.legend(loc='best', fontsize=9)
    
    plt.suptitle('Grid Convergence Study (Roache GCI, r=2.0)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plot_path = PLOTS_DIR / 'gci_convergence.png'
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {plot_path}")
    
    # Plot 2: GCI bar chart
    fig, ax = plt.subplots(figsize=(8, 5))
    
    quantities = []
    gci_fine = []
    gci_medium = []
    
    for key, data in results.items():
        if 'gci_fine' in data:
            short_name = key.split('(')[0].strip()
            quantities.append(short_name)
            gci_fine.append(data['gci_fine'])
    
    x = np.arange(len(quantities))
    bars = ax.bar(x, gci_fine, color='steelblue', edgecolor='black', alpha=0.8)
    
    # Threshold lines
    ax.axhline(y=2.0, color='green', linestyle='--', linewidth=1.5, label='Excellent (<2%)')
    ax.axhline(y=5.0, color='orange', linestyle='--', linewidth=1.5, label='Good (<5%)')
    
    ax.set_xlabel('Quantity', fontsize=12)
    ax.set_ylabel('GCI Fine Mesh [%]', fontsize=12)
    ax.set_title('Grid Convergence Index (95% Confidence)', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(quantities, rotation=10, ha='right', fontsize=10)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for bar, val in zip(bars, gci_fine):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
               f'{val:.2f}%', ha='center', fontsize=10, fontweight='bold')
    
    ax.set_ylim(0, max(gci_fine) * 1.3)
    
    plot_path = PLOTS_DIR / 'gci_uncertainty.png'
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {plot_path}")


if __name__ == "__main__":
    main()
