#!/usr/bin/env python3
"""
Modify blockMeshDict cell counts for GCI study.
Preserves grading ratios and geometry - only scales cell counts.
"""
import re
import sys
from pathlib import Path

def modify_blockmesh(input_file, output_file, scale_factor):
    """Modify cell counts in blockMeshDict by scale_factor."""
    
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Split into lines for processing
    lines = content.split('\n')
    modified_lines = []
    
    in_blocks_section = False
    
    for i, line in enumerate(lines):
        # Detect blocks section
        if line.strip() == 'blocks':
            in_blocks_section = True
        elif in_blocks_section and line.strip().startswith('edges'):
            in_blocks_section = False
        
        # Only modify cell count lines that look like "    (4 4 30)"
        # These come after hex definitions and contain only 3 integers
        if in_blocks_section:
            # Match lines with exactly 3 integers in parentheses (cell counts)
            match = re.match(r'^(\s+)\((\d+)\s+(\d+)\s+(\d+)\)(.*)$', line)
            if match:
                indent = match.group(1)
                n1 = max(2, int(int(match.group(2)) * scale_factor))
                n2 = max(2, int(int(match.group(3)) * scale_factor))
                n3 = max(2, int(int(match.group(4)) * scale_factor))
                suffix = match.group(5)
                line = f"{indent}({n1} {n2} {n3}){suffix}"
        
        modified_lines.append(line)
    
    with open(output_file, 'w') as f:
        f.write('\n'.join(modified_lines))
    
    print(f"Modified {input_file} -> {output_file} (scale={scale_factor})")

if __name__ == "__main__":
    base_dir = Path(__file__).parent
    
    medium_file = base_dir / "GCI/medium/system/blockMeshDict"
    
    # Create coarse (half cells)
    coarse_file = base_dir / "GCI/coarse/system/blockMeshDict"
    modify_blockmesh(medium_file, coarse_file, 0.5)
    
    # Create fine (double cells) 
    fine_file = base_dir / "GCI/fine/system/blockMeshDict"
    modify_blockmesh(medium_file, fine_file, 2.0)
    
    print("\nDone! Created coarse and fine blockMeshDict files.")
