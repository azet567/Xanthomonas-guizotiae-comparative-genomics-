#!/usr/bin/env python3
"""
Figure 3A: Orthogroup Frequency Distribution
CORRECTED DATA - Total = 7,140
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ============================================================
# CORRECT DATA - READ THIS CAREFULLY
# Total = 7,140 | Core = 1,043 | Soft-core = 230 | Accessory = 5,867
# ============================================================

orthogroups = [
    0,     # 1 genome
    151,   # 2 genomes
    1693,  # 3 genomes
    443,   # 4 genomes
    599,   # 5 genomes
    377,   # 6 genomes
    163,   # 7 genomes
    80,    # 8 genomes
    83,    # 9 genomes
    73,    # 10 genomes
    37,    # 11 genomes
    58,    # 12 genomes
    0,     # 13 genomes
    42,    # 14 genomes
    25,    # 15 genomes
    31,    # 16 genomes
    31,    # 17 genomes
    31,    # 18 genomes
    37,    # 19 genomes
    31,    # 20 genomes
    43,    # 21 genomes
    41,    # 22 genomes
    59,    # 23 genomes
    65,    # 24 genomes
    98,    # 25 genomes
    156,   # 26 genomes
    745,   # 27 genomes
    326,   # 28 genomes
    145,   # 29 genomes
    188,   # 30 genomes
    230,   # 31 genomes (soft-core)
    1043   # 32 genomes (core)
]

# ============================================================
# VERIFY THE DATA
# ============================================================

total = sum(orthogroups)
core = orthogroups[31]
soft = orthogroups[30]
accessory = total - core - soft

print("=" * 60)
print("FIGURE 3A - CORRECT DATA")
print("=" * 60)
print(f"Total orthogroups:          {total}")
print(f"Core (32/32):               {core} ({core/total*100:.1f}%)")
print(f"Soft-core (31/32):          {soft} ({soft/total*100:.1f}%)")
print(f"Accessory (2-30):           {accessory} ({accessory/total*100:.1f}%)")
print(f"Strain-specific (1/32):     {orthogroups[0]} (0.0%)")
print("=" * 60)

# ============================================================
# PLOTTING
# ============================================================

genomes = list(range(1, 33))

colors_core = '#D0021B'
colors_soft_core = '#F5A623'
colors_accessory = '#4A90D9'

bar_colors = []
for i, g in enumerate(genomes):
    if g == 32:
        bar_colors.append(colors_core)
    elif g == 31:
        bar_colors.append(colors_soft_core)
    else:
        bar_colors.append(colors_accessory)

fig, ax = plt.subplots(figsize=(14, 8))
fig.patch.set_facecolor('white')

bars = ax.bar(genomes, orthogroups, color=bar_colors,
              edgecolor='black', linewidth=0.8, alpha=0.85)

# Add value labels
for i, (bar, val) in enumerate(zip(bars, orthogroups)):
    if val > 100:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 20,
                f'{val:,}', ha='center', va='bottom', fontsize=8, fontweight='bold')

ax.set_xlabel('Number of Genomes Containing Orthogroup', fontsize=14, fontweight='bold')
ax.set_ylabel('Number of Orthogroups', fontsize=14, fontweight='bold')
ax.set_title('Orthogroup Frequency Distribution\nin 32 Xanthomonas Genomes',
             fontsize=16, fontweight='bold', pad=20)

ax.set_xticks(range(1, 33))
ax.set_xticklabels(range(1, 33), rotation=45, ha='right', fontsize=9)

ax.grid(True, alpha=0.25, linestyle='-', linewidth=0.5, axis='y')
ax.set_axisbelow(True)

ax.axvline(x=30.5, color='gray', linestyle='--', linewidth=1.5, alpha=0.6)
ax.axvline(x=31.5, color='gray', linestyle='--', linewidth=1.5, alpha=0.6)

ax.annotate(f'CORE GENOME\n{core:,} orthogroups\n({core/total*100:.1f}%)',
            xy=(32, core), xytext=(31.5, core + 350),
            arrowprops=dict(arrowstyle='->', color=colors_core, lw=2),
            fontsize=11, fontweight='bold', color=colors_core,
            ha='right', va='bottom')

ax.annotate(f'SOFT-CORE\n{soft:,} orthogroups\n({soft/total*100:.1f}%)',
            xy=(31, soft), xytext=(30.5, soft + 300),
            arrowprops=dict(arrowstyle='->', color=colors_soft_core, lw=2),
            fontsize=10, fontweight='bold', color=colors_soft_core,
            ha='right', va='bottom')

legend_elements = [
    mpatches.Patch(color=colors_core, label=f'Core (32/32): {core:,} ({core/total*100:.1f}%)'),
    mpatches.Patch(color=colors_soft_core, label=f'Soft-core (31/32): {soft:,} ({soft/total*100:.1f}%)'),
    mpatches.Patch(color=colors_accessory, label=f'Accessory (2-30): {accessory:,} ({accessory/total*100:.1f}%)'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=11, frameon=True,
          fancybox=True, edgecolor='black')

ax.text(0.02, 0.98, f'Total: {total:,} orthogroups',
        transform=ax.transAxes, fontsize=12, fontweight='bold',
        verticalalignment='top', horizontalalignment='left',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black'))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)

plt.tight_layout()

# ============================================================
# SAVE WITH NEW FILENAME TO CONFIRM
# ============================================================

plt.savefig('Figure_3A_CORRECTED.png', dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig('Figure_3A_CORRECTED.pdf', bbox_inches='tight', facecolor='white')
plt.savefig('Figure_3A_CORRECTED.svg', format='svg', bbox_inches='tight', facecolor='white')

# Also overwrite the old files
plt.savefig('Figure_3A_Orthogroup_Frequency.png', dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig('Figure_3A_Orthogroup_Frequency.pdf', bbox_inches='tight', facecolor='white')
plt.savefig('Figure_3A_Orthogroup_Frequency.svg', format='svg', bbox_inches='tight', facecolor='white')

print("\n" + "=" * 60)
print("✓ FIGURE 3A UPDATED WITH CORRECT DATA!")
print("=" * 60)
print("New files created:")
print("  - Figure_3A_CORRECTED.png")
print("  - Figure_3A_CORRECTED.pdf")
print("  - Figure_3A_CORRECTED.svg")
print("\nOld files overwritten:")
print("  - Figure_3A_Orthogroup_Frequency.png")
print("  - Figure_3A_Orthogroup_Frequency.pdf")
print("  - Figure_3A_Orthogroup_Frequency.svg")
print("=" * 60)

