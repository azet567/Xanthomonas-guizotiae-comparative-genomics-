#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import random
import warnings
warnings.filterwarnings('ignore')

print("========================================")
print("PANGENOME ACCUMULATION ANALYSIS (FINAL)")
print("========================================")

# Read the presence-absence matrix
df = pd.read_csv('presence_absence_matrix.csv', index_col=0)

print(f"Matrix dimensions: {df.shape}")

# Ensure correct orientation: genomes as columns, orthogroups as rows
# If rows > columns, we have orthogroups as rows (correct)
# If columns > rows, we need to transpose
if df.shape[0] > df.shape[1]:
    print("Orientation: Orthogroups as rows (correct)")
    binary = df.astype(int)
else:
    print("Orientation: Transposing to orthogroups as rows")
    binary = df.T.astype(int)

n_genomes = binary.shape[1]  # columns = genomes
n_orthogroups = binary.shape[0]  # rows = orthogroups

print(f"Genomes (columns): {n_genomes}")
print(f"Orthogroups (rows): {n_orthogroups}")

# Function to calculate pangenome accumulation
def pangenome_accumulation(matrix, n_permutations=1000):
    n_genomes = matrix.shape[1]  # columns
    n_orthogroups = matrix.shape[0]  # rows
    
    results = []
    
    for perm in range(n_permutations):
        if perm % 100 == 0:
            print(f"  Permutation {perm}/{n_permutations}")
        
        # Random order of genomes (columns)
        order = random.sample(range(n_genomes), n_genomes)
        cumulative = []
        current_set = set()
        
        for idx in range(n_genomes):
            genome_idx = order[idx]
            # Get orthogroups present in this genome (rows where column == 1)
            genes = set(matrix[matrix.iloc[:, genome_idx] == 1].index)
            current_set.update(genes)
            cumulative.append(len(current_set))
        
        results.append(cumulative)
    
    results = np.array(results)
    mean = np.mean(results, axis=0)
    ci_lower = np.percentile(results, 2.5, axis=0)
    ci_upper = np.percentile(results, 97.5, axis=0)
    
    return mean, ci_lower, ci_upper

# Run the accumulation
print("\nRunning pangenome accumulation with 1000 permutations...")
mean, ci_lower, ci_upper = pangenome_accumulation(binary, 1000)

# Print final values
final_mean = mean[-1]
print(f"\nFinal cumulative orthogroups: {final_mean:.0f}")

# Heaps' law fitting
x = np.arange(1, len(mean)+1)
log_x = np.log(x)
log_y = np.log(mean)

# Remove infinite values if any
finite_mask = np.isfinite(log_y)
log_x_finite = log_x[finite_mask]
log_y_finite = log_y[finite_mask]

model = LinearRegression()
model.fit(log_x_finite.reshape(-1, 1), log_y_finite.reshape(-1, 1))

alpha = model.coef_[0][0]
beta = np.exp(model.intercept_[0])
r2 = model.score(log_x_finite.reshape(-1, 1), log_y_finite.reshape(-1, 1))

print("\n========================================")
print("HEAPS' LAW PARAMETERS")
print("========================================")
print(f"α (alpha): {alpha:.6f}")
print(f"β (beta): {beta:.6f}")
print(f"R²: {r2:.6f}")

if alpha < 1:
    status = "OPEN pangenome"
else:
    status = "CLOSED pangenome"

print(f"\nINTERPRETATION: {status} (α {'<' if alpha < 1 else '>'} 1)")
if alpha < 1:
    print("The gene repertoire continues to expand as new genomes are sequenced.")

# Save statistics
with open('heaps_law_params_final.txt', 'w') as f:
    f.write(f"alpha\t{alpha:.6f}\n")
    f.write(f"beta\t{beta:.6f}\n")
    f.write(f"R2\t{r2:.6f}\n")
    f.write(f"Interpretation\t{status}\n")
    f.write(f"Permutations\t1000\n")
    f.write(f"Genomes\t{n_genomes}\n")
    f.write(f"Orthogroups\t{n_orthogroups}\n")
    f.write(f"Final_cumulative_orthogroups\t{final_mean:.0f}\n")

print("\nStatistics saved to: heaps_law_params_final.txt")

# Plot
print("\nGenerating figure...")
fig, ax = plt.subplots(figsize=(12, 8))

# Plot mean with confidence intervals
x_plot = range(1, len(mean)+1)
ax.plot(x_plot, mean, 'b-', linewidth=2.5, label='Mean (1,000 permutations)')
ax.fill_between(x_plot, ci_lower, ci_upper, alpha=0.25, color='blue', label='95% CI')

# Plot Heaps' law fit
x_fit = np.linspace(1, len(mean), 100)
y_fit = beta * (x_fit ** alpha)
ax.plot(x_fit, y_fit, 'r--', linewidth=2, label=f'Heaps\' law fit\nα = {alpha:.4f}, R² = {r2:.4f}')

ax.set_xlabel('Number of Genomes Sampled', fontsize=14, fontweight='bold')
ax.set_ylabel('Cumulative Number of Orthogroups', fontsize=14, fontweight='bold')
ax.set_title(f'Pangenome Accumulation Curve (1,000 Permutations)\nα = {alpha:.4f}, β = {beta:.4f}, R² = {r2:.4f}', fontsize=14)

if alpha < 1:
    status_label = "OPEN PANGENOME"
else:
    status_label = "CLOSED PANGENOME"
ax.text(0.98, 0.02, f"Status: {status_label}", transform=ax.transAxes, fontsize=12, fontweight='bold',
        verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax.legend(loc='lower right', fontsize=11)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('pangenome_accumulation_final.png', dpi=300, bbox_inches='tight')
plt.savefig('pangenome_accumulation_final.pdf', dpi=300, bbox_inches='tight')
print("Figure saved to: pangenome_accumulation_final.png and .pdf")

print("\n========================================")
print("COMPLETE")
print("========================================")
print(f"Final cumulative orthogroups: {final_mean:.0f}")
print(f"α (alpha): {alpha:.6f}")
print(f"β (beta): {beta:.6f}")
print(f"R²: {r2:.6f}")
print(f"Pangenome status: {status}")
