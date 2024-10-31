import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

cr_al_data = pd.read_csv('deseq2_mm10_cr-al.csv', header=None)
lipten_wt_data = pd.read_csv('deseq2_mm10_mut-wt.csv', header=None)

# Convert Log2FC and p-values to numeric, coercing errors to NaN
cr_al_data[2] = pd.to_numeric(cr_al_data[2], errors='coerce')
cr_al_data[5] = pd.to_numeric(cr_al_data[5], errors='coerce')
lipten_wt_data[2] = pd.to_numeric(lipten_wt_data[2], errors='coerce')
lipten_wt_data[5] = pd.to_numeric(lipten_wt_data[5], errors='coerce')

cr_al_sorted = cr_al_data.sort_values(by=[2, 5])
lipten_wt_sorted = lipten_wt_data.sort_values(by=[2, 5])

top_1000_cr_al = cr_al_sorted.head(1000)
bottom_1000_cr_al = cr_al_sorted.tail(1000)

# Identify downregulated genes (Group A) in CR vs. AL
group_a_genes = top_1000_cr_al[top_1000_cr_al[2] < 0][0].tolist()

# Identify upregulated genes (Group B) in LiPTEN vs. WT
group_b_genes = lipten_wt_sorted[lipten_wt_sorted[2] > 0][0].tolist()

group_a_set = set(group_a_genes)
group_b_set = set(group_b_genes)

plt.figure(figsize=(4, 3))
venn_diagram = venn2([group_a_set, group_b_set], ('CR vs AL downregulated', 'LiPTEN vs WT upregulated'))

# Set font properties for the diagram
for text in venn_diagram.set_labels:
    text.set_fontsize(6)
    text.set_family('Arial')

for text in venn_diagram.subset_labels:
    text.set_fontsize(6)
    text.set_family('Arial')

plt.savefig("venn_diagram.png", dpi=400, bbox_inches='tight')
plt.close()

overlapping_genes = group_a_set.intersection(group_b_set)
print("Overlapping Genes:")
for gene in overlapping_genes:
    print(gene)