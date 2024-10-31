import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# Ensure the user provides a file as argument
if len(sys.argv) < 2:
    print("Usage: python script.py <input_file.csv>")
    sys.exit(1)

# Read the CSV file passed via command-line argument
input_file = sys.argv[1]
data = pd.read_csv(input_file)

# Check if the required columns are present in the input CSV
required_columns = ['log2FoldChange', 'pvalue', 'baseMean']
if not all(col in data.columns for col in required_columns):
    print(f"Missing required columns in the input file. Expected columns: {required_columns}")
    sys.exit(1)

# Add a column for -log10(pvalue)
data['-log10pvalue'] = -np.log10(data['pvalue'])

# Define thresholds
fc_threshold = 1.0  # Log2 fold change threshold
log_pvalue_threshold = -np.log10(0.05) # Adjusted p-value threshold

# Set up color-coding based on thresholds
data['color'] = np.where((data['log2FoldChange'] >= fc_threshold) & (data['-log10pvalue'] >= log_pvalue_threshold), 'red',
                np.where((data['log2FoldChange'] <= -fc_threshold) & (data['-log10pvalue'] >= log_pvalue_threshold), 'blue', 'grey'))

# Count upregulated and downregulated genes
num_upregulated = data[(data['log2FoldChange'] >= fc_threshold) & (data['-log10pvalue'] >= log_pvalue_threshold)].shape[0]
num_downregulated = data[(data['log2FoldChange'] <= -fc_threshold) & (data['-log10pvalue'] >= log_pvalue_threshold)].shape[0]


# Create the volcano plot
plt.figure(figsize=(10, 8))
plt.scatter(data['log2FoldChange'], data['-log10pvalue'], c=data['color'], alpha=0.75)

# Set manual xmax value if needed.
xmax = 10; plt.xlim(-xmax, xmax)
#xmax = max(abs(data['log2FoldChange'])); plt.xlim(-xmax, xmax)
ymax = max(data['-log10pvalue']); plt.ylim(0, ymax + 0.5)

# Add lines for thresholds
plt.axhline(y=log_pvalue_threshold, color='lightgray', linestyle='--', lw=0.50)
plt.axvline(x=fc_threshold, color='lightgray', linestyle='--', lw=0.50)
plt.axvline(x=-fc_threshold, color='lightgray', linestyle='--', lw=0.50)

# Labeling the plot
plt.xlabel('Log2 Fold Change', fontsize=10)
plt.ylabel('-Log10 p-value', fontsize=10)

# Annotate the number of upregulated and downregulated genes on the plot
plt.text(xmax * 0.65, ymax, f'Upregulated: {num_upregulated}', fontsize=8, color='red')
plt.text(-xmax * 0.95, ymax, f'Downregulated: {num_downregulated}', fontsize=8, color='blue')

# Filter for significant genes based on thresholds
significant_genes = data[(abs(data['log2FoldChange']) >= fc_threshold) & (data['-log10pvalue'] >= log_pvalue_threshold)]

# Select top 50 by fold change in each direction (up and down)
top_log2fc_up = significant_genes[significant_genes['log2FoldChange'] > 0].nlargest(50, 'log2FoldChange')
top_log2fc_down = significant_genes[significant_genes['log2FoldChange'] < 0].nsmallest(50, 'log2FoldChange')

# Select top 50 by -log10(pvalue) (most significant p-values)
top_pvalue = significant_genes.nlargest(50, '-log10pvalue')

# Combine the unique genes from both criteria
top_genes = pd.concat([top_log2fc_up, top_log2fc_down, top_pvalue]).drop_duplicates()

# Label the selected top genes
for i, row in top_genes.iterrows():
    plt.text(row['log2FoldChange'], row['-log10pvalue'], row['Unnamed: 0'], fontsize=8)

# Show the plot
plt.tight_layout()
filename = input_file.strip(".csv")
plt.savefig(f"{filename}_top50.png", dpi=400)
plt.close()