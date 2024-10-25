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
fc_threshold = 0.5  # Log2 fold change threshold
log_pvalue_threshold = -np.log10(0.05) # Adjusted p-value threshold

# Set up color-coding based on thresholds
data['color'] = np.where((data['log2FoldChange'] >= fc_threshold) & (data['-log10pvalue'] >= log_pvalue_threshold), 'red',
                np.where((data['log2FoldChange'] <= -fc_threshold) & (data['-log10pvalue'] >= log_pvalue_threshold), 'blue', 'grey'))

# Count upregulated and downregulated genes
num_upregulated = data[(data['log2FoldChange'] >= fc_threshold) & (data['-log10pvalue'] >= log_pvalue_threshold)].shape[0]
num_downregulated = data[(data['log2FoldChange'] <= -fc_threshold) & (data['-log10pvalue'] >= log_pvalue_threshold)].shape[0]

# Create the volcano plot
plt.figure(figsize=(8, 6))
plt.scatter(data['log2FoldChange'], data['-log10pvalue'], c=data['color'], alpha=0.75)
xmax = max(abs(data['log2FoldChange'])); plt.xlim(-xmax, xmax)
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

# Optionally, highlight genes with significant fold change and p-value
top_genes = data[(data['-log10pvalue'] > 1.3) & (abs(data['log2FoldChange']) > 0.5)]
for i, row in top_genes.iterrows():
    plt.text(row['log2FoldChange'], row['-log10pvalue'], row['Unnamed: 0'], fontsize=8)

# Show the plot
plt.tight_layout()
filename = input_file.strip(".csv")
plt.savefig(f"{filename}.png", dpi=400)
plt.close()
