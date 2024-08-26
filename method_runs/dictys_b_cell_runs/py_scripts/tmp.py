import pandas as pd
import numpy as np
import gzip

def compute_qc_metrics(expression_file):
    """
    Computes and prints QC metrics for a given expression matrix.

    Parameters:
    -----------
    expression_file : str
        Path to the input gzipped TSV file containing the expression matrix.
    """

    # Load the expression matrix
    with gzip.open(expression_file, 'rt') as f:
        expression_matrix = pd.read_csv(f, sep='\t', index_col=0)

    # Calculate metrics
    gene_totals = expression_matrix.sum(axis=1)  # Total read depth per gene (summing across cells)
    cell_totals = expression_matrix.sum(axis=0)  # Total coverage per cell (summing across genes)
    
    # Number of cells where each gene is expressed (non-zero)
    cells_per_gene = (expression_matrix > 0).sum(axis=1)
    
    # Proportion of cells where each gene is expressed
    cells_per_gene_percentage = (cells_per_gene / expression_matrix.shape[1]) * 100
    
    # Number of genes expressed in each cell (non-zero)
    genes_per_cell = (expression_matrix > 0).sum(axis=0)
    
    # Proportion of genes expressed in each cell
    genes_per_cell_percentage = (genes_per_cell / expression_matrix.shape[0]) * 100

    # Compute QC metrics
    min_total_read_depth_per_gene = gene_totals.min()
    min_cells_expressed_per_gene = cells_per_gene.min()
    min_cells_expressed_per_gene_percentage = cells_per_gene_percentage.min()
    min_total_coverage_per_cell = cell_totals.min()
    min_genes_expressed_per_cell = genes_per_cell.min()
    min_genes_expressed_per_cell_percentage = genes_per_cell_percentage.min()

    # Print QC metrics
    print(f"QC Metrics for {expression_file}:")
    print(f"Min total read depth per gene: {min_total_read_depth_per_gene}")
    print(f"Min number of cells expressed per gene: {min_cells_expressed_per_gene}")
    print(f"Min percentage of cells expressed per gene: {min_cells_expressed_per_gene_percentage:.2f}%")
    print(f"Min total coverage per cell: {min_total_coverage_per_cell}")
    print(f"Min number of genes expressed per cell: {min_genes_expressed_per_cell}")
    print(f"Min percentage of genes expressed per cell: {min_genes_expressed_per_cell_percentage:.2f}%")
    print("")

# Example usage for a single subset
subset_file = "tmp_static/Subset13/expression0.tsv.gz"
compute_qc_metrics(subset_file)
