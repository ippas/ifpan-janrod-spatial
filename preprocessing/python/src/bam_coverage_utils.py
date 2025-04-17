import pysam
import pandas as pd

def get_coverage(bam_path, chrom, start, end):
    """
    Returns the per-base coverage for a specified region (chrom:start-end) from a BAM file.
    
    Args:
        bam_path (str): Path to the BAM file.
        chrom (str): Chromosome name, e.g., "chr1".
        start (int): Start position of the region (0-based).
        end (int): End position of the region (0-based, exclusive).
    
    Returns:
        pd.DataFrame: DataFrame with columns ['chrom', 'position', 'depth'].
    """
    # Open the BAM file
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    
    # Retrieve coverage information for the specified region
    coverage = bamfile.count_coverage(chrom, start, end)
    
    # Sum the coverage across all nucleotide bases (A, C, G, T)
    total_coverage = [sum(bases) for bases in zip(*coverage)]
    
    # Close the BAM file
    bamfile.close()

    # Generate 1-based positions for the specified range
    positions = list(range(start + 1, end + 1))
    
    # Create a DataFrame with the coverage information
    df = pd.DataFrame({
        'chrom': chrom,
        'position': positions,
        'depth': total_coverage
    })
    
    return df


def get_per_base_coverage(bam_path, chrom, start, end):
    """
    Zwraca pokrycie (liczbę odczytów) dla każdej pozycji w zadanym zakresie (chrom:start-end) z pliku BAM.
    
    Args:
        bam_path (str): Ścieżka do pliku BAM.
        chrom (str): Nazwa chromosomu, np. "chr1".
        start (int): Początek zakresu (0-based inclusive).
        end (int): Koniec zakresu (0-based exclusive).
    
    Returns:
        pd.DataFrame: DataFrame z kolumnami ['chrom', 'position', 'depth'].
                      'position' jest 1-based.
    """
    # Otwórz plik BAM
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        # Pobierz pokrycie dla zadanego regionu
        coverage = bamfile.count_coverage(chrom, start, end)
        
        # Suma pokrycia dla wszystkich nukleotydów (A, C, G, T)
        total_coverage = [sum(base_counts) for base_counts in zip(*coverage)]
        
        # Utwórz DataFrame z wynikami
        positions = list(range(start + 1, end + 1))  # Konwersja na 1-based
        df = pd.DataFrame({
            'chrom': chrom,
            'position': positions,
            'depth': total_coverage
        })
    
    return df

def extract_max_coverage_for_genes(genes, files, gene_annotation):
    """
    Extracts the maximum per-base coverage for a given list of genes across multiple BAM files.

    Parameters:
    -----------
    genes : list of str
        List of gene names (symbols) to extract coverage for.

    files : list of str
        List of paths to BAM files (one per sample), for which the coverage should be extracted.

    gene_annotation : pandas.DataFrame
        DataFrame containing gene annotation information.
        Must contain at least the following columns: 'gene_name', 'start', 'end', 'chromosome'.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame with columns: 'sample_id', 'gene_symbol', 'max_coverage',
        where each row corresponds to the maximum coverage value for a given gene in a given sample.
    """
    results = []

    for bam_path in files:
        # Extract sample ID from file path (e.g. from "S13839Nr3.bam" to "S13839Nr3")
        sample_id = bam_path.split('/')[-1].replace('.bam', '')
        print(f"\nProcessing: {sample_id}")
        
        for gene in genes:
            # Find gene annotation row for current gene
            row = gene_annotation[gene_annotation['gene_name'] == gene]
            
            if row.empty:
                print(f"WARNING - Gene {gene} not found in annotation.")
                results.append({
                    "sample_id": sample_id,
                    "gene_symbol": gene,
                    "max_coverage": None
                })
                continue
            
            # Extract genomic coordinates
            gene_start = int(row['start'].values[0])
            gene_end = int(row['end'].values[0])
            gene_chrom = str(row['chromosome'].values[0])
            
            # Add 'chr' prefix if missing
            if not gene_chrom.startswith('chr'):
                gene_chrom = f"chr{gene_chrom}"

            try:
                # Get per-base coverage for the gene region
                gene_coverage = get_per_base_coverage(
                    bam_path=bam_path,
                    chrom=gene_chrom,
                    start=gene_start,
                    end=gene_end
                )
                # Extract the maximum depth value
                max_coverage = gene_coverage['depth'].max()

                results.append({
                    "sample_id": sample_id,
                    "gene_symbol": gene,
                    "max_coverage": max_coverage
                })

            except Exception as e:
                # Handle errors (e.g. if gene not found in BAM)
                print(f"ERROR - {gene} in {sample_id}: {e}")
                results.append({
                    "sample_id": sample_id,
                    "gene_symbol": gene,
                    "max_coverage": None
                })

    return pd.DataFrame(results)

