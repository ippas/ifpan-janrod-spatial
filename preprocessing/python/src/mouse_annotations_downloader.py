from pybiomart import Dataset
import pandas as pd

def get_mouse_annotations(biomart_version="v102"):
    """
    Returns mouse gene annotations (GRCm38/mm10 or GRCm39/mm39) 
    for a given Ensembl BioMart version, e.g., biomart_version = "v102", "v105", "v110".

    Supported versions: v99–v110 (easily extendable).

    Returns:
        pd.DataFrame: columns: gene_id, gene_name, chromosome, start, end, strand, biotype, description
    """
    # Map of supported Ensembl BioMart versions and their archive hosts
    ensembl_hosts = {
        "v99": "http://jan2020.archive.ensembl.org",
        "v100": "http://apr2020.archive.ensembl.org",
        "v101": "http://aug2020.archive.ensembl.org",
        "v102": "http://nov2020.archive.ensembl.org",
        "v103": "http://feb2021.archive.ensembl.org",
        "v104": "http://may2021.archive.ensembl.org",
        "v105": "http://dec2021.archive.ensembl.org",
        "v106": "http://apr2022.archive.ensembl.org",
        "v107": "http://jul2022.archive.ensembl.org",
        "v108": "http://oct2022.archive.ensembl.org",
        "v109": "http://feb2023.archive.ensembl.org",
        "v110": "http://jul2023.archive.ensembl.org"
    }

    # Check if the requested version is supported
    if biomart_version not in ensembl_hosts:
        raise ValueError(f"Unknown BioMart version: {biomart_version}. Supported versions: {list(ensembl_hosts.keys())}")

    # Get the corresponding host URL
    host = ensembl_hosts[biomart_version]

    # Connect to the BioMart dataset for mouse (Mus musculus)
    dataset = Dataset(name='mmusculus_gene_ensembl', host=host)

    # Query extended gene annotations
    df = dataset.query(attributes=[
        'ensembl_gene_id',
        'external_gene_name',
        'chromosome_name',
        'start_position',
        'end_position',
        'strand',
        'gene_biotype',
        'description'
    ])

    # Rename columns for clarity
    df.columns = ['gene_id', 'gene_name', 'chromosome', 'start', 'end', 'strand', 'biotype', 'description']

    print(f"✅ Data retrieved from Ensembl {biomart_version} ({host})")
    return df