# workflow_utils.py

"""
Utility functions to support workflow steps such as file discovery and data organization.
"""

import os

def find_all_bam_files(directory):
    """
    Finds all .bam files in the given directory.

    Parameters:
    -----------
    directory : str
        Path to the directory to search in.

    Returns:
    --------
    list of str
        List of full paths to .bam files found in the directory.
    """
    bam_files = []
    for filename in os.listdir(directory):
        if filename.endswith(".bam"):
            full_path = os.path.join(directory, filename)
            bam_files.append(full_path)
    
    return sorted(bam_files)
