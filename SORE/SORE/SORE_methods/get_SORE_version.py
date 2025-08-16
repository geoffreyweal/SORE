"""
get_SORE_version.py, Geoffrey Weal, 26/3/22

get_SORE_version will return the version of the SORE program.
"""
import os

def get_SORE_version():
    """
    This method will grab the version of this SORE program without starting up the SORE process in full

    Returns
    -------
    The version of this SORE program
    """

    # First, get the path to the __init__file.
    path_to_init_file = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'../../__init__.py'))

    # Second, read through the __init__.py file and search for the __version__ line.
    version_no = str(None)
    with open(path_to_init_file, 'r') as FILE:
        for line in FILE:
            if '__version__ =' in line:
                line = line.rstrip().split()
                version_no = line[2].replace("'",'')
                break

    # Third, return the version number
    return version_no