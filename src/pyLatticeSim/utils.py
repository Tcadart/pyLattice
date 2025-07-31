import os

def clear_directory(directoryPath):
    """
    Clear all files in a directory

    Parameters:
    ------------
    directoryPath: string
        Path to the directory to clear
    """
    # List all files in the directory
    files = [file for file in os.listdir(directoryPath) if os.path.isfile(os.path.join(directoryPath, file))]

    # Loop through and delete each file
    for file in files:
        os.remove(os.path.join(directoryPath, file))