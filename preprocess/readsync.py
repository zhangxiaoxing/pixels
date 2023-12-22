import sys
import numpy as np
from struct import unpack


def readsync(filename):
    """
    Reads the last 2 bytes of each 770-byte chunk in a file as an int16,
    saves them in a NumPy array, and saves the array to output.npy.

    Args:
      filename: The path to the file to read.
    """
    chunk_size = 770
    data = []

    with open(filename, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            # Extract the last 2 bytes as int16
            value = unpack("<h", chunk[-2:])[0]
            data.append(value)

    # Convert data to NumPy array and save
    data = np.array(data, dtype=np.int8)
    np.save("sync_raw.npy", data)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        print(f"Processing file: {filename}")
        readsync(filename)
    else:
        print("Please provide a filename as an argument.")
