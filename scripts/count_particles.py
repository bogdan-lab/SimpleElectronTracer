import numpy as np


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()
    
    
    total = 0
    for fname in args.files:
        data = np.loadtxt(fname)
        print("%s --> %i particles" % (fname, len(data)))
        total += len(data)
    print("TOTAL PARTICLES --> %i" % total)