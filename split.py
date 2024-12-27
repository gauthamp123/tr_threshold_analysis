import pandas as pd
import numpy as np
import random

def define_splits(bed_file, train_percent, test_percent, val_percent):
    with open(bed_file, 'r') as infile, open('data/prepared.bed', "w") as outfile:
        for line in infile:
            # Skip empty lines
            if line.strip():
                rand_val = random.random()
                if rand_val < train_percent:
                    split = "train"
                elif rand_val < train_percent + test_percent:
                    split = "test"
                else:
                    split = "valid"
                # Write the line with the new column
                outfile.write(f"{line.strip()}\t{split}\n")
        

define_splits('data/hg38_ver13.bed', 0.7, 0.2, 0.1)