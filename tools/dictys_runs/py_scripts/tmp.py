import os
import sys
import shutil
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.io import mmread

# Load the data
data_dir = "/ocean/projects/cis240075p/asachan/datasets/B_Cell/multiome_1st_donor_UPMC_aggr/outs/filtered_feature_bc_matrix/"

# load tsv files features and barcodes
features = pd.read_csv(data_dir + "features.tsv.gz", sep="\t", header=None)
barcodes = pd.read_csv(data_dir + "barcodes.tsv.gz", sep="\t", header=None)

# print the features head and the column 2 unique values
#print(features.head())
print(features[2].unique())

#select the Peaks rows and print their head
peaks = features[features[2] == "Peaks"]
print(peaks.head())