#to check if png path exists and return a dataframe/csv with a new column for this

import glob
import logging
import os
import random
import subprocess
import time
import re
import numpy as np
import pandas as pd

file_path = "/scratch/ojomoleye/mirai-validation/magnified_mammos.csv"
#change to my metadata path
data = pd.read_csv(file_path)
df=pd.DataFrame(data)

print("number of magnified pngs is {}".format(len(df)))

#check if png file exists
df["path_exists"] = df["png_path"].apply(lambda x: os.path.isfile(x))

df.to_csv(path_or_buf="/scratch/ojomoleye/mirai-validation/magnified_mammos_with_paths.csv") 

print(" of {} total magnified pngs, {} were converted".format(len(df), len(df[df.path_exists==True])))