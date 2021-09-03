#generating image means into a new dataframe
#Subomi Omoleye

import os
from argparse import ArgumentParser
from typing import Any, Callable, Dict, Literal, Optional, Tuple, Union
import cv2
import nibabel as nib
import numpy as np
import pandas as pd
from PIL import Image
from torchvision import transforms
print("all libs imported")

metadata_path= "/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/chimec_mammo_retry.csv" #corrected png paths
metadata = pd.read_csv(metadata_path)

print("loaded {} cases".format(
    metadata[metadata.case == True].study_id.nunique()))
print("loaded {} controls".format(
    metadata[metadata.case == False].study_id.nunique()))

columns_to_keep= ["study_id", "exam_id", "BurnedInAnnotation", "PresentationIntentType", "BreastImplantPresent",
"ContentDate", "ImageType", "RequestedProcedureDescription", "Laterality", "Modality", "SeriesDescription", 
"ImageLaterality", "ViewPosition", "ProtocolName", 'case','png_path', 'CompressionForce', "Modality", 
'EstimatedRadiographicMagnificationFactor', "FilterType", 'ReasonForStudy']

df_slim = metadata.filter(items= columns_to_keep, axis= 1)

#filterig for existing paths
df_slim["path_exists"] = df_slim["png_path"].apply(lambda x: os.path.isfile(x))
print("of {} total png paths, {} exist". format(len(df_slim), len(df_slim.loc[df_slim['path_exists']==True])))

df_true = df_slim[df_slim.path_exists==True]
print("created dataframe with only existing pngs")

df_true["image_means"]= df_true['png_path'].apply(lambda x: np.mean(np.asarray(Image.open(x))))
print("completed generation of image means")
print("converting dataframe to csv file")

df_true.to_csv("/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/data/true_pngs_mean.csv")