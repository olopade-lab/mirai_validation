# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import glob
import logging
import os
import random
import subprocess
import time
import re
import numpy as np
import pandas as pd


# %%
metadata_path = "/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/chimec_mammo_series_metadata.csv"
#change to my metadata path
metadata = pd.read_csv(metadata_path)

print("loaded {} cases".format(
    metadata[metadata.case == True].study_id.nunique()))
print("loaded {} controls".format(
    metadata[metadata.case == False].study_id.nunique()))


# %%
#11000/54353

df= pd.DataFrame(metadata)
# df_sample= df.truncate(after=11012) #truncate based on stage at png creation for testing

df
#list columns
list_of_col= []
for col in df.columns:
    list_of_col.append(col)
# print(list_of_col, "this is the number of colums:", len(list_of_col))

# %% [markdown]
# Yala CSV input
# columns: patient_id, exam_id, laterality, view, file_path(png), years_to_cancer, years_to_last_followup, split_group

# %%
#selecting columns to keep
columns_to_keep= ["study_id", "exam_id", "BurnedInAnnotation", "PresentationIntentType", "BreastImplantPresent",
"ContentDate", "ContentTime", "ImageType", "Laterality", "Modality", "PatientName", "SeriesDescription", 
"ImageLaterality", "ViewPosition", "ProtocolName", 'Study DateTime',
'DNA_available', 'time_since_exam', 'case', 'date_dx', 'dx_time_after_exam', 'dx_years_after_exam', 'date_of_last_contact', 
'follow_up_time', 'follow_up_years', 'path',]

df_s2 = df[columns_to_keep]


#create png_paths(file_path) as column
png_dir= "/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/chimec_mammo_pngs/"
df_s2['file_path'] = df_s2['path'].replace(to_replace= "/gpfs/data/huo-lab/Image/ChiMEC/", value= png_dir, regex=True)
df_s2['file_path'] = df_s2['file_path'] + '.png'
df_s2['file_path'].head(20)

#ImageLaterality most times stated but sometimes NaN (and Laterality available in other column). Laterality on the other hand most times missing
#ViewPosition almost always stated and when missing ProtocolName also missing most times?
#Mirai takes "FOR PRESENTATION"    
#SeriesDescription not uniform
#PatientID==study_id
#AccessionNumber ==exam_id
#Screen out patients with breast implants


# %%
df_s2.head(30)

# %% [markdown]
# check if path exists
# remove burned-in annotations

# %%
df_s3 =df_s2[df_s2["BurnedInAnnotation"]== 'NO'] #filtering based on annotations
df_s3 = df_s3[df_s3["BreastImplantPresent"]=='NO'] # filter out patients with breast implants

columns_to_keep2= ["study_id", "exam_id",
"ImageLaterality", "ViewPosition", 'case','dx_years_after_exam','follow_up_years', 'file_path',]
df_s3= df_s3[columns_to_keep2]

key_dict= {"study_id": "patient_id", "exam_id": "exam_id", "ImageLaterality": "laterality", "ViewPosition": "view", "file_path":"file_path", "dx_years_after_exam": 
"years_to_cancer" , "follow_up_years": "years_to_last_followup"}

df_s3= df_s3.rename(mapper= key_dict, axis=1)
df_s3.head(20)


# %%
#check if png file exists
df_s3["path_exists"] = df_s3["file_path"].apply(lambda x: os.path.isfile(x))


# %%
df_s3.loc[df_s3["case"]== False, 'years_to_cancer'] = 100
#set years_to_cancer as 100 for controls
df_s3['years_to_cancer']=df_s3['years_to_cancer'].astype(int)

df_s3.loc[df_s3['case']==True, 'years_to_last_followup']= df_s3['years_to_cancer']
# replace year_to_last_followup from NaN to years_to_cancer for cases
# & set year_to_cancer as integer type data 
df_s3['years_to_last_followup']=df_s3['years_to_last_followup'].astype(int)
 
df_s3['split_group']= 'test' #making split_group column for validation 
df_s3.head(20)

len(df_s3)


# %%
#create interim csv
df_s3.to_csv(path_or_buf= "/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/interim_full_mammo_metadata.csv")

# %%
