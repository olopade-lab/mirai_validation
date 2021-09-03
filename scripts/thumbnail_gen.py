import glob 
import os
import sys
import pandas as pd
from p_tqdm import p_umap
from PIL import Image 
print ("all libs imported")
size = (120, 120)

png_dir = "/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/chimec_mammo_pngs"
thumb_dir = "/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/chimec_mammo_pngs/thumbnails"
series_metadata_path = "/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/data/chimec_mammo_retry.csv"

def png_path_to_thumbnail_path (png_path, png_dir, thumb_dir):
    """Converts a DICOM path to a PNG path by replacing dicom_dir with png_dir in the path.

    Arguments:
        png_dir(str): Path to a directory containing PNG images
        thumb_dir(str): Path to the directory where thumbnail versions of PNGS
        will be stored

    Returns:
        The same path as png_path but with a subfolder for thumbnails
    """
    png_path_after_dir = png_path.replace(png_dir, "", 1).strip("/")
    thumb_path = os.path.join(thumb_dir, png_path_after_dir)
    return thumb_path 

metadata = pd.read_csv(series_metadata_path)
df = pd.DataFrame(metadata)
print ("checking if png paths exist")
total_n = len(df)
df["true"]= df.png_path.apply(lambda x: os.path.isfile(x))
df = df[df.true == True]

print ("of {} png paths, {} exist and will be converted to thumbnails".format(
    total_n, len(df)
))

png_paths= df.png_path.tolist()
thumb_paths = [png_path_to_thumbnail_path(p, png_dir, thumb_dir) for p in png_paths]
print(f"Converting {len(png_paths)} PNGs to thumbnails")

for image in png_paths:
    try:
        im = Image.open(image)
    except IOError as e:
        #report error and then skip to the next argument
        print("Problem opening", image, ":", e)
        continue
    im.resize(size)
    path = png_path_to_thumbnail_path(image, png_dir, thumb_dir)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    im.save(path)