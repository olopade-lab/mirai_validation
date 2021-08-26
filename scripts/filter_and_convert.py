# Modified from the original version by Adam Yala
import glob
import os

import pandas as pd
from p_tqdm import p_umap

from oncodata.dicom_metadata.get_dicom_metadata import get_dicom_metadata
from oncodata.dicom_to_png.dicom_to_png import dicom_to_png_dcmtk

dicom_dir = "/gpfs/data/huo-lab/Image/ChiMEC"
png_dir = "/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/chimec_mammo_pngs"
series_metadata_path = "/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/chimec_mammo_retry.csv"


def dicom_path_to_png_path(dicom_path, dicom_dir, png_dir, dicom_ext):
    """Converts a DICOM path to a PNG path by replacing dicom_dir with png_dir in the path.

    Arguments:
        dicom_path(str): Path to a DICOM file.
        dicom_dir(str): Path to a directory containing DICOM files.
        png_dir(str): Path to a directory where PNG version of the
            DICOM images will be saved.
        dicom_ext(str): The extension of the dicom files (should include the dot)

    Returns:
        The same path as dicom_path but with dicom_dir replaced by
        png_dir and with dicom_ext replaced with '.png'.
    """

    dicom_path_after_dir = dicom_path.replace(dicom_dir, "", 1).strip("/")
    if dicom_ext != "":
        png_path_after_dir = dicom_path_after_dir.replace(dicom_ext, "")
    png_path_after_dir = dicom_path_after_dir + ".png"
    png_path = os.path.join(png_dir, png_path_after_dir)

    return png_path


metadata = pd.read_pickle(
    "/gpfs/data/huo-lab/Image/annawoodard/maicara/data/interim/downloaded_imaging_metadata.pkl"
)
metadata = metadata[~pd.isnull(metadata.exam_id)]
metadata = metadata[metadata["Modality"].str.contains("MG")]

print("loaded {} cases".format(metadata[metadata.case == True].study_id.nunique()))
print("loaded {} controls".format(metadata[metadata.case == False].study_id.nunique()))
print(f"will filter {len(metadata)} exams for presentation type and view position")


def get_series_metadata(dicom_dir, exam):
    dicom_files = glob.glob(
        os.path.join(dicom_dir, str(exam["study_id"]), exam["exam_id"], "*/*dcm")
    )
    series_metadata = []
    for f in dicom_files:
        dicom_metadata = {}
        try:
            dicom_metadata = get_dicom_metadata(f)
        except Exception as e:
            print(f"problem getting dicom metadata: {e}")
        if (not "ViewPosition" in dicom_metadata) or (
            not "PresentationIntentType" in dicom_metadata
        ):
            continue
        passes_position_cut = dicom_metadata["ViewPosition"] in ["CC", "MLO"]
        passes_presentation_cut = (
            dicom_metadata["PresentationIntentType"] == "FOR PRESENTATION"
        )
        if passes_position_cut and passes_presentation_cut:
            dicom_metadata.update(exam)
            dicom_metadata["dicom_path"] = f
            dicom_metadata["png_path"] = f.replace("/gpfs/data/huo-lab/Image/ChiMEC/", "/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/chimec_mammo_pngs/", ) + '.png'
            series_metadata.extend([dicom_metadata])
    return series_metadata


series_metadata = p_umap(
    get_series_metadata,
    [dicom_dir for _ in range(len(metadata))],
    metadata.to_dict("records"),
    num_cpus=2
)
series_metadata = sum(series_metadata, [])

df = pd.DataFrame(series_metadata)
df.to_csv(series_metadata_path)

# dicom_paths = df.path.tolist()
# png_paths = [dicom_path_to_png_path(p, dicom_dir, png_dir, ".dcm") for p in dicom_paths]
# print(f"Converting {len(dicom_paths)} dicoms to PNGs")
# p_umap(dicom_to_png_dcmtk, dicom_paths, png_paths, [[] for _ in dicom_paths], num_cpus=2)
