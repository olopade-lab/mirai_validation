# mirai_validation
Repository for validating Mirai, a DL mammography based breast cancer risk prediction tool, in the ChiMEC cohort.
`chimec_data_describe.ipynb` looks at the image metadata generated from preprocessing, describes the data, and splits it into groups for subset analysis.
`results_explore.ipynb` merges input sets (e.g. full set, black cohort etc) with the respective outputted predictions, and calls on functions 
in utils/metadata_inspect.py to compute AUCs and c-indices with their confidence intervals. 
