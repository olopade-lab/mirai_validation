import pandas as pd
import numpy as np
from scipy import misc
import matplotlib.pyplot as plt
import seaborn as sns 
from PIL import Image
from io import BytesIO
from IPython.display import HTML
import re
import os
import sys
from p_tqdm import p_umap

series_metadata_path= "/gpfs/data/huo-lab/Image/annawoodard/maicara/data/interim/mammo_v9/series_metadata.pkl"
data_path = "/gpfs/data/phs/groups/Projects/Huo_projects/SPORE/ojomoleye/data/CRDW_Registry_IndexDx_2020_Jul23.csv"
data = pd.read_csv(data_path)

mrn_to_study_id = pd.read_csv(
    "/gpfs/data/phs/groups/Projects/Huo_projects/SPORE/ojomoleye/data/mrn_to_study_id.csv",
    names=["mrn", "study_id"],
)

downloaded_images = pd.read_pickle(
    "/gpfs/data/phs/groups/Projects/Huo_projects/SPORE/ojomoleye/data/downloaded_imaging_metadata.pkl"
)

cases_and_controls = pd.read_csv(
    "/gpfs/data/phs/groups/Projects/Huo_projects/SPORE/ojomoleye/data/dr_7934_pats.txt",
    sep="|",
)
spore_registration = pd.read_csv(
    "/gpfs/data/phs/groups/Projects/Huo_projects/SPORE/ojomoleye/data/SPORERegistrationDat_DATA_2021-06-21_0927.csv"
)

race_data= pd.read_excel("/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/data/MRN_to_Race.xls")

mrn_to_tumortype = pd.read_excel("/gpfs/data/phs/groups/Projects/Huo_projects/SPORE/ojomoleye/data/mrn_tumortype.xlsx")



def print_summary(metadata):
    '''
    prints a summary of metadata table, ouputs: total number of files(pngs);
    total number of cases and case exams; total number of controls and control exams
    '''
    total_n = len(metadata) #total number of dicom files
    total_p = metadata.study_id.nunique() #total number of patients ##put try and except for if names replaced (patient_id)
    total_e = metadata.exam_id.nunique() #total number of exams
    total_CE = metadata[metadata.case==True].exam_id.nunique() #total number of case exams
    total_NE = metadata[metadata.case==False].exam_id.nunique() # total number of control exams
    total_cases = metadata[metadata.case==True].study_id.nunique() # total number of cases
    total_controls = metadata[metadata.case==False].study_id.nunique() # total number of controls
    return print("{}AVAILABLE/LOADED{} \n {} patients \n {} case exams for {} cases \n {} control exams for {} controls \n {} total entries".format(
        "*"*5, "*"*5,
        total_p,
        total_CE, total_cases,
        total_NE, total_controls,
        total_n))

def print_summary2(metadata):
    '''
    prints a summary of metadata table, ouputs: total number of files(pngs);
    total number of cases and case exams; total number of controls and control exams
    vers2: study_id replaced with patient_id
    '''
    total_n = len(metadata) #total number of dicom files
    total_p = metadata.patient_id.nunique() #total number of patients ##put try and except for if names replaced (patient_id)
    total_e = metadata.exam_id.nunique() #total number of exams
    total_CE = metadata[metadata.case==True].exam_id.nunique() #total number of case exams
    total_NE = metadata[metadata.case==False].exam_id.nunique() # total number of control exams
    total_cases = metadata[metadata.case==True].patient_id.nunique() # total number of cases
    total_controls = metadata[metadata.case==False].patient_id.nunique() # total number of controls
    return print("{}AVAILABLE/LOADED{} \n {} patients \n {} case exams for {} cases \n {} control exams for {} controls \n {} total entries".format(
        "*"*5, "*"*5,
        total_p,
        total_CE, total_cases,
        total_NE, total_controls,
        total_n))


def case_control_prop(metadata):
    ''' returns the case and control proportions'''
    case_prop = \
        (metadata[metadata.case==True].study_id.nunique() /metadata.study_id.nunique()) * 100
    control_prop = \
        (metadata[metadata.case==False].study_id.nunique() /metadata.study_id.nunique()) * 100
    return print(f'there are {case_prop :.2f}% cases and {control_prop :.2f}% controls')

def exist_filter(metadata):
    """
    takes in dataframe with png paths as argument
    filters out rows with non-exitent corresponding pngs
    returns filtered dataframe
    """
    total_n = len(metadata) #total number of dicom files
    total_CE = metadata[metadata.case==True].exam_id.nunique() #total number of case exams
    total_NE = metadata[metadata.case==False].exam_id.nunique() # total number of control exams
    total_cases = metadata[metadata.case==True].study_id.nunique() # total number of cases
    total_controls = metadata[metadata.case==False].study_id.nunique() # total number of controls
    metadata["png_true"]= metadata.png_path.apply(lambda x: os.path.isfile(x))
    metadata = metadata[metadata.png_true == True]
    print(
        "will filter out {} files, {} case exams for {} unique cases and {} control exams for {} unique controls".format(
        total_n - len(metadata), total_CE - metadata[metadata.case==True].exam_id.nunique(), 
        total_cases - metadata[metadata.case==True].study_id.nunique(),
        total_NE - metadata[metadata.case==False].exam_id.nunique(),
        total_controls - metadata[metadata.case==False].study_id.nunique())
    )
    return metadata 

def prep_full_set(metadata_path=series_metadata_path):
    import os, sys 
    ''' 
    this function takes in the metadata pkl file and prepares for Mirai input (full_set)
    full_set can be split into subsets by donwnstream functions
    metadata contains more columns and useful for description in the notebook chimec_data_describe.ipynb
    argument passed in should be the path to the metadata generated from preprocessing
    '''
    series_metadata= pd.read_pickle(metadata_path)
    metadata = series_metadata[series_metadata['png_path'].notnull()] #drop non-existent png_paths
    metadata = metadata[pd.isnull(metadata['filter'])] # drop if filter tag is present
    metadata['true_path'] = metadata['png_path'].apply(lambda x: os.path.exists(x))
    metadata = metadata[metadata['true_path']==True]
    columns_to_keep= ["study_id", "unique_exam_id",
    "ImageLaterality", "ViewPosition", 'case','years_from_exam_to_diagnosis','follow_up_years', 'png_path'] #keep only columns useful for Mirai input
    df = metadata[columns_to_keep]
    key_dict= {"study_id": "patient_id", "unique_exam_id": "exam_id", "ImageLaterality": "laterality", "ViewPosition": "view", "png_path":"file_path", 'years_from_exam_to_diagnosis': 
    "years_to_cancer" , "follow_up_years": "years_to_last_followup"}
    df= df.rename(mapper= key_dict, axis=1)
    df.loc[df["case"]== False, 'years_to_cancer'] = 100 #set years_to_cancer as 100 for controls
    df['years_to_cancer']= df['years_to_cancer'].astype(int) #convert years from float to integer
    # https://stackoverflow.com/questions/49161120/pandas-python-set-value-of-one-column-based-on-value-in-another-column
    df['years_to_last_followup'] = np.where(df.case == True,df.years_to_cancer,df.years_to_last_followup)
    # where it is a case, set value of years to follow up as years to cancer
    # otherwise, leave to retain current value
    # & set year_to_cancer as integer type data 
    df['years_to_last_followup'] = df['years_to_last_followup'].astype(int)
    df['split_group']= 'test' #making split_group column for validation 
    return df, metadata 

def prep_by_race(df):
    ''' where df is the dataframe of the full image set,
    this function splits into race groups based on available data
    and returns three separate dataframes for white, black, hispanic and asian
    '''
    # study id and mrn both integer type in mrn_to_study_id
    df["study_id"] = df["patient_id"]
    df = df.merge(mrn_to_study_id, how="left", on="study_id", indicator=True)
    df._merge.value_counts()  # all patients in CRDW represented in mrn_to_studyid map #all merged
    df.drop(columns="_merge", inplace=True) #drop merge validation column
    race_data.columns = race_data.columns.map(str.lower)
    race_data.columns = race_data.columns.map(str.strip) #merge did not work before stripping (white space?)
    df = df.merge(race_data, how="left", on="mrn", indicator=True)
    df._merge.value_counts()  # all patients in CRDW represented in mrn_to_studyid map #all merged
    df.race.value_counts(dropna=False)
    df = df[df._merge =="both"]
    df.drop(columns="_merge", inplace=True)
    print(
    "{} White - {} cases and {} controls \n{} Black - {} cases and {} controls \n{} Hispanic - {} cases and {} controls \n{} Asian - {} cases and {} controls \n{} Native Ame - {} cases and {} controls ".format(
        df[df['race'] == "White"].patient_id.nunique(), df[(df['race'] == "White") & \
        (df['years_to_cancer']!=100)].patient_id.nunique(), df[(df['race'] == "White") & \
        (df['years_to_cancer']==100)].patient_id.nunique(), df[df['race'] == "Black"].patient_id.nunique(), \
        df[(df['race'] == "Black") & (df['years_to_cancer']!=100)].patient_id.nunique(), \
        df[(df['race'] == "Black") & (df['years_to_cancer']==100)].patient_id.nunique(), \
        df[df['race'] == "Hispanic"].patient_id.nunique(), df[(df['race'] == "Hispanic") & \
        (df['years_to_cancer']!=100)].patient_id.nunique(), df[(df['race'] == "Hispanic") & \
        (df['years_to_cancer']==100)].patient_id.nunique(),df[df['race'] == "Asian"].patient_id.nunique(), \
        df[(df['race'] == "Asian") & (df['years_to_cancer']!=100)].patient_id.nunique(), \
        df[(df['race'] == "Asian") & (df['years_to_cancer']==100)].patient_id.nunique(),
        df[df['race'] == 'Native American'].patient_id.nunique(), df[(df['race'] == "Native American") & \
        (df['years_to_cancer']!=100)].patient_id.nunique(), df[(df['race'] == "Native American") & \
        (df['years_to_cancer']==100)].patient_id.nunique()
                ))
    #create pie chart showing race distribution
    plt.rcdefaults()
    ax= df['race'].value_counts().plot(kind='pie', autopct="%1.1f%%")
    ax.set_title("Race Distribution in ChiMEC")
    plt.rcParams.update({'font.size':5})
    plt.show()
    white = df[df.race =="White"]
    black = df[df.race =="Black"]
    hispanic = df[df.race =="Hispanic"]
    asian_and_native = df[(df.race =="Asian") | (df.race=="Native American")]
    white.drop(columns=['study_id', 'mrn', 'race'], axis=1, inplace=True)
    black.drop(columns=['study_id', 'mrn', 'race'], axis=1, inplace=True)
    hispanic.drop(columns=['study_id', 'mrn', 'race'], axis=1, inplace=True)
    asian_and_native.drop(columns=['study_id', 'mrn', 'race'], axis=1, inplace=True)
    return white, black, hispanic, asian_and_native 

def prep_by_tumor_subtype(df):
    ''' 
    takes in the dataframe of the full set and outputs dataframe split by tumor 
    subtype
    '''
    df["study_id"] = df["patient_id"]
    df = df.merge(mrn_to_study_id, how="left", on="study_id", indicator=True) 
    df.drop(columns="_merge", inplace=True) 
    mrn_to_tumortype.columns = mrn_to_tumortype.columns.map(str.lower)
    mrn_to_tumortype.columns = mrn_to_tumortype.columns.map(str.strip)
    df = df.merge(mrn_to_tumortype, how="left", on="mrn", indicator=True)
    cases= df.loc[df['years_to_cancer']!=100]
    cases_w_type = cases.loc[cases['_merge']=='both'] 
    unique_case_exams_w_type= cases_w_type['exam_id'].nunique()
    unique_case_exams = cases['exam_id'].nunique()
    percent_with_type= (unique_case_exams_w_type / unique_case_exams) * 100 
    print(
    'there are {} unique case exams \n{} unique case exams with tumor type\
         \n{} percent of case exams with tumor type available'.format(
             unique_case_exams, unique_case_exams_w_type,percent_with_type)
            )   
    df2 = cases_w_type 
    df2['HR_positive'] = None 
    df2['HER2_positive'] = None 
    df2.loc[(df2['er_status1']=='Pos') | (df2['pr_status1']=='Pos'),'HR_positive'] = True
    df2.loc[(df2['er_status1']=='Neg') & (df2['pr_status1']=='Neg'),'HR_positive'] = False
    df2.loc[df2['her2_status1']=='Pos', 'HER2_positive'] = True
    df2.loc[df2['her2_status1']=='Neg', 'HER2_positive'] = False
    columns_to_keep=['patient_id', 'exam_id', 'laterality', 'view', 'case', 'years_to_cancer', 'years_to_last_followup', 'file_path', 'split_group']
    hr_pos = df2[df2['HR_positive'] == True][columns_to_keep]
    hr_neg = df2[df2['HR_positive'] == False][columns_to_keep]
    her2_pos = df2[df2['HER2_positive'] == True][columns_to_keep]
    her2_neg = df2[df2['HER2_positive'] == False][columns_to_keep]
    triple_neg = df2[(df2['HR_positive']==False) & (df2['HER2_positive']==False)][columns_to_keep]
    return hr_pos, hr_neg, her2_pos, her2_neg, triple_neg


def prep_results_input(results_path, input_path):
    ''' 
    takes in mirai prediction path from results_path, and the input csv path
    from input_path
    outputs dataframe to be used in survival analysis to combined_path
    '''
    input = pd.read_csv(input_path)
    results = pd.read_csv(results_path)
    input['study_id'] = input['patient_id']
    results['exam_id']= results['patient_exam_id'].str.split("\t", n=1, expand = True)[1] #split out exam_id to merge with input sheet
    results = results.merge(input, how='left', on='exam_id', indicator=True)
    columns_to_keep= ["exam_id", "1_year_risk", "2_year_risk", "3_year_risk", "4_year_risk", "5_year_risk", "case", "years_to_cancer", "years_to_last_followup", "patient_id"]
    results = results[columns_to_keep]
    results["study_id"] = results["patient_id"]
    print("there are {} duplicated rows and {} unduplicated rows \n{} total rows".\
    format(results.duplicated().sum(), (~results.duplicated()).sum(), len(results)))
    print('*'*20)
    print("dropping duplicate rows; so only one entry per exam")
    results = results.drop_duplicates() 
    return results 


def auc_calculator(df, bootstraps):
    ''' 
    calculates AUC for 1-5 years,
    also outputs case & total exams per year
    '''
    from sklearn import metrics
    import scipy.stats
    # https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_curve.html
    import numpy as np
    df['event'] = True #seetting all events to true at first
    df.loc[df['years_to_cancer'] == 100, 'event'] = False #years_to_cancer for controls is 100; set event to false for controls
    df['time'] = df['years_to_cancer'] 
    # time is years to cancer for cases and years to last follow up (censoring time) for controls
    df.loc[df['years_to_cancer'] == 100, 'time'] = df["years_to_last_followup"]
    #df.drop(columns= ['years_to_cancer', 'years_to_last_followup', 'patient_id', 'study_id' ], inplace = True)
    df['case_1yr'] = False
    df['case_2yr'] = False
    df['case_3yr'] = False
    df['case_4yr'] = False
    df['case_5yr'] = False
    df.loc[(df['case'] == True) & (df['time'] <= 1), 'case_1yr'] = True
    df.loc[(df['case'] == True) & (df['time'] <= 2), 'case_2yr'] = True
    df.loc[(df['case'] == True) & (df['time'] <= 3), 'case_3yr'] = True
    df.loc[(df['case'] == True) & (df['time'] <= 4), 'case_4yr'] = True
    df.loc[(df['case'] == True) & (df['time'] <= 5), 'case_5yr'] = True
    
    case_exams_at_1yr =  len(df[df['case_1yr'] == True])
    case_exams_at_2yr = len(df[df['case_2yr'] == True])
    case_exams_at_3yr = len(df[df['case_3yr'] == True])
    case_exams_at_4yr = len(df[df['case_4yr'] == True])
    case_exams_at_5yr = len(df[df['case_5yr'] == True])

    def get_column(df, years, column):
        return np.array(df[(df['case'] == True) | (df["years_to_last_followup"] >= years)][column])
    
    def mean_confidence_interval(data, confidence=0.95):
        a = 1.0 * np.array(data)
        n = len(a)
        m, se = np.mean(a), scipy.stats.sem(a)
        h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
        return m, m-h, m+h
    # https://stackoverflow.com/questions/15033511/compute-a-confidence-interval-from-sample-data
    
    y = get_column(df, 1, 'case_1yr')
    pred = get_column(df, 1, '1_year_risk')
    exam_ids = get_column(df, 1, 'exam_id')
    fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=True)
    print('1 year auc is {} \n{} case exams at 1 year, out of {} total exams'.format(
        metrics.auc(fpr, tpr), case_exams_at_1yr, len(set(exam_ids))))
    num_exams = len(y)
    sample_size = int(len(set(exam_ids)) * 0.15)
    aucs1yr = []
    for _ in range(bootstraps):
        rows = np.random.randint(num_exams, size=sample_size)
        if len(np.unique(y[rows]) > 1):
            fpr, tpr, thresholds = metrics.roc_curve(y[rows], pred[rows], pos_label=True)
            aucs1yr.append(metrics.auc(fpr,tpr))

    x= mean_confidence_interval(data=aucs1yr, confidence=0.95)
    print(f'mean auc: {x[0]}; confidence interval: {x[1:3]}')
    print('*'*30)

    y = get_column(df, 2, 'case_2yr')
    pred = get_column(df, 2, '2_year_risk')
    exam_ids = get_column(df, 2, 'exam_id')
    fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=True)
    print('2 year auc is {} \n{} case exams at 2 year, out of {} total exams'.format(
        metrics.auc(fpr, tpr), case_exams_at_2yr, len(set(exam_ids))))
    num_exams = len(y)
    sample_size = int(len(set(exam_ids)) * 0.15)
    aucs2yr = []
    for _ in range(bootstraps):
        rows = np.random.randint(num_exams, size=sample_size)
        if len(np.unique(y[rows]) > 1):
            fpr, tpr, thresholds = metrics.roc_curve(y[rows], pred[rows], pos_label=True)
            aucs2yr.append(metrics.auc(fpr,tpr))
        
    x= mean_confidence_interval(data=aucs2yr, confidence=0.95)
    print(f'mean auc: {x[0]}; confidence interval: {x[1:3]}')
    print('*'*30)

    y = get_column(df, 3, 'case_3yr')
    pred = get_column(df, 3, '3_year_risk')
    exam_ids = get_column(df, 3, 'exam_id')
    fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=True)
    print('3 year auc is {} \n{} case exams at 3 year, out of {} total exams'.format(
        metrics.auc(fpr, tpr), case_exams_at_3yr, len(set(exam_ids))))
    num_exams = len(y)
    sample_size = int(len(set(exam_ids)) * 0.15)
    aucs3yr = []
    for _ in range(bootstraps):
        rows = np.random.randint(num_exams, size=sample_size)
        if len(np.unique(y[rows]) > 1):
            fpr, tpr, thresholds = metrics.roc_curve(y[rows], pred[rows], pos_label=True)
            aucs3yr.append(metrics.auc(fpr,tpr))
    x= mean_confidence_interval(data=aucs3yr, confidence=0.95)
    print(f'mean auc: {x[0]}; confidence interval: {x[1:3]}')
    print('*'*30)

    y = get_column(df, 4, 'case_5yr')
    pred = get_column(df, 4, '3_year_risk')
    exam_ids = get_column(df, 4, 'exam_id')
    fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=True)
    print('4 year auc is {} \n{} case exams at 4 year, out of {} total exams'.format(
        metrics.auc(fpr, tpr), case_exams_at_4yr, len(set(exam_ids))))
    num_exams = len(y)
    sample_size = int(len(set(exam_ids)) * 0.15)
    aucs4yr = []
    for _ in range(bootstraps):
        rows = np.random.randint(num_exams, size=sample_size)
        if len(np.unique(y[rows]) > 1):
            fpr, tpr, thresholds = metrics.roc_curve(y[rows], pred[rows], pos_label=True)
            aucs4yr.append(metrics.auc(fpr,tpr))
    x= mean_confidence_interval(data=aucs4yr, confidence=0.95)
    print(f'mean auc: {x[0]}; confidence interval: {x[1:3]}')
    print('*'*30)

    y = get_column(df, 5, 'case_5yr')
    pred = get_column(df, 5, '5_year_risk')
    exam_ids = get_column(df, 5, 'exam_id')
    fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=True)
    print('5 year auc is {} \n{} case exams at 5 year, out of {} total exams'.format(
        metrics.auc(fpr, tpr), case_exams_at_5yr, len(set(exam_ids))))
    num_exams = len(y)
    sample_size = int(len(set(exam_ids)) * 0.15)
    aucs5yr = []
    for _ in range(bootstraps):
        rows = np.random.randint(num_exams, size=sample_size)
        if len(np.unique(y[rows]) > 1):
            fpr, tpr, thresholds = metrics.roc_curve(y[rows], pred[rows], pos_label=True)
            aucs5yr.append(metrics.auc(fpr,tpr))
    x= mean_confidence_interval(data=aucs4yr, confidence=0.95)
    print(f'mean auc: {x[0]}; confidence interval: {x[1:3]}')
    print('*'*30)
    return print('AUC values computed for df')

def compute_concordance_index(df,bootstraps):
    ''' 
    computes the concordance index for the full set
    returns the concordance index per year, and an average of c_indices across years
    returns concordance index using average risk as well as the mean and 
    confidence interval of this c-index
    '''
    from sksurv.metrics import concordance_index_censored as cic
    import scipy.stats
    df['average_risk']= df[['1_year_risk', '2_year_risk', '3_year_risk', \
    '4_year_risk', '5_year_risk']].mean(axis=1)
    df['event'] = True #seetting all events to true at first
    df.loc[df['years_to_cancer'] == 100, 'event'] = False #years_to_cancer for controls is 100; set event to false for controls
    df['time'] = df['years_to_cancer'] 
    # time is years to cancer for cases and years to last follow up (censoring time) for controls
    df.loc[df['years_to_cancer'] == 100, 'time'] = df["years_to_last_followup"]
    ci_dict = dict()
    event_indicator = np.array(df['event'])
    event_time= np.array(df['time'])
    risk = np.array(df['average_risk']) 
    c_index_full = cic(event_indicator, event_time, estimate= risk)[0]
    c_index_list = []

    def mean_confidence_interval(data, confidence=0.95):
        a = 1.0 * np.array(data)
        n = len(a)
        m, se = np.mean(a), scipy.stats.sem(a)
        h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
        return m, m-h, m+h

    for _ in range(bootstraps):
        sample = df.sample(frac=0.15, replace=False)
        event_indicator = np.array(sample['case'])
        event_time= np.array(sample['time'])
        risk = np.array(sample['average_risk']) 
        c_index = cic(event_indicator, event_time, estimate=risk)[0]
        c_index_list.append(c_index)
    ci_dict['c-index full'] = c_index_full 
    ci_dict['c-index mean'] = mean_confidence_interval(data=c_index_list, confidence=0.95)[0]
    ci_dict['95% CI'] = mean_confidence_interval(data=c_index_list, confidence=0.95)[1:3]
    return print(ci_dict)


def get_old_paths(df):
    '''
    takes in a parepared Mirai_input dataframe (from mammoV8)
    replaces file_path with older 'DCMTK' converted image paths
    Useful for debugging/comparing
    '''
    import glob 
    path_now = '/gpfs/data/huo-lab/Image/annawoodard/maicara/data/interim/mammo_v8/png/'
    path_update = '/gpfs/data/huo-lab/Image/ojomoleye/projects/mirai_validation/chimec_mammo_pngs/'
    df['file_path']= df['file_path'].str.replace(pat=path_now, repl=path_update)
    df['file_path']= df['file_path'].str.replace(pat='_L_MLO', repl='')
    df['file_path']= df['file_path'].str.replace(pat='_L_CC', repl='')
    df['file_path']= df['file_path'].str.replace(pat='_R_CC', repl='')
    df['file_path']= df['file_path'].str.replace(pat='_R_MLO', repl='')
    df['file_path1'] = df['file_path'].str.rsplit('/', n=1, expand=True)[0]
    df['file_path2']= df['file_path'].str.rsplit('/', n=1, expand=True)[1]
    df['file_path2']= df['file_path'].str.rsplit('_', n= 1, expand=True)[1]
    df['file_path2']= df['file_path2'].str.replace(pat='.png', repl='')
    df['path']= df['file_path1'] + '/' + df['file_path2'] + '/*'
    df['png_path']= df.path.apply(lambda x: glob.glob(x))
    print(f'df length before exploding is {len(df)}')
    df= df.explode('png_path') #glob returns a list; explode to bring out the string
    print(f'df length after exploding is {len(df)}')
    df.drop(columns=['file_path1', 'file_path2', 'path', 'file_path'], inplace=True)
    df.rename(mapper={'png_path':'file_path'}, axis=1, inplace=True)
    df = df[df['file_path'].notnull()] #non existent paths- glob will not be able to retrive path
    return df 