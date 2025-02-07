# Because of the design of the algorithm, this calculation process is very slow
# This code is used to evaluate the conservation on each lncRNA genome using a sliding window algorithm
# In the future we will update an efficient algorithm based on convolutional operations
# If you want to implement a similar analysis for sliding window calculations, use that algorithmic implementation 
# (will update the address in this repository when released)

# Note that the implementation of this analysis is based on our server
# If you want to reproduce the code, change the address of the input data 
# the corresponding data is in the data folder



# import package
from Conservertion_Tools import calc
import pandas as pd
import pyBigWig


# import file
# PhastCon-based Conservativeness Score Downloaded from UCSC Genome Browser
all_trans_bed = pd.read_csv(
    '/data/jxwang_data/WMDS_lncRNA/bed_file/all.transcripts.bed', sep='\t', header=None)
Con_bw_file = pyBigWig.open(
    '/data/jxwang_data/WMDS_lncRNA/Genome_ref_file/hg38.phastCons100way.bw')
# WMDS driver result
all_LNC_df = pd.read_excel(
    '/data/jxwang_data/WMDS_lncRNA/data_input/WMDS_data_latest/LNC_name.xlsx', header=None)
# FDR 0.01
all_14_driver = pd.read_excel(
    '/data/jxwang_data/WMDS_lncRNA/data_input/WMDS_data_latest/fdr0.01_logfc0.5/gene_14.xlsx', header=None)

all_cancer_DriverList = []
for cancer in range(all_14_driver.shape[1]):
    tmp = [x for x in list(all_14_driver.iloc[:, cancer]) if pd.isna(x) != True]
    all_cancer_DriverList += tmp
all_cancer_DriverSet = set(all_cancer_DriverList)

count_dict = {}
for gene in all_cancer_DriverSet:
    count_dict[gene] = all_cancer_DriverList.count(gene)
    

def pan_cancer_driveranalysis(cancer_number=-1, time_record=True, path='/data/jxwang_data/WMDS_lncRNA/data_output/PP_result/'):

    target_lst = []
    for _, gene in enumerate(count_dict):
        if int(count_dict[gene]) >= cancer_number:
            target_lst.append(gene)

    lst = target_lst
    
    if cancer_number == -1:
        print('calculate all lncRNAs')
        lst = list(all_LNC_df.iloc[:, 0])

    # Conservation Analysis

    calc_1 = calc(all_trans_bed, Con_bw_file, P_bw_file, lst)
    # trans gene symbol into ENST_ID
    ENST_ID = calc_1.get_trans_ID_from_gene()
    print('Trans number: {}'.format(len(ENST_ID)))


    # phastCon
    best = calc_1.best_200_phastCon()


    # Save result
    if cancer_number != -1:

        filename3 = str(path) + 'best_more' + str(cancer_number) + '.csv'

        best.to_csv(filename3, index=True)
    else:
        filename3 = str(path) + 'best_more' + '_all.csv'

        best.to_csv(filename3, index=True)
    
    print(calc_1.time_recordor())


# all drivers
pan_cancer_driveranalysis(1, path='/data/jxwang_data/WMDS_lncRNA/data_output/PP_result/fdr001_logfc0.5/')
# pan-cancer drivers
pan_cancer_driveranalysis(2, path='/data/jxwang_data/WMDS_lncRNA/data_output/PP_result/fdr001_logfc0.5/')
pan_cancer_driveranalysis(5, path='/data/jxwang_data/WMDS_lncRNA/data_output/PP_result/fdr001_logfc0.5/')
pan_cancer_driveranalysis(7, path='/data/jxwang_data/WMDS_lncRNA/data_output/PP_result/fdr001_logfc0.5/')
pan_cancer_driveranalysis(10, path='/data/jxwang_data/WMDS_lncRNA/data_output/PP_result/fdr001_logfc0.5/')
# all lncRNAs
pan_cancer_driveranalysis(path='/data/jxwang_data/WMDS_lncRNA/data_output/PP_result/fdr001_logfc0.5/')