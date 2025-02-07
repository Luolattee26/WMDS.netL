# This code is used to preprocess the csv data for input to the training code for model training
# Please note that since the training process takes place on our servers, 
# if you want to reproduce the code, be sure to change the path to the file beforehand!



from utils import *


pancancer_dataset_prepare(read_path='/data/jxwang_data/WMDS_lncRNA/data_input/ml_model/combine_GTEx/useTopFeatures',
                          save_path='/data/jxwang_data/WMDS_lncRNA/data_input/ml_model/combine_GTEx/useTopFeatures',
                          auto_sub_str=(0, 4),
                          extension='.csv',
                          random_state=12)