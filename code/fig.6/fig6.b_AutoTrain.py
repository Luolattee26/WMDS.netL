# This code is used to automate the training of multiple machine learning models
# Please note that since the training process takes place on our servers, 
# if you want to reproduce the code, be sure to change the path to the file beforehand!



from utils import *
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import loguniform


# define model
clf_SVC = SVC(kernel='rbf', probability=True)
clf_LR = LogisticRegression(n_jobs=-1)
clf_RF = RandomForestClassifier(oob_score=False, random_state=26, n_jobs=-1)
# define hyperparams
params_SVC = {
    'C': loguniform.rvs(1e-1, 1e3, size=500),
    'gamma': loguniform.rvs(1e-4, 1e0, size=500)
}
params_LR = {
    'C': loguniform.rvs(1e-1, 1e3, size=500)
}
params_RF = {
    'n_estimators': np.arange(50, 1000, 10),
    'min_samples_leaf': np.arange(1, 20, 1),
    'max_features': np.random.uniform(0.3, 0.8, size=5000),
    'max_depth': np.arange(1, 50, 1),
    'min_samples_split': np.arange(3, 10, 1)
}

# auto training
auto_train(dataset_path='/data/jxwang_data/WMDS_lncRNA/data_input/ml_model/combine_GTEx/useTopFeatures',
           hyperparam_save_path='/data/jxwang_data/WMDS_lncRNA/data_output/hyperparam',
           hyperparam_list=[params_SVC, params_LR, params_RF],
           models=[clf_SVC, clf_LR, clf_RF],
           model_names=['SVC', 'LR', 'RF'],
           trainer_type='random',
           endswith='GTEx_useTopFeatures',
           dataset_extension='.pkl',
           auto_sub_str=(11, 15))