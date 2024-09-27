import pandas as pd
import numpy as np
from sklearn.feature_selection import VarianceThreshold, RFECV, SelectKBest, f_classif, mutual_info_classif
import matplotlib.pyplot as plt
from sklearn.preprocessing import RobustScaler, StandardScaler
from mlxtend.feature_selection import SequentialFeatureSelector as SFS
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
from sklearn.model_selection import StratifiedKFold
from tqdm import tqdm
from sklearn.model_selection import train_test_split
import pickle
import os
import sys
import re
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression


class DataSet:
    
    def __init__(self, path, scale = True, Scaler = 'Standard', test_size=0.3, random_state=26) -> None:
        self._expMat = pd.read_csv(path, sep=',', index_col=0).T
        self._X = self._expMat.iloc[:, 0:(self._expMat.shape[1]-1)]
        self._y = self._expMat.iloc[:, -1]
        self.scaled = False
        self.Scaler = Scaler
        self.scale = scale
        self._test_size = test_size
        self._random_state = random_state
    
    @property
    def X(self):
        if self.scale:
            if self.scaled == False:
                print('Performing scaling: {}'.format(self.Scaler))
                if self.Scaler == 'Standard':
                    self._X = StandardScaler().fit_transform(self._X)
                    self.scaled = True
                elif self.Scaler == 'Robust':
                    self._X = RobustScaler().fit_transform(self._X)
                    self.scaled = True
                else:
                    print('Select scaler form Standard and Robust \n And run it again')
                    self.scaled = False
                    return
            else:
                print('The exp data has been scaled')
        else:
            print('Choose not to scale')
        return self._X
            
    @property
    def y(self):
        return self._y
    
    @property
    def df_autogluon(self):
        return self._expMat
    
    def save(self, path, name):
        X_train, X_test, y_train, y_test = train_test_split(
                self.X, self.y, test_size=self._test_size, random_state=self._random_state, stratify=self.y)
        file_name = path + str('/data_split_{}.pkl'.format(name))
        with open(file_name, 'wb') as file:
            pickle.dump((X_train, X_test, y_train, y_test), file)
        return
    
def ROC_curve(hyperparams, X, y, model_names=[], full_data=False, models=None, test_size=0.3, random_state=26, trained=False, tittle=None,
              width=6, height=6, dpi=300, format='jpg', path=False, use_defult_hyperparams=False, show=False):
    from sklearn.metrics import roc_curve, auc, accuracy_score
    import matplotlib.pyplot as plt
    
    if full_data:
        stratify = y
        X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=test_size, random_state=random_state, stratify=stratify)
    else:
        if trained == False:
            raise ValueError
        X_test = X
        y_test = y

    plt.figure(figsize=(10, 8))
    models_figured = []
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    linestyles = ['-', '--', '-.', ':']
    for index in range(len(models)):
        model = models[index]
        if use_defult_hyperparams:
            model.set_params()
        else:
            model.set_params(**hyperparams[index])
        if trained == False:
            model.fit(X_train, y_train)
        models_figured.append(model)
    
    df = pd.DataFrame(columns=['AUC', 'Accuracy'], index=model_names)
    
    index = -1
    for model in models_figured:
        index += 1
        probs = model.predict_proba(X_test)
        # print(y_test)
        # print(probs[:, 1])
        fpr, tpr, _ = roc_curve(np.array(y_test), np.array(probs[:, 1]))
        roc_auc = auc(fpr, tpr)
        acc = accuracy_score(np.array(y_test), model.predict(X_test))
        # plt.plot(fpr, tpr, label='{}: AUC = {:.3f} ACC. = {:.3f}'.format(model_names[index], roc_auc, acc),
        #          color=colors[index % len(colors)], linestyle=linestyles[index % len(linestyles)], linewidth=5)
        plt.plot(fpr, tpr,
                 color=colors[index % len(colors)], linestyle=linestyles[index % len(linestyles)], linewidth=5)
        
        df.loc[model_names[index]] = [roc_auc, acc]
    
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([-0.05, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('FPR', fontsize=50)
    plt.ylabel('TPR', fontsize=50)
    if tittle == None:
        plt.title('ROC', fontsize=50)
    else:
        plt.title('{}'.format(tittle), fontsize=50)
    # plt.legend(loc="lower right", prop={'size': 25})
    
    df = df.map('{:.3f}'.format)
    table = plt.table(cellText=df.values,
                    colLabels=df.columns,
                    rowLabels=df.index,
                    cellLoc = 'center', 
                    bbox=[0.52, 0, 0.475, 0.4])
    table.scale(1, 1.5)
    table.auto_set_font_size(False)
    table.set_fontsize(25)
    table.auto_set_column_width(col=list(range(len(df.columns))))
    
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    for i in range(len(df.index)):
        for j in range(len(df.columns)):
            cell = table.get_celld()[(i+1, j)]
            cell.get_text().set_color(colors[i % len(colors)])
    
    plt.subplots_adjust(right=0.8)
    
    if path != False:
        fig = plt.gcf()
        fig.set_dpi(dpi)
        filename = path + str('/ROC_{}.{}'.format(tittle, format))
        plt.savefig(filename, format=format, bbox_inches='tight')
    
    if show:
        plt.show()

def pancancer_dataset_prepare(read_path, save_path, auto_sub_str=False, extension='.csv', random_state=26):
    save_path = save_path + str('/')
    file_names = [f for f in os.listdir(read_path) if f.endswith(extension)]
    if auto_sub_str != False:
        cancer_names = [f[auto_sub_str[0]:auto_sub_str[1]] for f in file_names]
    else:
        cancer_names = [f for f in file_names]
        
    file_index = -1
    for file in file_names:
        file_index += 1
        path = read_path + str('/{}'.format(file))
        dataset = DataSet(path=path, random_state=random_state)
        dataset.save(save_path, cancer_names[file_index])
    return

def auto_train(dataset_path, hyperparam_save_path, hyperparam_list, 
               models, model_names, trainer_type='random', 
               endswith='', dataset_extension='.pkl', auto_sub_str=False,
               trainer_params=''):
    log_file_path = hyperparam_save_path + str('/log_') + endswith + str('.txt')

    if os.path.exists(log_file_path):
        open(log_file_path, 'w').close()
    else:
        open(log_file_path, 'a').close()
    
    log_file = open(log_file_path, 'w')
    sys.stdout = log_file
        
    hyperparam_save_path = hyperparam_save_path + str('/')
    file_names = [f for f in os.listdir(dataset_path) if f.endswith(dataset_extension)]
    
    if auto_sub_str != False:
        cancer_names = [f[auto_sub_str[0]:auto_sub_str[1]] for f in file_names]
    else:
        cancer_names = [f for f in file_names]
    
    index = -1
    for file in file_names:
        print('****************************************')
        print('use dateset: {}'.format(file))
        print('****************************************')
        index += 1
        path = dataset_path + str('/{}'.format(file))
        with open(path, 'rb') as file:
            X_train, _, y_train, _ = pickle.load(file)
            
        for model_index in range(len(models)):
            model = models[model_index]
            
            if trainer_type == 'random':
                trainer = RandomizedSearchCV(estimator=model, 
                                             param_distributions=hyperparam_list[model_index],
                                             n_iter=500, scoring='roc_auc', 
                            n_jobs=-1, refit=True, cv=5, verbose=0, 
                            random_state=66, return_train_score=True)
            elif trainer_type == 'grid':
                trainer = GridSearchCV(estimator=model, 
                                             param_distributions=hyperparam_list[model_index],
                                             n_iter=500, scoring='roc_auc', 
                            n_jobs=-1, refit=True, cv=5, verbose=0, 
                            random_state=66, return_train_score=True)
            else:
                print('only suppert random or grid trainer')
                raise ValueError
            
            print('training {} model on {} dataset'.format(model_names[model_index], 
                                                              cancer_names[index]))
            X_train = np.array(X_train)
            y_train = np.array(y_train)
            trainer.fit(X_train, y_train)
            print('above training end')
            print('the best params is: {}'.format(trainer.best_params_))
            print('the best score is: {}'.format(trainer.best_score_))
            
            save_path = hyperparam_save_path + str('{}_{}_{}.pkl'.format(model_names[model_index],
                                                                         cancer_names[index],
                                                                         endswith))
            with open(save_path, 'wb') as file:
                pickle.dump((trainer.best_score_, trainer.best_params_), file)
                
    sys.stdout.flush()
    sys.stdout = sys.__stdout__
    log_file.close()
    return

def pancancer_ROC_curve(input_file_path, out_path, only_TCGA=False, customEndwith=False, custom_input_file_path=None):
    file_names = [f for f in os.listdir(input_file_path) if f.endswith('.pkl')]
    cancer_name = []
    for file_name in file_names:
        matched = re.search(r'data_split_(\w+).pkl', file_name)
        if matched:
            name = matched.group(1)
            cancer_name.append(name)
            
    for cancer in cancer_name:
        if customEndwith == False:
            if only_TCGA:
                path1 = '/data/jxwang_data/WMDS_lncRNA/data_output/hyperparam/' + str('RF_{}_onlyTCGA.pkl'.format(cancer))
                path2 = '/data/jxwang_data/WMDS_lncRNA/data_output/hyperparam/' + str('SVC_{}_onlyTCGA.pkl'.format(cancer))
                path3 = '/data/jxwang_data/WMDS_lncRNA/data_output/hyperparam/' + str('LR_{}_onlyTCGA.pkl'.format(cancer))
                path4 = '/data/jxwang_data/WMDS_lncRNA/data_input/ml_model/' + str('data_split_{}.pkl'.format(cancer))
                with open(path1, 'rb') as file:
                    score, param_RF = pickle.load(file)
                with open(path2, 'rb') as file:
                    score, param_SVM = pickle.load(file)
                with open(path3, 'rb') as file:
                    score, param_LR = pickle.load(file)
                # X, y
                with open(path4, 'rb') as file:
                    X_train, X_test, y_train, y_test = pickle.load(file)

                X_train = np.array(X_train)
                y_train = np.array(y_train)
                
                clf_SVM = SVC(kernel='rbf', probability=True).fit(X_train, y_train)
                clf_RF = RandomForestClassifier(random_state=26, n_jobs=-1).fit(X_train, y_train)
                clf_LR = LogisticRegression(n_jobs=-1).fit(X_train, y_train)
                
                ROC_curve(hyperparams=[param_SVM, param_RF, param_LR],
                            X=X_test, y=y_test,
                            model_names=['SVM', 'RF', 'LR'],
                            full_data=False,
                            models=[clf_SVM, clf_RF, clf_LR],
                            tittle=cancer,
                            trained=True,
                            path=out_path,
                            use_defult_hyperparams=False)
            else:
                path1 = '/data/jxwang_data/WMDS_lncRNA/data_output/hyperparam/' + str('RF_{}_withGTEx.pkl'.format(cancer))
                path2 = '/data/jxwang_data/WMDS_lncRNA/data_output/hyperparam/' + str('SVC_{}_withGTEx.pkl'.format(cancer))
                path3 = '/data/jxwang_data/WMDS_lncRNA/data_output/hyperparam/' + str('LR_{}_withGTEx.pkl'.format(cancer))
                path4 = '/data/jxwang_data/WMDS_lncRNA/data_input/ml_model/combine_GTEx/' + str('data_split_{}.pkl'.format(cancer))
                with open(path1, 'rb') as file:
                    score, param_RF = pickle.load(file)
                with open(path2, 'rb') as file:
                    score, param_SVM = pickle.load(file)
                with open(path3, 'rb') as file:
                    score, param_LR = pickle.load(file)
                # X, y
                with open(path4, 'rb') as file:
                    X_train, X_test, y_train, y_test = pickle.load(file)
                    
                X_train = np.array(X_train)
                y_train = np.array(y_train)

                clf_SVM = SVC(kernel='rbf', probability=True).fit(X_train, y_train)
                clf_RF = RandomForestClassifier(random_state=26, n_jobs=-1).fit(X_train, y_train)
                clf_LR = LogisticRegression(n_jobs=-1).fit(X_train, y_train)
                
                ROC_curve(hyperparams=[param_SVM, param_RF, param_LR],
                            X=X_test, y=y_test,
                            model_names=['SVM', 'RF', 'LR'],
                            full_data=False,
                            models=[clf_SVM, clf_RF, clf_LR],
                            tittle=cancer,
                            trained=True,
                            path=out_path,
                            use_defult_hyperparams=False)
        else:
            if custom_input_file_path == None:
                print('this model must use custom input_file_path')
                raise ValueError
            path1 = '/data/jxwang_data/WMDS_lncRNA/data_output/hyperparam/' + str('RF_{}_{}.pkl'.format(cancer, customEndwith))
            path2 = '/data/jxwang_data/WMDS_lncRNA/data_output/hyperparam/' + str('SVC_{}_{}.pkl'.format(cancer, customEndwith))
            path3 = '/data/jxwang_data/WMDS_lncRNA/data_output/hyperparam/' + str('LR_{}_{}.pkl'.format(cancer, customEndwith))
            path4 = custom_input_file_path + str('/') + str('data_split_{}.pkl'.format(cancer))
            print(path4)
            with open(path1, 'rb') as file:
                score, param_RF = pickle.load(file)
            with open(path2, 'rb') as file:
                score, param_SVM = pickle.load(file)
            with open(path3, 'rb') as file:
                score, param_LR = pickle.load(file)
            # X, y
            with open(path4, 'rb') as file:
                X_train, X_test, y_train, y_test = pickle.load(file)
                
            X_train = np.array(X_train)
            y_train = np.array(y_train)

            clf_SVM = SVC(kernel='rbf', probability=True).fit(X_train, y_train)
            clf_RF = RandomForestClassifier(random_state=26, n_jobs=-1).fit(X_train, y_train)
            clf_LR = LogisticRegression(n_jobs=-1).fit(X_train, y_train)
            
            ROC_curve(hyperparams=[param_SVM, param_RF, param_LR],
                        X=X_test, y=y_test,
                        model_names=['SVM', 'RF', 'LR'],
                        full_data=False,
                        models=[clf_SVM, clf_RF, clf_LR],
                        tittle=cancer,
                        trained=True,
                        path=out_path,
                        use_defult_hyperparams=False)
    return