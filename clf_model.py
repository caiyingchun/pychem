import time
import os
import pickle
import joblib
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import SGDClassifier, LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.neighbors import KNeighborsClassifier
# from sklearn.naive_bayes import GaussianNB
# from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import GridSearchCV, cross_val_predict
from sklearn.metrics import confusion_matrix, accuracy_score
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import warnings
warnings.filterwarnings('ignore')

def MACCSFP(smi):
    mol = Chem.MolFromSmiles(smi)
    fp = AllChem.GetMACCSKeysFingerprint(mol)
    return [fp[i+1] for i in range(len(fp)-1)]  # 第一位为占位符

def RDKitFP(smi, fpSize=1024):
    mol = Chem.MolFromSmiles(smi)
    fp = Chem.RDKFingerprint(mol, fpSize=fpSize)
    return [bit for bit in fp]

def MorganFP(smi, radius=3, nBits=1024):
    mol = Chem.MolFromSmiles(smi)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits)
    return [bit for bit in fp]

def AtomPairFP(smi, maxLength=30, nBits=1024):
    mol = Chem.MolFromSmiles(smi)
    fp = AllChem.GetHashedAtomPairFingerprintAsBitVect(mol, maxLength=maxLength, nBits=nBits)
    return [bit for bit in fp]

def TTFP(smi, nBits=1024):
    mol = Chem.MolFromSmiles(smi)
    fp = AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=nBits)
    return [bit for bit in fp]

def MolDesc(smi):
    mol = Chem.MolFromSmiles(smi)
    nms=[x[0] for x in Descriptors._descList]
    # nms = ['MolWt', 'NumAliphaticCarbocycles', 'NumAliphaticHeterocycles', 
           # 'NumAromaticHeterocycles', 'NumAromaticRings', 'NumHAcceptors', 
           # 'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 
          # 'RingCount', 'MolLogP', 'TPSA']
    # nms = ['EState_VSA1', 'EState_VSA10', 'EState_VSA11', 
            # 'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 
            # 'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 
            # 'EState_VSA8', 'EState_VSA9', 'VSA_EState1', 
            # 'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 
            # 'VSA_EState4', 'VSA_EState5', 'VSA_EState6', 
            # 'VSA_EState7', 'VSA_EState8', 'VSA_EState9']
            
    #nms = ['NumHAcceptors', 'NumHDonors']
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
    return list(calc.CalcDescriptors(mol))

def split_train_test(data, test_size):
    shuffled_indices = np.random.permutation(len(data))
    train_set_size = int((1 - test_size) * len(data))
    train_indices = shuffled_indices[:train_set_size]
    test_indices = shuffled_indices[train_set_size:]
    return data.iloc[train_indices], data.iloc[test_indices]

def feature_preprocessor(data):
    features = data.iloc[:, 2:]

    # 去除低方差的特征
    feat_list = []
    feat_length = features.shape[-1]
    for i in range(feat_length):
        feature = features.iloc[:, i].values
        if np.mean(feature) != 0:
            if np.std(feature) > np.mean(feature) * 0.1:
                feat_list.append(i)
    feat_1 = features.iloc[:, feat_list]

    # # 去除与目标值相关系数小于0.05的特征
    # y = data.iloc[:, 1]
    # feat_1['target'] = y
    # corr_matrix = feat_1.corr()
    # feat_name_list = corr_matrix['target'].index[:-1][np.abs(corr_matrix['target'].values[:-1]) < 1].tolist()
    # feat_corr_list = np.abs(corr_matrix['target'].values[:-1][np.abs(corr_matrix['target'].values[:-1]) < 1].tolist())
    # feat_dict = dict(zip(feat_name_list, feat_corr_list))
    # feat_name_sort = sorted(feat_dict.keys(), reverse=True)
    # feat_2 = feat_1.loc[:, feat_name_sort]

    # # 去除冗余特征，保留与目标值相关性高的特征
    # def remove_redundant(list_input=feat_name_sort, feat_new=[]):
    #     if len(list_input) < 2:
    #         return (feat_new)
    #     else:
    #         feat0 = list_input[0]
    #         feat_df = features.loc[:, list_input]
    #         corr_matrix = feat_df.corr()
    #         feat_name_list = corr_matrix[feat0].index[1:][np.abs(corr_matrix[feat0].values[1:]) > 0.95].tolist()
    #         list_input = [item for item in list_input if item not in feat_name_list]
    #         feat_new.append(list_input.pop(0))
    #         return remove_redundant(list_input, feat_new=feat_new)
    # feat_name_list = remove_redundant()
    feat_name_list = feat_1.columns.tolist()
    return feat_name_list

def data_preprocessor(filename, smi_col_name, class_col_name, test_size=0.2, set_col_name=False, Xscaler=False, descriptor='MolDesc'):
    try:
        global data
        type(data) == type(pd.DataFrame())
    except:
        print('Loading data................................')
        data = pd.read_csv(filename)
    smi  = data[smi_col_name].tolist()
    if descriptor == 'AtomPairFP':
        des = list(map(lambda mol: AtomPairFP(mol), smi))
        feat_name_list = [f'AtomPairFP_{str(i+1)}' for i in range(len(des[0]))]
    elif descriptor == 'MACCSFP':
        des = list(map(lambda mol: MACCSFP(mol), smi))
        feat_name_list = [f'MACCSFP_{str(i)}' for i in range(len(des[0]))]
    elif descriptor == 'RDKitFP':
        des = list(map(lambda mol: RDKitFP(mol), smi))
        feat_name_list = [f'RDKitFP_{str(i)}' for i in range(len(des[0]))]
    elif descriptor == 'MorganFP':
        des = list(map(lambda mol: MorganFP(mol), smi))
        feat_name_list = [f'MorganFP_{str(i+1)}' for i in range(len(des[0]))]
    elif descriptor == 'TTFP':
        des = list(map(lambda mol: TTFP(mol), smi))
        feat_name_list = [f'TTFP_{str(i+1)}' for i in range(len(des[0]))]
    elif descriptor == 'MolDesc':
        des = list(map(lambda mol: MolDesc(mol), smi))
        feat_name_list = [x[0] for x in Descriptors._descList]
    else:
        print(f'The descriptor you assigned ({descriptor}) is not defined!')
        return
    data_with_des = pd.DataFrame({'SMILES': smi, class_col_name: data[class_col_name]})
    for i in range(len(feat_name_list)):
        data_with_des[feat_name_list[i]] = [des[j][i] for j in range(len(smi))]
    col_null_list = data_with_des.notnull().all(axis=0).tolist()
    data_with_des = data_with_des.iloc[:, col_null_list]
    feat_name_list = data_with_des.columns.tolist()[2:]
    trainData, testData = split_train_test(data_with_des, test_size)
    if set_col_name:
        trainData = data_with_des[data[set_col_name].isin(['Train', 'Training', 'train', 'training'])]
        testData  = data_with_des[data[set_col_name].isin(['Test', 'test'])]
    feat_for_model = feature_preprocessor(trainData)
    trainData = trainData[['SMILES', class_col_name] + feat_for_model]
    testData  = testData[['SMILES', class_col_name] + feat_for_model]
    trainData.to_csv('Train.csv', index=False)
    testData.to_csv('Test.csv', index=False)
    pd.DataFrame({'feat_for_model': feat_for_model}).to_csv('feat_for_model.csv', index=False)

    X_train = trainData.loc[:, feat_for_model]
    y_train = trainData.loc[:, class_col_name]
    X_test = testData.loc[:, feat_for_model]
    y_test = testData.loc[:, class_col_name]
    if Xscaler:
        X_scaler = StandardScaler().fit(X_train)
        X_train = X_scaler.transform(X_train)
        X_test = X_scaler.transform(X_test)
        joblib.dump(X_scaler, 'X_scaler.joblib')

    return X_train, y_train, X_test, y_test, feat_for_model

def classification_modeling(filename, smi_col_name, class_col_name, set_col_name=False, Xscaler=False, k_fold=10, random_training_num=3, n_jobs=1, \
    FP_list=['AtomPairFP', 'MorganFP', 'TTFP', 'MolDesc'], clf_name=['log_clf', 'sgd_clf', 'rnd_clf', 'svm_clf', 'knn_clf']):
    
    for i in range(1, random_training_num+1):
        random_state = np.random.randint(0, 1000)
        os.mkdir(f'{str(i)}_random_state_{random_state}')
        os.chdir(f'{str(i)}_random_state_{random_state}')

        # classification_methods
        log_clf = LogisticRegression(random_state=random_state)
        sgd_clf = SGDClassifier(random_state=random_state)
        rnd_clf = RandomForestClassifier(random_state=random_state, max_features='auto')
        svm_clf = SVC(probability=True, random_state=random_state)
        knn_clf = KNeighborsClassifier()

        conf_mx_dict = {}

        for FP in FP_list:
            conf_mx_dict[FP] = {}
            os.mkdir(FP)
            os.chdir(FP)
            X_train, y_train, X_test, y_test, feat_for_model = data_preprocessor('../../../'+filename, smi_col_name=smi_col_name, \
            class_col_name=class_col_name, set_col_name=set_col_name, Xscaler=Xscaler, descriptor=FP)
            estimators = []

            with open('result.csv', 'w') as f:
                f.write('Model,TN,FP,FN,TP,SP,SE,CA\n')

            trainPredData = pd.DataFrame()
            testPredData  = pd.DataFrame()

            for clfName in clf_name:
                print('-----Training_%s_%s-----\n' % (clfName, FP))
                try:
                    if clfName == 'rnd_clf':
                        params = {'n_estimators': [300]}
                        clf = GridSearchCV(eval(clfName), params, n_jobs=n_jobs, cv=k_fold)
                        clf.fit(X_train, y_train)
                        clf = clf.best_estimator_
                    elif clfName == 'knn_clf':
                        params = {'n_neighbors': [5, 8, 10]}
                        clf = GridSearchCV(eval(clfName), params, n_jobs=n_jobs, cv=k_fold)
                        clf.fit(X_train, y_train)
                        clf = clf.best_estimator_
                    elif clfName == 'svm_clf':
                        params = {'C': [1000]}
                        clf = GridSearchCV(eval(clfName), params, n_jobs=n_jobs, cv=k_fold)
                        clf.fit(X_train, y_train)
                        clf = clf.best_estimator_
                    else:
                        clf = eval(clfName)

                    clf.fit(X_train, y_train)
                    train_pred = cross_val_predict(clf, X_train, y_train, cv=k_fold)
                    test_pred = clf.predict(X_test)
                    conf_mx_dict[FP][clfName + f'_Train({k_fold}-CV)'] = confusion_matrix(y_train, train_pred)
                    conf_mx_dict[FP][clfName + '_Test'] = confusion_matrix(y_test, test_pred)
                    joblib.dump(clf, '%s.joblib' % clfName)
                    trainPredData[f'{clfName}_cv_pred'] = train_pred
                    testPredData[f'{clfName}_test_pred'] = test_pred
                    if accuracy_score(y_train, train_pred) >= 0.5:
                        estimators.append((clfName, clf))

                    with open('result.csv', 'a') as f:
                        for key in [clfName + f'_Train({k_fold}-CV)', clfName + '_Test']:
                            cof_matrix = conf_mx_dict[FP][key]
                            tn, fp, fn, tp = cof_matrix.ravel()
                            sp = tn * 100 / (tn + fp)
                            se = tp * 100 / (tp + fn)
                            ca = (tn + tp) * 100 / cof_matrix.sum()
                            f.write('%s,%d,%d,%d,%d,%.2f%%,%.2f%%,%.2f%%\n' % (FP+'_'+key, tn, fp, fn, tp, sp, se, ca))
                except Exception as e:
                    print('*******Training_%s_%s********Error!\n' % (clfName, FP))
                    print(e)


            if len(estimators) > 1:
                print(f'-----Training_Voting_{FP}-----\n')
                voting_clf = VotingClassifier(estimators=estimators, voting='hard')
                voting_clf.fit(X_train, y_train)
                train_pred = cross_val_predict(voting_clf, X_train, y_train, cv=k_fold)
                test_pred = voting_clf.predict(X_test)
                conf_mx_dict[FP]['voting_clf' + f'_Train({k_fold}-CV)'] = confusion_matrix(y_train, train_pred)
                conf_mx_dict[FP]['voting_clf' + '_Test'] = confusion_matrix(y_test, test_pred)
                joblib.dump(voting_clf, 'voting_clf.joblib')
                trainPredData['voting_clf_cv_pred'] = train_pred
                testPredData['voting_clf_test_pred'] = test_pred

                with open('result.csv', 'a') as f:
                    for key in ['voting_clf' + f'_Train({k_fold}-CV)', 'voting_clf' + '_Test']:
                        cof_matrix = conf_mx_dict[FP][key]
                        tn, fp, fn, tp = cof_matrix.ravel()
                        sp = tn * 100 / (tn + fp)
                        se = tp * 100 / (tp + fn)
                        ca = (tn + tp) * 100 / cof_matrix.sum()
                        f.write('%s,%d,%d,%d,%d,%.2f%%,%.2f%%,%.2f%%\n' % (FP+'_'+key, tn, fp, fn, tp, sp, se, ca))

            with open('conf_mx.pkl', 'wb') as f:
                pickle.dump(conf_mx_dict, f)

            trainPredData.to_csv('trainPredData.csv', index=False)
            testPredData.to_csv('testPredData.csv', index=False)
            os.chdir('..')

        with open('result_all.csv', 'w') as f:
            f.write('Model,TN,FP,FN,TP,SP,SE,CA\n')
            for FP in conf_mx_dict.keys():
                for key in conf_mx_dict[FP].keys():
                    cof_matrix = conf_mx_dict[FP][key]
                    tn, fp, fn, tp = cof_matrix.ravel()
                    sp = tn * 100 / (tn + fp)
                    se = tp * 100 / (tp + fn)
                    ca = (tn + tp) * 100 / cof_matrix.sum()
                    f.write('%s,%d,%d,%d,%d,%.2f%%,%.2f%%,%.2f%%\n' % (FP+'_'+key, tn, fp, fn, tp, sp, se, ca))

        os.chdir('..')

def main():
    time_start = time.time()
    filename = 'BBB.csv'  # 数据文件名
    smi_col_name = 'Structure'  # SMILES 列名
    class_col_name = 'Class'  # 类别 列名
    set_col_name = None  # 列表中是否有明确的列划分训练集和测试集，如有，指定该列名
    Xscaler = 0  # 特征是否标准化
    k_fold = 10  # 交叉验证倍数
    random_training_num = 5  # 随机训练次数
    n_jobs = 1  # 训练核数
    FP_list = ['AtomPairFP', 'MorganFP', 'TTFP', 'MolDesc', 'MACCSFP', 'RDKitFP']  # 指定描述符或指纹: ['AtomPairFP', 'MorganFP', 'TTFP', 'MolDesc', 'MACCSFP', 'RDKitFP']
    clf_name = ['log_clf', 'sgd_clf', 'rnd_clf', 'svm_clf', 'knn_clf']  # 指定分类器: ['log_clf', 'sgd_clf', 'rnd_clf', 'svm_clf', 'knn_clf']
    

    os.mkdir('classification')
    os.chdir('classification')
    classification_modeling(filename=filename, smi_col_name=smi_col_name, class_col_name=class_col_name, set_col_name=set_col_name, Xscaler=Xscaler, \
        k_fold=k_fold, random_training_num=random_training_num, n_jobs=n_jobs, FP_list=FP_list, clf_name=clf_name)
    os.chdir('..')

    time_end = time.time()
    print('--------Totol running time is: %s s.---------' % str(int(time_end - time_start)))

if __name__ == '__main__':
    main()
