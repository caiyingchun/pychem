import os
import time
import joblib
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV, cross_val_predict
from sklearn.feature_selection import RFECV, VarianceThreshold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVR
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

def data_preprocessor(filename, smi_col_name, reg_col_name, test_size=0.2, set_col_name=False, Xscaler=False, Yscaler=False, descriptor='MolDesc'):
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
    data_with_des = pd.DataFrame({'SMILES': smi, reg_col_name: data[reg_col_name]})
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
    trainData = trainData[['SMILES', reg_col_name] + feat_for_model]
    testData  = testData[['SMILES', reg_col_name] + feat_for_model]
    trainData.to_csv('Train.csv', index=False)
    testData.to_csv('Test.csv', index=False)
    pd.DataFrame({'feat_for_model': feat_for_model}).to_csv('feat_for_model.csv', index=False)

    X_train = trainData.loc[:, feat_for_model]
    y_train = trainData.loc[:, reg_col_name]
    X_test = testData.loc[:, feat_for_model]
    y_test = testData.loc[:, reg_col_name]
    if Xscaler:
        X_scaler = StandardScaler().fit(X_train)
        X_train = X_scaler.transform(X_train)
        X_test = X_scaler.transform(X_test)
        joblib.dump(X_scaler, 'X_scaler.joblib')
    if Yscaler:
        with open('Yscaler.txt', 'w') as f:
            f.write('min: %s, max: %s' %(str(y_train.min()), str(y_train.max())))
        y_test = (y_test - y_train.min()) / (y_train.max() - y_train.min())
        y_train = (y_train - y_train.min()) / (y_train.max() - y_train.min())

    return X_train, y_train, X_test, y_test, feat_for_model

def modeling(filename, smi_col_name, reg_col_name, random_state=123, set_col_name=None, Xscaler=False, Yscaler=False, FP_list = ['MolDesc'], \
    k_fold=10, n_jobs=1, reg_name=['lin_reg', 'rnd_reg', 'knn_reg', 'svm_reg', 'gdb_reg']):

    lin_reg = Ridge()
    rnd_reg = RandomForestRegressor(random_state=random_state, max_features='auto')
    knn_reg = KNeighborsRegressor(weights='distance')
    svm_reg = SVR()
    gdb_reg = GradientBoostingRegressor(random_state=random_state, max_features='auto')

    resultDict = {}
    for FP in FP_list:
        resultDict[FP] = {}
        os.mkdir(FP)
        os.chdir(FP)
        X_train, y_train, X_test, y_test, feat_for_model = data_preprocessor(f'../../../{filename}', smi_col_name=smi_col_name, \
            reg_col_name=reg_col_name, set_col_name=set_col_name, Xscaler=Xscaler, Yscaler=Yscaler, descriptor=FP)
        with open('result.csv', 'w') as f:
            f.write('Model,cv_r2,cv_rmse,train_r2,train_rmse,test_r2,test_rmse\n')
        trainPredData = pd.DataFrame()
        testPredData  = pd.DataFrame()
        for regName in reg_name:
            try:
                resultDict[FP][regName] = {}
                print('-----Training_%s_%s-----\n' % (regName, FP))
                if regName == 'lin_reg':
                    reg = RFECV(estimator=eval(regName), step=int(0.1*len(feat_for_model)), cv=k_fold)
                    reg.fit(X_train, y_train)
                    with open('lin_feat_select.txt', 'w') as f:
                        i = 0
                        for each in reg.support_:
                            if each:
                                f.write(feat_for_model[i] + '\n')
                            i += 1
                
                elif regName == 'rnd_reg':
                    n = len(y_train)
                    param_rnd = {'n_estimators': [300]}
                    reg = GridSearchCV(eval(regName), param_rnd, cv=k_fold)
                    reg.fit(X_train, y_train)
                    reg = reg.best_estimator_
                    reg.fit(X_train, y_train)
                    feature_importances = reg.feature_importances_
                    feat_importance_list = list(sorted(zip(feature_importances, feat_for_model), reverse=True))
                    with open('rnd_feat_importances.txt', 'w') as f:
                        for feat_name, feat_score in feat_importance_list[:20]:
                            f.write(str(feat_name) + '\t' + str(feat_score) + '\n')

                elif regName == 'knn_reg':
                    params_knn = {'n_neighbors': [3, 5, 7]}
                    reg = GridSearchCV(eval(regName), params_knn, cv=k_fold)
                    reg.fit(X_train, y_train)
                    reg = reg.best_estimator_
                    reg.fit(X_train, y_train)

                elif regName == 'svm_reg':
                    params_svm = {'C': [1, 10, 100, 1000]}
                    reg = GridSearchCV(eval(regName), params_svm, cv=k_fold)
                    reg.fit(X_train, y_train)
                    reg = reg.best_estimator_
                    reg.fit(X_train, y_train)
                    # reg = eval(regName)
                    # reg.fit(X_train, y_train)

                else:
                    params_gdb = {'n_estimators': [50, 100], 'min_samples_leaf': [0.008, 0.005, 0.003]}
                    reg = GridSearchCV(eval(regName), params_gdb, cv=k_fold)
                    reg.fit(X_train, y_train)
                    reg = reg.best_estimator_
                    reg.fit(X_train, y_train)
                    feature_importances = reg.feature_importances_
                    feat_importance_list = list(sorted(zip(feature_importances, feat_for_model), reverse=True))
                    with open('gdb_feat_importances.txt', 'w') as f:
                        for feat_name, feat_score in feat_importance_list[:20]:
                            f.write(str(feat_name) + '\t' + str(feat_score) + '\n')

                cv_pred = cross_val_predict(reg, X_train, y_train, cv=k_fold)
                cv_r2 = r2_score(y_train, cv_pred)
                resultDict[FP][regName]['cv_r2'] = cv_r2
                cv_rmse = np.sqrt(mean_squared_error(y_train, cv_pred))
                resultDict[FP][regName]['cv_rmse'] = cv_rmse
                train_pred = reg.predict(X_train)
                train_r2 = r2_score(y_train, train_pred)
                resultDict[FP][regName]['train_r2'] = train_r2
                train_rmse = np.sqrt(mean_squared_error(y_train, train_pred))
                resultDict[FP][regName]['train_rmse'] = train_rmse
                test_pred = reg.predict(X_test)
                test_r2 = r2_score(y_test, test_pred)
                resultDict[FP][regName]['test_r2'] = test_r2
                test_rmse = np.sqrt(mean_squared_error(y_test, test_pred))
                resultDict[FP][regName]['test_rmse'] = test_rmse
                joblib.dump(reg, '%s.joblib' % regName)
                trainPredData[f'{regName}_cv_pred'] = cv_pred
                trainPredData[f'{regName}_train_pred'] = train_pred
                testPredData[f'{regName}_test_pred'] = test_pred

                with open('result.csv', 'a') as f:
                    cv_r2 = resultDict[FP][regName]['cv_r2']
                    cv_rmse = resultDict[FP][regName]['cv_rmse']
                    train_r2 = resultDict[FP][regName]['train_r2']
                    train_rmse = resultDict[FP][regName]['train_rmse']
                    test_r2 = resultDict[FP][regName]['test_r2']
                    test_rmse = resultDict[FP][regName]['test_rmse']
                    f.write('%s,%.3f,%.4f,%.3f,%.4f,%.3f,%.4f\n' % (FP+'_'+regName, cv_r2, cv_rmse, train_r2, train_rmse, test_r2, test_rmse))
            except Exception as e:
                print('*******Training_%s_%s********Error!\n' % (regName, FP))
                print(e)

        trainPredData.to_csv('trainPredData.csv', index=False)
        testPredData.to_csv('testPredData.csv', index=False)
        os.chdir('..')

    with open('result_all.csv', 'w') as f:
        f.write('Model,cv_r2,cv_rmse,train_r2,train_rmse,test_r2,test_rmse\n')
        for FP in list(resultDict.keys()):
            regNameList = list(resultDict[FP].keys())
            for model in regNameList:
                try:
                    cv_r2 = resultDict[FP][model]['cv_r2']
                except:
                    cv_r2 = 0
                try:
                    cv_rmse = resultDict[FP][model]['cv_rmse']
                except:
                    cv_rmse = 0
                try:
                    train_r2 = resultDict[FP][model]['train_r2']
                except:
                    train_r2 = 0
                try:
                    train_rmse = resultDict[FP][model]['train_rmse']
                except:
                    train_rmse = 0
                try:
                    test_r2 = resultDict[FP][model]['test_r2']
                except:
                    test_r2 = 0
                try:
                    test_rmse = resultDict[FP][model]['test_rmse']
                except:
                    test_rmse = 0
                f.write('%s,%.3f,%.4f,%.3f,%.4f,%.3f,%.4f\n' % (FP+'_'+model, cv_r2, cv_rmse, train_r2, train_rmse, test_r2, test_rmse))

def main():
    time_start = time.time()
    filename = 'delaney-processed.csv'
    smi_col_name = 'smiles'
    reg_col_name = 'measured log solubility in mols per litre'
    set_col_name = None
    Xscaler = 1
    Yscaler = 1
    k_fold = 10
    FP_list = ['AtomPairFP', 'MorganFP', 'TTFP', 'MolDesc', 'MACCSFP', 'RDKitFP']
    n_jobs = 1
    reg_name = ['rnd_reg', 'svm_reg', 'gdb_reg']  #['rnd_reg', 'knn_reg', 'svm_reg', 'lin_reg', 'gdb_reg']

    os.mkdir('regression')
    os.chdir('regression')
    for i in range(1, 4):
        random_state = np.random.randint(0, 1000)
        os.mkdir(f'{str(i)}_random_state_{random_state}')
        os.chdir(f'{str(i)}_random_state_{random_state}')
        modeling(filename=filename, smi_col_name=smi_col_name, reg_col_name=reg_col_name, set_col_name=set_col_name, Xscaler=Xscaler, Yscaler=Yscaler, k_fold=k_fold, \
            random_state=random_state, FP_list = FP_list, n_jobs=n_jobs, reg_name=reg_name)
        os.chdir('..')

    time_end = time.time()
    print('--------Totol running time is: %s s.---------' % str(int(time_end - time_start)))

if __name__ == '__main__':
    main()
