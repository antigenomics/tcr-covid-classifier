import warnings

import pandas as pd
from imblearn.pipeline import Pipeline
from sklearn.metrics import f1_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier

warnings.filterwarnings('ignore')


def get_model(model):
    '''
    Get the model from six state-of-the-art machine learning models.
    '''
    if model == 'svm':
        from sklearn.svm import SVC
        names = ["Linear SVM"]
        classifiers = [
            SVC()
        ]
    elif model == 'ab':
        from sklearn.ensemble import AdaBoostClassifier
        names = ["AdaBoost"]
        classifiers = [
            AdaBoostClassifier()
        ]
    elif model == 'knn':
        from sklearn.neighbors import KNeighborsClassifier
        names = ["K-Nearest Neighbors"]
        classifiers = [
            KNeighborsClassifier()
        ]
    elif model == 'rfc':
        from sklearn.ensemble import RandomForestClassifier
        names = ["Random Forest"]
        classifiers = [
            RandomForestClassifier()
        ]
    elif model == 'xgboost':
        from xgboost import XGBClassifier
        names = ["XGBoost"]
        classifiers = [
            XGBClassifier()
        ]
    elif model == 'mlpclassifier':
        from sklearn.neural_network import MLPClassifier
        names = ["MLPClassifier"]
        classifiers = [
            MLPClassifier()
        ]
    else:
        raise RuntimeError('Unknown classifier')

    return classifiers


def get_parameters():
    return {
        'svm': {'model__C': (1, 5, 10, 50, 100),
                'model__kernel': ('linear', 'poly', 'rbf', 'sigmoid'),
                'model__probability': [True],
                'model__random_state': [42]},
        'ab': {'model__n_estimators': (10, 25, 50, 100, 125, 150, 200),
               'model__base_estimator': (DecisionTreeClassifier(max_depth=1), DecisionTreeClassifier(max_depth=3)),
               'model__random_state': [42]},
        'knn': {'model__n_neighbors': (3, 5, 10, 50, 75, 100),
                'model__leaf_size': (1, 2, 3, 5, 10, 15, 20),
                'model__weights': ['uniform', 'distance']},
        'rfc': {'model__max_depth': (1, 2, 3),
                'model__n_estimators': (50, 75, 100, 125, 150, 200),
                'model__min_samples_leaf': (8, 10),
                'model__oob_score': (False, True),
                'model__random_state': [42]},
        'mlpclassifier': {'model__hidden_layer_sizes': (
            (100, 60, 30, 10),
            (150, 100, 50, 25, 10),
            (100, 50, 25, 10)),
            'model__alpha': (0.0001, 0.001, 0.00001, 0.01),
            'model__learning_rate': ['constant', 'adaptive'],
            'model__max_iter': [1000],
            'model__random_state': [42]
        },
        'xgboost': {'model__n_estimators': (10, 25, 50, 75, 100),
                    'model__subsample': (0.25, 0.5, 0.75, 1),
                    'model__random_state': [42]}
    }


def prepare_run_column(df):
    df['run'] = df['run'].apply(lambda x: x.split('.')[0])
    return df


def prepare_data(run_to_number_of_clones_path,
                 desc_path,
                 clonotype_matrix_path,
                 hla_keys_path=None,
                 make_freq=True,
                 use_hla_clones=True,
                 use_hla_bool=True,
                 make_all_features_bool=False,
                 use_standardization=True,
                 raw_target_column='COVID_status',
                 raw_target_clumn_success_label='COVID',
                 final_target_column='covid'):
    run_to_number_of_clones = prepare_run_column(pd.read_csv(run_to_number_of_clones_path))
    columns_to_save = []
    desc = prepare_run_column(pd.read_csv(desc_path))
    if 'folder' in desc.columns:
        columns_to_save = ['run', raw_target_column, 'folder']
    else:
        columns_to_save = ['run', raw_target_column]
    desc = desc[columns_to_save]
    desc = desc[desc[raw_target_column] != 'unknown']
    desc[raw_target_column] = desc[raw_target_column].apply(lambda x: 1 if raw_target_clumn_success_label == x else 0)

    df = prepare_run_column(pd.read_csv(clonotype_matrix_path).drop(columns=['Unnamed: 0']))
    if use_hla_clones:
        if hla_keys_path is not None:
            hla_keys = pd.read_csv(hla_keys_path)['0']
            for hla in hla_keys:
                hla_df = prepare_run_column(pd.read_csv(
                    f'data/hla_sign_clone_matrix/hla_covid_clonotype_matrix_500k_top_1_mismatch_hla_{hla}.csv').drop(
                    columns=['Unnamed: 0']))
                df = df.merge(hla_df)
        else:
            raise Exception('HLA keys path essential!')

    df = df.merge(run_to_number_of_clones)
    for col in df.columns:
        if col != 'run':
            if make_freq:
                df[col] = df[col] / df['number_of_clones']
            elif make_all_features_bool:
                df[col] = df[col].apply(lambda x: x > 0)
    desc = desc.merge(df.drop(columns=['number_of_clones']))

    if use_standardization:
        data = pd.DataFrame(
            data=StandardScaler().fit_transform(desc.drop(columns=columns_to_save)),
            columns=desc.columns[len(columns_to_save):])
    else:
        data = desc.drop(columns=columns_to_save)

    if use_hla_bool:
        if hla_keys_path is not None:
            hla_keys = pd.read_csv(hla_keys_path)['0']
            for hla in hla_keys:
                hla_runs = prepare_run_column(pd.read_csv(f'data/hla_desc/fmba_desc_hla_{hla}.csv'))
                data[hla] = desc.run.isin(hla_runs.run)
        else:
            raise Exception('HLA keys path essential!')

    data[final_target_column] = desc[raw_target_column]
    if 'folder' in columns_to_save:
        data['folder'] = desc['folder']
    return data


def evaluate_models(X_train, y_train, X_test, y_test, model_params, scoring_function='f1', debug=False):
    # TODO add smote support here
    skf = StratifiedKFold(n_splits=5, random_state=42, shuffle=True)
    best_clfs_encode_fmba = {}
    scores = []
    params = []
    model_names = list(model_params.keys())
    for name, param in model_params.items():
        if debug:
            print(f'Started evaluating {name}')
        model = get_model(name)[0]

        steps = [('model', model)]

        pipeline = Pipeline(steps=steps)
        clf = GridSearchCV(pipeline, param, scoring=scoring_function, cv=skf, n_jobs=-1)
        clf.fit(X_train, y_train)

        pipe = clf.best_estimator_['model']

        scores.append(f1_score(y_test, pipe.predict(X_test)))
        params.append(clf.best_params_)
        if debug:
            print(f'Best params for {name}:', clf.best_params_)
            print('Test f1-score for the best model %.2f' % f1_score(y_test, pipe.predict(X_test)))
            print()

        best_clfs_encode_fmba[name] = clf.best_estimator_
    print(pd.DataFrame({'classifier': model_names,
                        'f1-score': scores
                        }).set_index('classifier').T)
    best_model_idx = max(range(len(scores)), key=lambda i: scores[i])
    print(f'Best model is {model_names[best_model_idx]} with params: {params[best_model_idx]}')

    return best_clfs_encode_fmba, scores[best_model_idx], model_names[best_model_idx]


def split_data_by_batch(data, test_batches, y_column, batch_column='folder'):
    train_data = data[~data[batch_column].isin(test_batches)]
    test_data = data[data[batch_column].isin(test_batches)]

    y_train = train_data[y_column]
    X_train = train_data.drop(columns=[batch_column, y_column])

    y_test = test_data[y_column]
    X_test = test_data.drop(columns=[batch_column, y_column])

    return X_train, y_train, X_test, y_test


def cross_validation_between_batches(model, data, valid_batches, metrics=f1_score, metrics_name='f1_score', return_metrics_results=False,
                                     y_column='covid', batch_column='folder', debug=False):
    f1_scores = []
    for folder in valid_batches:
        X_train, y_train, X_test, y_test = split_data_by_batch(data, [folder], y_column, batch_column=batch_column)

        model.fit(X_train, y_train)
        if debug:
            print(folder)
            print(f'Valid {metrics_name} for the model %.2f' % metrics(y_test, model.predict(X_test)))
            print()
        f1_scores.append(metrics(y_test, model.predict(X_test)))
    f1_scores = pd.Series(f1_scores)
    if not return_metrics_results:
        print(f'Final validation {metrics_name} is %.2f±%.2f' % (f1_scores.mean(), f1_scores.std()))
        return f1_scores.mean(), f1_scores.std()
    else:
        return pd.DataFrame({'folder': valid_batches, metrics_name: f1_scores})


def make_hla_predictor(hla, run_to_number_of_clones_path, desc_path, clonotype_matrix_path, debug=False):
    print(f'{hla} predictor')
    data = prepare_data(run_to_number_of_clones_path=run_to_number_of_clones_path,
                        desc_path=desc_path,
                        clonotype_matrix_path=clonotype_matrix_path,
                        hla_keys_path=None,
                        make_freq=True,
                        use_hla_clones=False,
                        use_hla_bool=False,
                        make_all_features_bool=False,
                        use_standardization=True,
                        raw_target_column=hla,
                        raw_target_clumn_success_label='present',
                        final_target_column='allele')
    allele_percent = data.allele.sum() / (data.shape[0]) * 100
    features_count = data.shape[1]
    if allele_percent < 1:
        return allele_percent, features_count, None, 0, 0, 0
    X_train, y_train, X_test, y_test = split_data_by_batch(data=data,
                                                           test_batches=['2020/10_FMBA_NovaSeq6'],
                                                           y_column='allele',
                                                           batch_column='folder')

    clfs, best_score, model_name = evaluate_models(X_train, y_train, X_test, y_test, get_parameters(),
                                                   scoring_function='f1')

    cv_f1, cv_std = cross_validation_between_batches(clfs[model_name], data, data.folder.unique(), y_column='allele',
                                                     batch_column='folder')
    if debug:
        print(f'Features count is {data.shape[1]}')
        print(f'There are {data.allele.sum()} samples with {hla}')
        print(f'There are {data.shape[0] - data.allele.sum()} samples without {hla}')
        print(f'Percent of data with {hla} is {round(allele_percent, 2)}%')
    print()
    return allele_percent, features_count, clfs[model_name], best_score, cv_f1, cv_std