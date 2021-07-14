#!/usr/bin/python3
# pylint: disable=missing-docstring, no-self-use
"""
Builds and tests a classifier
Uses only 35 features , selected using rf feature importance
"""
import sys
import argparse

import joblib
import pandas as pd
import numpy as np

from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, matthews_corrcoef
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix

import biograph.tools.log as log

def parse_args(clargs):
    """ Command line UI """
    parser = argparse.ArgumentParser(prog="build_classifier", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--train", required=True,
                        help="Table on which to train.")
    parser.add_argument("--test", required=True,
                        help="Table on which to test.")
    parser.add_argument("--out", required=True,
                        help="Output file to save the model.")
    parser.add_argument("--stats", default=None,
                        help="File to write model benchmarking stats. (stdout)")
    parser.add_argument("--svonly", action="store_true",
                        help="Build a model specificially for SVs")
    args = parser.parse_args(clargs)
    log.setup_logging()
    return args


class build_classifier:

    def __init__(self):
        pass

    def set_feat_tar_estimator(self):
        """
        Set your features and target variable
        """
        features = ['VARLEN', 'NUMASM', 'OV', 'DP', 'AD_r', 'AD_a', 'PDP', 'PAD_a', 'UC_r', 'UC_a',
                    'DC_a', 'UDC_r', 'UDC_a', 'UCC_r', 'UCC_a', 'DCC_a', 'DMO_a', 'NR_r', 'NR_a', 'MO_a', 'XC_r',
                    'XC_a', 'AC_a', 'MC_a', 'EC_a', 'PL_ref', 'PL_hom', 'SCORE', 'REFSPAN', 'ASMLEN', 'LANCH',
                    'RANCH', 'REFGC', 'ALTGC', 'POP_True']

        target = 'state_TP'

        estimator = RandomForestClassifier(n_jobs=-1, random_state=0, n_estimators=500,
                                           bootstrap=True, max_features=0.3, max_depth=35)

        return features, target, estimator

    def one_hot_encoding(self, ret, categorical_features):
        for col in categorical_features:
            for_dummy = ret.pop(col)
            ret = pd.concat([ret, pd.get_dummies(for_dummy, prefix=col)], axis=1)
        return ret

    def prepare_input(self, df, sv_only=False):

        ret = df.copy()
        ret.drop_duplicates(inplace=True)

        ## Identify categorical variables and encode them.

        state_d = {'FP':0, 'TP':1}
        POP_d = {False:0, True:1}
        ret['POP_True'] = ret['POP'].map(POP_d).fillna(ret['POP'])
        ret['state_TP'] = ret['state'].map(state_d).fillna(ret['state'])

        if sv_only:
            ret = ret[abs(ret["VARLEN"]) >= 50]
        else:
            ret = ret[abs(ret["VARLEN"]) < 50]

        return(ret)


    def build_model(self, train, features, target, clf):
        """
        Return the model
        """
        clf.fit(train[features], train[target])
        return clf

    def get_performance_metrics(self, y_test, y_predicted):
        """
        Basic performance metrics over test data
        """
        accuracy = accuracy_score(y_test, y_predicted)
        precision = precision_score(y_test, y_predicted)
        recall = recall_score(y_test, y_predicted)
        f1 = f1_score(y_test, y_predicted)
        mcc = matthews_corrcoef(y_test, y_predicted)
        tn, fp, fn, tp = confusion_matrix(y_test, y_predicted).ravel()
        fpr = fp/(fp+tn)
        return (accuracy, precision, recall, mcc, fpr, f1, tn, fp, fn, tp)

    def get_cv_metric(self, train, features, target, estimator, n_fold=5): # pylint: disable=too-many-locals
        """
        Performs n-fold cross-validation on training data
        """

        ret_data = []

        X = train[features].values
        y = train[target].values

        skf = StratifiedKFold(n_splits=n_fold, random_state=1)
        s1 = ("\nResults of {} fold cross-validation (Stratified KFold) using threshold 0.5\n".format(n_fold))
        s2 = ('{:16s}{:8s}{:8s}{:8s}{:8s}{:8s}{:8s}{:8s}'.format('Fold', 'TPR', 'PPV', 'FPR', 'F1', 'MCC', 'ACC', 'ACC_train'))
        ret_data.append(s1)
        ret_data.append(s2)

        train_sum = 0
        test_sum = 0
        MCC_sum = 0
        F1_score_sum = 0
        prec_sum = 0
        recall_sum = 0
        AUC_sum = 0
        FPR_sum = 0
        k = 0

        for k, (train_index, test_index) in enumerate(skf.split(X, y)):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            clf = estimator.fit(X_train, y_train)

            accuracy_train = clf.score(X_train, y_train)
            train_sum = train_sum + accuracy_train

            y_pred = clf.predict(X_test)
            y_prob = clf.predict_proba(X_test)

            accuracy_test, PPV, TPR, MCC, FPR, F1_score, tn, fp, fn, tp = self.get_performance_metrics(y_test, y_pred) # pylint: disable=unused-variable
            AUC = roc_auc_score(y_test, y_prob[:, 1])

            test_sum = test_sum + accuracy_test
            F1_score_sum += F1_score
            MCC_sum += MCC
            prec_sum += PPV
            recall_sum += TPR
            AUC_sum += AUC
            FPR_sum += FPR

            s3 = ('{:<16d}{:<8.3f}{:<8.3f}{:<8.3f}{:<8.3}{:<8.3f}{:<8.3}{:<8.3f}'.format(
                k+1, TPR, PPV, FPR, F1_score, MCC, accuracy_test, accuracy_train
            ))

            ret_data.append(s3)


        s4 = ('{:16s}{:<8.3f}{:<8.3f}{:<8.3f}{:<8.3}{:<8.3f}{:<8.3}{:<8.3}'.format(
            "Avg",
            recall_sum / (k + 1),
            prec_sum / (k + 1),
            FPR_sum / (k + 1),
            F1_score_sum / (k + 1),
            MCC_sum / (k + 1),
            test_sum / (k + 1),
            train_sum / (k + 1)
        ))

        ret_data.append(s4)
        return ret_data


    def run_classifier(self, classifier, data, features=None):
        """
        Run a classfier on data with appropriate features
        """
        return classifier.predict_proba(data[features])

    def benchmark_stats(self, predicted_probs, test, target, thresholds=np.linspace(0.15, 0.5, num=35), svonly=False): # pylint: disable=too-many-locals

        ret_data = []
        header = ('{:<14s}{:<8s}{:<8s}{:<8s}{:<8s}{:<8s}{:<8s}{:<8s}{:<10s}{:<10s}{:<10s}{:<10s}'.format(
            'type', 'thresh', 'TPR', 'PPV', 'FPR', 'F1', 'MCC', 'ACC', 'tp', 'fn', 'tn', 'fp'
        ))

        ret_data.append(header)

        for threshold in thresholds:

            predicted = predicted_probs[:, 1] >= threshold
            actual = test[target].array

            if svonly:
                category = 'SV'
            else:
                category = 'VARLEN < 50'

            ACC, PPV, TPR, MCC, FPR, F1, tn, fp, fn, tp = self.get_performance_metrics(actual, predicted)
            stats = ('{:<14s}{:<8.3f}{:<8.3f}{:<8.3f}{:<8.3f}{:<8.3}{:<8.3f}{:<8.3f}{:<10d}{:<10d}{:<10d}{:<10d}'.format
                     (category, threshold, TPR, PPV, FPR, F1, MCC, ACC, tp, fn, tn, fp))
            ret_data.append(stats)

        return(ret_data)

    def show_feature_importance(self, clf, train, features):
        """
        print("# View a list of the features and their importance scores")
        """
        ret = ["Feature Importance\nFeature\tRank\n"]
        data = list(zip(train[features], clf.feature_importances_))
        data.sort(key=lambda x: x[1], reverse=True)
        for feat, score in data:
            ret.append("%s\t%.4f\n" % (feat, score))
        return "".join(ret) + '\n'


    def benchmark(self, probs, model, train, test, features, target, svonly=False, output=None):
        """
        Given a model, its training data, its test data, its features, and its target,
        output benchmark-information
        """
        with open(output, 'w') as fh:
            fh.write('\n{:s}\n'.format("#"*100))
            fh.writelines("%s\n" % place for place in self.get_cv_metric(train, features, target, model, n_fold=3))
            fh.write('\n{:s}\n'.format("#"*100))
            fh.write('\nPerformance on the test set\n\n')
            fh.writelines("%s\n" % place for place in self.benchmark_stats(probs, test, target, svonly=svonly))
            fh.write('\n{:s}\n'.format("#"*100))
            fh.write(self.show_feature_importance(model, train, features))

    def write_model(self, model, out):
        """
        Write the model to a file
        """
        jl = {}
        jl['version'] = "biograph_sdk-5.0.0"
        jl['model'] = model
        joblib.dump(jl, out)


def main(args):
    """
    Main Runner
    """
    args = parse_args(args)
    log.info("Loading data")
    #print("Loading data")

    train = joblib.load(args.train)
    train = train['data']

    test = joblib.load(args.test)
    test = test['data']

    bc = build_classifier()

    train = bc.prepare_input(train, args.svonly)
    test = bc.prepare_input(test, args.svonly)

    #joblib.dump(train, "train.ml")
    #joblib.dump(test, "test.ml")

    log.info("Building")
    features, target, clf = bc.set_feat_tar_estimator()
    model = bc.build_model(train, features, target, clf)

    log.info("Testing")
    probs = bc.run_classifier(model, test, features)

    log.info("Benchmarking")
    bc.benchmark(probs, model, train, test, features, target, args.svonly, args.stats)

    log.info("Writing")
    bc.write_model(model, args.out)


if __name__ == '__main__':
    main(sys.argv[1:])
