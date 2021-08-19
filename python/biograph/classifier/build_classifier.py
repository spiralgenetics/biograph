#!/usr/bin/python3
# pylint: disable=missing-docstring, no-self-use
"""
Builds qual classifier.
"""
import sys
import argparse

import joblib
from sklearn.ensemble import RandomForestClassifier
#import biograph.tools.log as log

def parse_args(clargs):
    """ Command line UI """
    parser = argparse.ArgumentParser(prog="build_classifier", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--train", required=True,
                        help="Table & features on which to train.")
    parser.add_argument("--out", required=True,
                        help="Output file to save the model.")
    parser.add_argument("--svonly", action="store_true",
                        help="Build a model specificially for SVs")
    args = parser.parse_args(clargs)
#    log.setup_logging()
    return args

class build_classifier:

    def __init__(self):
        pass

    def convert_object2integer(self, data):
        """
        Converts object type columns from grm.jl to numeric
        """
        dct = {'-':-1, '+':1}
        cols = ['rup_max_strand', 'rup_min_strand', 'rdn_max_strand', 'rdn_min_strand', 'aup_max_strand',
                'aup_min_strand', 'adn_max_strand', 'adn_min_strand']
        for col in cols:
            data[col] = data[col].map(dct).fillna(0)
        return(data)

    def norm_rl(self, df):
        """
        Normalizes read length dependent features
        """
        rl_feat = ["US_r", "US_a", "DS_a", "DS_r", "UXO_r", "UXO_a", "DXO_r", "DXO_a", "UMO_r", "UMO_a", "DMO_r", "DMO_a", "MO_r", "MO_a", "XO_r", "XO_a"]
        rl = df['MO_r'].max()
        df[rl_feat] = df[rl_feat] / rl
        return df

    def prepare_input(self, ret, sv_only=False):
        """
        Prepares data for training
        """
        type_d = {'DEL': -1, 'INS': 1, 'SUBSDEL': -1, 'SUBSINS': 1}
        ret['VARTYPE'] = ret['var_type'].map(type_d).fillna(0)
        ret['VAR_LEN_TYPE'] = ret['VARLEN'] * ret['VARTYPE']
        ret = ret[ret['RC'].notna()]
        if sv_only:
            ret = ret[abs(ret["VARLEN"]) >= 50]
        else:
            ret = ret[abs(ret["VARLEN"]) < 50]
        ret = self.convert_object2integer(ret)
        ret = self.norm_rl(ret)
        return(ret)

    def build_model(self, train, features, target, clf):
        """
        Return the model
        """
        clf.fit(train[features], train[target])
        return clf

    def write_model(self, model, out):
        """
        Write the model to a file
        """
        jl = {}
        jl['version'] = "biograph-7.1.0"
        jl['model'] = model
        joblib.dump(jl, out)


def main(args):
    """
    Main Runner
    """
    args = parse_args(args)
#    log.info("Loading data")
    #print("Loading data")

    train_fld = joblib.load(args.train)
    train = train_fld['data']
    features = train_fld['features']

    bc = build_classifier()

    train = bc.prepare_input(train, args.svonly)

#    log.info("Building")

    target = 'state_TP'
    clf = RandomForestClassifier(n_jobs=-1, random_state=0, n_estimators=500,
                                 bootstrap=True, max_features="sqrt", max_depth=30)
    model = bc.build_model(train, features, target, clf)

#   log.info("Writing")
    bc.write_model(model, args.out)

if __name__ == '__main__':
    main(sys.argv[1:])
