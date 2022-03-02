import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from sklearn.datasets import load_wine
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import FunctionTransformer
from sklearn.pipeline import Pipeline
from sklearn.pipeline import make_pipeline

from sklearn.model_selection import train_test_split
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import make_scorer
from sklearn.ensemble import GradientBoostingClassifier

import seaborn as sns
import sklearn

import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import sys


#dat_file = sys.argv[1]
#bp = sys.argv[2]
#outdir = sys.argv[3]
dat_file = '/stanley/levin_dr/kwanho/projects/Amanda/Atlas/share-seq_URD/RNA/downstream/branchpoint_DE/data_varGenesLasso+DEG_bp1.tsv'
bp = "bp1"
outdir = '/stanley/levin_dr/kwanho/projects/Amanda/Atlas/share-seq_URD/RNA/downstream/branchpoint_DE/classifier/'

print(dat_file)
print(bp)

dat = pd.read_csv(dat_file,'\t')

print(dat.head())
X=dat.iloc[:,1:-1] # remove cell name column

dat['class'] = "segment" + dat['class'].astype(str)
y=dat['class'].astype('category').cat.codes
print(dict( enumerate(dat['class'].astype('category').cat.categories ) ))
y=y.to_frame()
y.columns=['label']
print(y.head())
print(y.label.value_counts())
print(X.shape,y.shape)

pipe=make_pipeline(GradientBoostingClassifier())
print(pipe.named_steps)
print(X.values.shape)

nfolds=10
param_grid = {'gradientboostingclassifier__n_estimators': [25, 50, 75, 100], 'gradientboostingclassifier__max_depth': [3,4,5]}
search = GridSearchCV(pipe, param_grid, iid=False, cv=StratifiedKFold(nfolds), n_jobs=8, verbose=10)
search.fit(X.values,y.values.ravel())

print(search.cv_results_)
print(search.best_score_, search.best_params_)

classes=dict( enumerate(dat['class'].astype('category').cat.categories ) )
print(classes)

max_depth=search.best_params_['gradientboostingclassifier__max_depth']
n_estimators=search.best_params_['gradientboostingclassifier__n_estimators']
n_feats= X.values.shape[1]
#optpipe=make_pipeline(StandardScaler(with_std=False),GradientBoostingClassifier(max_depth=4, n_estimators=50))
optpipe=make_pipeline(GradientBoostingClassifier(max_depth=max_depth, n_estimators=n_estimators))

cvres=cross_validate(pipe,X,y.squeeze(),cv=nfolds,scoring='f1_micro', return_estimator=True, return_train_score=True)
print(cvres['test_score'])
print(cvres['train_score'])

feat_imps=np.zeros((10,len(classes),n_feats))

for fold in range(nfolds):
    for c in range(len(classes)):
        for i in range(n_estimators):
            feat_imps[fold,c]+=cvres['estimator'][fold].named_steps.gradientboostingclassifier.estimators_[i,c].tree_.compute_feature_importances(normalize=False)
        feat_imps[fold,c]/=feat_imps[fold,c].sum()


plt.figure(num=None, figsize=(24, 6), dpi=80, facecolor='w', edgecolor='k');
plt.boxplot(feat_imps[:,0,:]);
plt.xticks(range(1,n_feats+1),dat.columns[1:-1],rotation='vertical')
plt.grid();
plt.title(classes[0]);
plt.savefig(outdir + bp + '_' + classes[0] + '.png')

plt.figure(num=None, figsize=(24, 6), dpi=80, facecolor='w', edgecolor='k');
plt.boxplot(feat_imps[:,1,:]);
plt.xticks(range(1,n_feats+1),dat.columns[1:-1],rotation='vertical')
plt.grid();
plt.title(classes[1]);
plt.savefig(outdir + bp + '_' + classes[1] + '.png')

plt.figure(num=None, figsize=(24, 6), dpi=80, facecolor='w', edgecolor='k');
plt.boxplot(feat_imps[:,2,:]);
plt.xticks(range(1,n_feats+1),dat.columns[1:-1],rotation='vertical')
plt.grid();
plt.title(classes[2]);
plt.savefig(outdir + bp + '_' + classes[2] + '.png')

print(feat_imps.shape)
print(n_feats)

df = dat.columns[1:-1]

for i in range(len(classes)):
    name = classes[i]
    np.savetxt(outdir + bp + '_' + name + '.csv', feat_imps[:,i,:], delimiter=",")



