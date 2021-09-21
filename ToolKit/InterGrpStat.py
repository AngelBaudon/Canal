# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 13:04:16 2021

@author: Angel.BAUDON
"""


# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 10:38:52 2021

@author: Angel.BAUDON


To do on this code :
    - Proportions
    - LME
    - Post-hocs
    - QQ plots (& unbiased methods to estimate the fit ?)
    - Modify normality testing
    - Add a random effect / slope ?
    - Export results in Pandas
    - Check + or * for LME
    - 2 comparisons for post-hocs ?
    - Add non parametric tests

"""


# def InterGrpStat(groups, Paired=False, Proportion=False):

import scipy
import pingouin as pg
from statsmodels.stats.multicomp import pairwise_tukeyhsd, MultiComparison
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
from scikit_posthocs import posthoc_dunn
import numpy as np
import pandas as pd
from scipy.stats import sem
import matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.contingency_tables import StratifiedTable
import scipy, numpy as np, pandas as pd
from statsmodels.stats.anova import AnovaRM
from scikit_posthocs import posthoc_dunn, posthoc_conover, posthoc_tukey
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import statsmodels.api as sm
from statsmodels.formula.api import ols




groups = [np.random.normal(3, 2.5, size=(2, 20)),
          np.random.normal(3, 2.5, size=(2, 20))]
Paired = True
Proportion = False

n_cell = [len(groups[0][0]), len(groups[1][0])]
trt_name, norm  = ['TrtA', 'TrtB'], []
grp_name = [f'Time {x}' for x,y in enumerate(groups[0])]


'''               Test the criteria for parametric tests         '''
if Proportion == True:
    rest1 = [[x-np.nanmean(group) for x in group] for group in groups]
    rest2 = [x for y in rest1 for x in y]
    
else:
    rest = [[[x-np.nanmean(time) for x in time] for time in trt] for trt in groups]
    rest1 = [[x for y in trt for x in y] for trt in rest]
    rest2 = [x for y in rest1 for x in y]

#Normality
S, p_val_s = scipy.stats.shapiro(rest2)
norm.append(False) if p_val_s < 0.05 else norm.append(True)

#Equality of variances
L, p_val_l = scipy.stats.levene(*rest1)
norm.append(False) if p_val_l < 0.05 else norm.append(True)

norm = [False]


'''                             Decision tree                              '''


'''                             Prportion test                            '''
if Proportion == True:    
    '''                  With a Khi2 test                    '''
    # p1, p2 = np.array([x*100 for x in groups[0]]), np.array([x*100 for x in groups[1]])
    # dif = len(p1)-len(p2)
    # for i in range(dif): p2 = np.append(p2, np.nan)
    # # caca, pipi = np.n((1,len(p1))), np.zeros((1,len(p2)))
    # obs = np.vstack((p1, p2))
    # # obs = np.array([[14452, 4073, 4287], [30864, 11439, 9887]])                
    # chi2, pval, dof, expected = scipy.stats.chi2_contingency(obs)
    # test = "Khi2 test"
    
    '''                  With a Fiscer's exact test                '''
    # prop1, prop2 = np.nanmean(groups[0])*100, np.nanmean(groups[1])*100
    # k, pval = scipy.stats.fisher_exact([[prop1, 100-prop1], [prop2, 100-prop2]])
    # test = "Fisher's exact test"
    
    '''                       With a Z-test                          '''
    # p1, p2 = [], 0.6
    # count, nobs = np.array([1, 2]), np.array([6, 12])
    # stat, pval = proportions_ztest(6, 10, value=p2)
    
    '''           With the Cochran-Mantel-Haenszel method           '''
    # x = StratifiedTable()
    
    '''           With a linear regression           '''
    # from sklearn.datasets import load_iris
    # from sklearn.linear_model import LogisticRegression
    # X, y = load_iris(return_X_y=True)
    # clf = LogisticRegression(random_state=0).fit(X, y)
    # clf.predict(X[:2, :])



else:
    # https://www.pybloggers.com/2016/03/three-ways-to-do-a-two-way-anova-with-python/
    # https://www.researchgate.net/post/Is_there_a_non-parametric_equivalent_of_a_2-way_ANOVA
    
    #Reorganize the data
    V = [[x for y in z for x in y] for z in groups]
    C = [list(np.arange(C))*len(grp_name) for C in n_cell]    
    Ti = [np.repeat(grp_name, N) for N in n_cell]
    Tr = [np.repeat(trt_name[i], N*len(grp_name)) for i, N in enumerate(n_cell)]

    
    df = pd.DataFrame({'Cell': [x for y in C for x in y],
                       'Values': [x for y in V for x in y],
                       'Time': [x for y in Ti for x in y],
                       'Trt': [x for y in Tr for x in y]})
    
    if Paired == False:
        model = ols('Values ~ Time + Trt + Time:Trt', data=df).fit()
        aov_table = sm.stats.anova_lm(model, typ=2)
        print(aov_table)

                
    elif Paired == True:

        lme = smf.mixedlm('Values ~ Time + Trt + Time:Trt', df, groups = df['Cell'])
        '''' Attention ! je sais pas si il faut mettre un + ou un * !!! '''
        lme_out = lme.fit().summary().tables[1]
        '''' Qu'est-ce que signifie "converged : no" ??? '''
        #Why is their 'Time[T.Time1] ??? Idem for Trt
        
        
        caca = ['AOV Time', 'AOV Trt', 'AOV Interaction']
        stat = list(lme_out['z']) #is that correct that the z value is the stat ?
        pval = list(lme_out['P>|z|'])
        
        aov_out = {'Test': 'LME & Pairwise Tukey HSD'}
        for x, y, z in zip(caca[:4], stat[:4], pval[:4]):
            aov_out.update({x: [float(y), float(z)]})

        #Post-Hoc
        tukey = pairwise_tukeyhsd(endog=df['Values'], groups=df['Time'], alpha=0.05)
        tukey = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
        
        tukey2 = pairwise_tukeyhsd(endog=df['Values'], groups=df['Trt'], alpha=0.05)        
        tukey2 = pd.DataFrame(data=tukey2._results_table.data[1:], columns=tukey2._results_table.data[0])
        
        aov_out.update({'PH Time': float(tukey['p-adj']), 'PH Trt': float(tukey2['p-adj'])})
        
        
        '''Do multple tests for post-hoc '''
        
        


# print(stat, '\n', pval, '\n', ph_out)
# # try: return(stat, pval, ph_out)


