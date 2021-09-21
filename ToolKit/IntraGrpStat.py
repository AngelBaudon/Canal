# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 13:05:09 2021

@author: Angel.BAUDON


To do on this code :
    - WARNING ! Every pval are corrected with yhe bonferoni method wich is the most conservative !
    - QQ plots (& unbiased methods to estimate the fit ?)
    - Test sphericity ?
    
"""
def IntraGrpStat(groups, Paired=False):
    
    import scipy, numpy as np, pandas as pd
    from statsmodels.stats.anova import AnovaRM
    from scikit_posthocs import posthoc_dunn, posthoc_conover, posthoc_tukey
    from statsmodels.stats.multicomp import pairwise_tukeyhsd

    group_name = [f'Time {x}' for x,y in enumerate(groups)]
    n_cell, norm = len(groups[0]), []
    
    
    '''               Test the criteria for parametric tests         '''
    rest1 = [[x-np.nanmean(group) for x in group] for group in groups]
    rest2 = [x for y in rest1 for x in y]
    
    
    #Normality
    S, p_val_s = scipy.stats.shapiro(rest2)
    norm.append(False) if p_val_s < 0.05 else norm.append(True)
    
    #Equality of variances
    L, p_val_l = scipy.stats.levene(*rest1)
    norm.append(False) if p_val_l < 0.05 else norm.append(True)
    
    
    
    '''                                   Decision tree                 '''
    #T-test familly
    if len(groups) == 2:
        #Parametric test
        if not False in norm:
            if Paired == True:
                stat, pval = scipy.stats.ttest_rel(*groups)
                test = 'Paired t-Test'
            else:
                stat, pval = scipy.stats.ttest_ind(*groups)
                test = 'Unpaired t-Test'
                
                
        #Non parametric test
        if False in norm:
            if Paired == True:
                stat, pval = scipy.stats.wilcoxon(*groups)
                test = 'Wilcoxon'
            else:
                stat, pval = scipy.stats.mannwhitneyu(*groups)
                test = 'Mann-Whitney'
          
    
    #Anova familly
    elif len(groups) > 2:
        #Reorganize the data
        df = pd.DataFrame({'Cell': [x for x in range(n_cell)]*len(groups),
                           'Values': [x for y in groups for x in y],
                           'Time': np.repeat(group_name, repeats=n_cell)})
        comp1, comp2 = [*[group_name[0]]*2, group_name[1]], [group_name[1], *[group_name[2]]*2]
        
        if Paired == False:
            if not False in norm:
                stat, pval = scipy.stats.f_oneway(*groups)
                ph = posthoc_tukey(df, val_col='Values', group_col='Time')
                
                ph_out = {'Test': 'One-Way ANOVA & Tukey'}
                for x, y in zip(comp1, comp2):
                    ph_out[f'{x} vs {y}'] = float(ph.loc[x, y])
                    
            else:
                stat, pval = scipy.stats.kruskal(*groups)
                ph = posthoc_dunn(groups, p_adjust = 'bonferroni')
                
                ph_out = {'Test': "Kruskal-Wallis & Dunn's MC"}
                for x, y, c1, c2 in zip([1,1,2], [2,3,3], comp1, comp2):
                    ph_out[f'{c1} vs {c2}'] = float(ph.loc[x, y])
                        
            
            
            
        elif Paired == True:
            if not False in norm:
                aovrm = AnovaRM(df, depvar='Values', subject='Cell', within=['Time'])
                res = aovrm.fit().summary().tables[0]
                stat, pval = float(res['F Value']), float(res['Pr > F'])
    
                # ph = pairwise_ttests(data=df, dv='Values', within='Time',
                #                      subject='Cell', padjust='bonf')
    
                ph = pairwise_tukeyhsd(df['Values'], df['Time'])
                
                ph = pd.DataFrame(data=ph._results_table.data[1:],
                                  columns=ph._results_table.data[0])
                ph_out = {'Test': 'RM ANOVA & paired Tukey'}
                for x, y, z in zip(*[list(ph[x]) for x in ['group1', 'group2', 'p-adj']]):
                    ph_out[f'{x} vs {y}'] = z
    
            
            else: 
                stat, pval = scipy.stats.friedmanchisquare(*groups)
                ph = posthoc_conover(groups, p_adjust = 'bonferroni')
                
                ph_out = {'Test': 'Friedman & Wilcoxon with Bonferoni correction'}
                for x, y, c1, c2 in zip([1,1,2], [2,3,3], comp1, comp2):
                    ph_out[f'{c1} vs {c2}'] = float(ph.loc[x, y])
    
    try: return(round(stat, 4), round(pval, 4), ph_out)
    except NameError: return(round(stat, 4), round(pval, 4), test)
