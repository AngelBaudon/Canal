# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 10:11:27 2021

@author: Angel.BAUDON
"""
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy import stats
from ToolKit.IntraGrpStat import IntraGrpStat
from matplotlib.patches import RegularPolygon
from matplotlib.path import Path
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D


def PlotMaker(data, drug, drug_time, file_name, analysis,
              analysis_cell, analysis_raw, show_fig = False):
     
    Paramz = ('Resp', 'Peak borders', 'Raw', 'RegLin',
              'dFF0', 'Indexes', 'Bursts', 'Resp border')
    
    resp, borders, raw, rl, dFF0, indexes, bursts, resp_border = (data[x] for x in Paramz)
    time = np.arange(len(rl[0]))
    
    n_cell = len(data['Cell n°'])
    r, c = (n_cell+1)//2, 2
    if r == 1: r, c = 2, 1
    
   
    #Final data plotting
    mosa = plt.figure()
    mosa.suptitle(f'{drug}, {file_name}, Final traces'), mosa.subplots_adjust(hspace=1, wspace=0.5)
    
    for cell, y in zip(data['Cell n°'], dFF0):
        index, peak_border = indexes[cell], borders[cell]
        fig_plot = plt.subplot(r,c,(cell+1))
        fig_plot.axvline(drug_time, c='gold')
        fig_plot.xaxis.set_visible(False)
        fig_plot.plot(time, y, c='b', lw=0.5, zorder=0)
        for i in index: plt.scatter(time[i], y[i], c='r', marker='x')
        for L, R in bursts[cell]: fig_plot.axvspan(L, R, alpha=0.25, color='b')
        if resp[cell]:
            fig_plot.axvspan(resp_border[cell][0], resp_border[cell][1], alpha=0.1, color='r')
            fig_plot.set_facecolor('honeydew')

        #Single figures
        plt.figure(), plt.plot(time, y, c='b', zorder=0)
        plt.title(f'{drug}, {file_name}, Cell {cell}')
        plt.axvline(drug_time, c='gold'), plt.axhline(np.nanmean(y), c='k')
        for i, ind in enumerate(index):
            plt.scatter(time[ind], y[ind], c='r', s=50, marker='x', zorder=1)
            mask1, mask2 = time>peak_border[i][0], time<peak_border[i][1]
            plt.fill_between(time, y, np.nanmean(y), where=(mask1 == mask2), zorder=0)
        for L, R in bursts[cell]: plt.axvspan(L, R, alpha=0.25, color='b')
        if resp[cell]:
            fig_plot.axvspan(resp_border[cell][0], resp_border[cell][1], alpha=0.1, color='r')
            fig_plot.set_facecolor('honeydew')
        plt.savefig(f'{analysis_cell}/Cell {cell}.pdf'), plt.close()
    mosa.savefig(f'{analysis}/Final traces.pdf')
    if not show_fig: plt.close()
    


class FinalPlots():
    def __init__(self, Output, drug, General, rec_duration, drug_time):
        self.drug_time = drug_time
        self.Output = Output
        self.drug = drug
        self.rec_duration = rec_duration
        self.analysis = '\\'.join(General.split('\\')[:-1])


    def HistoMaker(self, Response, General, show_fig = False):
        
        P = ('Baseline', self.drug, 'ratios', 'means', 'sems', 'pval')
        self.Output2 = pd.DataFrame(index = P)
        
    
        def Histo(i, name, j, data, show_fig):
            bl, dr = [[np.nanmean(x[y]) for x in data[name]] for y in (0,1)]
            
            ratio = [d/(b if b!=0 else 1) for b, d in zip(bl, dr)]
            
            meanz = (np.nanmean(bl), np.nanmean(dr))
            semz = (stats.sem([x for x in bl if str(x) != 'nan']),
                    stats.sem([x for x in dr if str(x) != 'nan']))
            if j: name += ' resp'
        
            plt.figure(), plt.ylabel(name), plt.xticks((0,1), ['Baseline', self.drug])
            plt.bar([0,1], meanz, yerr=semz, color='limegreen',
                    align='center', capsize=10, zorder=0)
        
            for b, d in zip(bl, dr):
                plt.scatter(0, b, s=20, color='w', marker='o', edgecolor='k', zorder=2)
                plt.scatter(1, d, s=20, color='w', marker='o', edgecolor='k', zorder=2)
                plt.plot([0, 1], [b, d], c='k', lw=0.5, zorder=1)
            bot, top = plt.ylim()
            try: stat, pval, test = IntraGrpStat([bl, dr], Paired=True)
            except ValueError: stat, pval, test = np.nan, np.nan, 'Error'
            if pval<0.001: plt.text(0.94, -(top/8), '***', size=20, weight='bold')
            elif 0.001<pval<0.01: plt.text(0.94, -(top/7), '**', size=20, weight='bold')
            elif 0.01<pval<0.05: plt.text(0.94, -(top/7), '*', size=20, weight='bold')
            else: plt.text(0.94, -(top/7), 'ns', size=20, weight='bold')
            
            plt.text(-0.4, top-(top/6),
                     f"Test: {test} \n Stat: {stat} \n p-val: {pval} \n n: {self.C} cells over {self.S} slices")
        
            if j: plt.savefig(f'{Response}\{name}.pdf')
            else: plt.savefig(f'{General}\{name}.pdf')
            if not show_fig: plt.close()

            self.Output2[name] = [bl, dr, ratio, meanz, semz, pval]
    
    
        P = ('AUC', 'mHz', 'Height', 'HMFW', 'Rise', 'Decay', 'Area')
        self.S, self.C = len(set(self.Output['Slice'])), len(self.Output.index)
        
        
        Out_resp = self.Output.loc[self.Output["Resp"] == True]
        for j, data in enumerate((self.Output, Out_resp)):
            for n, name in enumerate(P): Histo(n, name, j, data, show_fig)

        return self.Output2

        
        
    def MosaMaker(self, General, show_fig = False):
            
        resp = []
        for s in list(set(self.Output['Slice'])):
            sl = self.Output[self.Output['Slice'] == s]
            resp.append(sum(sl['Resp'])/len(sl))
        self.Output2['Resp'] = (*(np.nan,)*2, resp, np.nanmean(resp), stats.sem(resp), np.nan)
        
        
        
        for r in ('Resp duration', 'Resp delay'):
            R = [x for x in list(self.Output[r]) if str(x) != 'nan']
            self.Output2[r] = (*(np.nan,)*2, R, np.nanmean(R), stats.sem(R), np.nan)
            
        mosa = plt.figure(figsize=(11, 7))
        mosa.suptitle(f"{self.drug}, n: {self.C} cells over {self.S} slices")
        mosa.subplots_adjust(hspace=0.5, wspace=1)
    
        for i, name in enumerate(self.Output2):
            row, col = 3, (len(self.Output2.columns)+2)//3
            subplot = plt.subplot(row, col, i+1)
            
            if name in ('Resp', 'Resp duration', 'Resp delay'):
                subplot.set_ylabel(name)
                *_, data, means, sems, _ = self.Output2[name]
                subplot.bar(0, means, yerr=sems, color='limegreen', align='center', capsize=10)
                for x in data: subplot.scatter(0, x, s=20, c='w', marker='o', edgecolor='k', zorder=2)
                subplot.xaxis.set_visible(False)
            else:
                
                bl, dr, _, meanz, semz, pval = self.Output2[name]
                subplot.set_ylabel(name), subplot.xaxis.set_visible(False)
                subplot.bar((0,1), meanz, yerr=semz, color='limegreen',
                            align='center', capsize=5, zorder=0)
                for b, d in zip(bl, dr):
                    subplot.scatter(0, b, s=10, color='w', marker='o', edgecolor='k', zorder=2)
                    subplot.scatter(1, d, s=10, color='w', marker='o', edgecolor='k', zorder=2)
                    subplot.plot([0, 1], [b, d], c='k', lw=0.5, zorder=1)
        
            if pval<=0.001: subplot.set_title('***')
            elif 0.001<pval<=0.01: subplot.set_title('**')
            elif 0.01<pval<=0.05: subplot.set_title('*')
            else: subplot.set_title('ns')    
        
        mosa.savefig(f"{self.analysis}/Mosaic sum up.pdf")
        if not show_fig: plt.close()
        return self.Output2
        
    
    def TimeCoursesMaker(self, exp_duration, bin_len = 10, show_fig = False):
        sampling_Hz = self.rec_duration//exp_duration
        for P, label in zip(('Binary', 'AUC time course'), ('n Events', 'AUC')):
            plt.figure(), plt.ylabel(label), plt.xlabel(f'Time (bin = {bin_len}s)')
            data = [[np.mean(y) for y in np.split(x, exp_duration//bin_len)]
                                for x in self.Output[P]]
            norm = [[(i/max(x) if max(x) else 0) for i in x] for x in data]
            mean, sem = np.mean(np.asarray(norm), axis=0), stats.sem(norm)
            plt.fill_between(np.arange(len(sem)), mean+sem, mean-sem,
                          facecolor = 'violet', zorder=1)
            plt.plot(mean, c='purple', zorder=1)
            [plt.plot(p, lw=.5, alpha=.5, zorder=0) for p in norm]
            plt.axvline(self.drug_time/(sampling_Hz*bin_len), c='gold', lw=2)
            plt.savefig(f'{self.analysis}\{P}.pdf')
            if not show_fig: plt.close()
            
        
    def RasterMaker(self, show_fig = False):
        plt.figure(), plt.title('sorted raster plot {}'.format(self.drug))
        plt.ylabel('Cells'), plt.xlabel('Time')
        plt.xlim(0, self.rec_duration)
        for y, cb in enumerate(self.Output['Binary']):
            [plt.scatter(x, y, s=10, c='g', marker='s') for x in np.where(cb>=1)[0]]
        plt.axvline(self.drug_time, color='gold', lw = 2)
        plt.savefig(f'{self.analysis}\Raster plot.pdf')
        if not show_fig: plt.close()

  

    def RadarMaker(self, normalizer, show_fig = False):
        names = self.Output2.columns
        meanz = self.Output2.loc['means']
        ratios = [x[0]/x[1] for x in meanz[:-3]]
        [ratios.append(x) for x in list(meanz[-3:])]
        print(ratios)
        
        def Radar(number_of_variable, frame):
            theta = np.linspace(0, 2*np.pi, number_of_variable, endpoint=False)
        
            class RadarAxes(PolarAxes):
                name, RESOLUTION = 'radar', 1
    
                def __init__(self, *args, **kwargs):
                    super().__init__(*args, **kwargs)
                    self.set_theta_zero_location('N')
        
                def fill(self, *args, closed=True, **kwargs):
                    return super().fill(closed=closed, *args, **kwargs)
        
                def plot(self, *args, **kwargs):
                    lines = super().plot(*args, **kwargs)
                    for line in lines: self._close_line(line)
        
                def _close_line(self, line):
                    x, y = line.get_data()
                    if x[0] != x[-1]: line.set_data(np.concatenate((x, [x[0]])),
                                                    np.concatenate((y, [y[0]])))
        
                def set_varlabels(self, labels):
                    self.set_thetagrids(np.degrees(theta), labels)
        
                def _gen_axes_patch(self):
                    return RegularPolygon((0.5, 0.5), number_of_variable, radius=.5, edgecolor="k")
        
                def _gen_axes_spines(self):
                    spine = Spine(axes=self, spine_type='circle',
                                  path=Path.unit_regular_polygon(number_of_variable))
                    spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                        + self.transAxes)
                    return {'polar': spine}
        
            register_projection(RadarAxes)
            return theta

        
        theta = Radar(len(ratios), 'polygon')
        fig, axes = plt.subplots(figsize=(9, 9), subplot_kw=dict(projection='radar'))
        plt.title(self.drug), axes.set_rgrids(np.linspace(0, 2, 9))
        axes.plot(theta, ratios, c='g')
        axes.fill(theta, ratios, facecolor='limegreen', alpha=0.25)
        axes.set_varlabels(names), axes.set_ylim(0,2)
        plt.savefig(f'{self.analysis}/Radar sum up.pdf')
        if not show_fig: plt.close()
        
        
            
