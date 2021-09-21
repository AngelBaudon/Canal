# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 09:28:01 2021

@author: Angel.BAUDON
"""

import numpy as np, scipy
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter
from scipy.interpolate import interp1d


class Analysor():
    def __init__(self, data, c, rec_duration, sampling_Hz, drug_time):
        self.c = c
        self.data = data
        self.rec_duration = rec_duration
        self.sampling_Hz = sampling_Hz
        self.drug_time = drug_time
    
    
    def envlp(self, data, chunk_range=(5,15)):
        y_new = []
        for chunk in range(*chunk_range):
            lmin = (np.diff(np.sign(np.diff(data)))>0).nonzero()[0]+1
            low = lmin[[i+np.argmin(data[lmin[i:i+chunk]])
                              for i in range(0,len(lmin),chunk)]]
            try: interp = interp1d(low, data[low], fill_value='extrapolate')
            except ValueError: pass
            y_new.append(interp(np.arange(len(data))))
        return np.nanmean(np.asarray(y_new), axis=0)


    def CellReader(self, RegLin = True, Envlp = True, baseline_rl = False):
        self.Raw = self.data[:,self.c]
        self.raw = self.Raw[:]
        self.raw = self.raw[:self.rec_duration]
        
        #Define the coeficient
        self.time = np.arange(len(self.raw))
        self.coef = self.drug_time/(len(self.raw)-self.drug_time)
        
        #Apply a linear regression
        a, b = np.polyfit(self.time, self.raw, 1)
        if RegLin: self.rl = self.raw[:] - np.asarray([a*x+b for x in self.time])
        else:
            if baseline_rl: 
                a, b = np.polyfit(self.time[:self.drug_time], self.raw[:self.drug_time], 1)
                self.rl = self.raw[:] - np.asarray([a*x+b for x in self.time])
            else: self.rl = self.raw
        
        #Apply a Polynomial local filter
        self.cel = savgol_filter(self.rl, 11, 3)
        hyperfiltr = savgol_filter(self.cel, 101, 3)
        
        #Estimate the lower envelop
        envlp_forward = self.envlp(hyperfiltr)
        envlp_reverse = np.flip(self.envlp(np.flip(hyperfiltr)))
        
        self.envlp = [x if x>y else y for x,y in zip(envlp_forward, envlp_reverse)]

        if Envlp: self.cell = self.cel - self.envlp
        else: self.cell = self.cel
        
        #Evaluate F0
        self.mode = savgol_filter(self.cell, 31, 3)
        self.F0 = np.std((self.cell-self.mode).ravel())
    
        #ΔF/F0
        self.dFF0 = (self.rl-self.F0)/self.F0

        return self.Raw, self.rl, self.dFF0, self.cell, self.coef
    
    
    def PeakDetector(self, analysis_raw, show_fig = False):
        
        #Peak detection
        self.indexes, properties = find_peaks(self.cell, height=8*self.F0,
                                              distance=10, prominence=5*self.F0)
        self.heights = list(properties.get('peak_heights'))
    
        #Plot the traces and the detected events
        plt.figure().suptitle(f'Cell n°{self.c} \n SD = {self.F0}')
        fig = plt.subplot(3,1,1)
        fig.plot(self.raw, c='b', label='Raw')
        fig.plot(self.rl, c='r', label='RegLin')
        fig.plot(self.cel, c='k', label='Filtered'), plt.legend()
        fig = plt.subplot(3,1,2)
        fig.plot(self.cel, c='b', label='Filtered')
        fig.plot(self.envlp, c='r', label='Envlp'), plt.legend()
        fig = plt.subplot(3,1,3)
        fig.plot(self.cell, c='b', label='Cell peak')
        fig.plot(self.mode, c='r', lw=.5, label='Noise estimation'), plt.legend()
        for i in self.indexes: fig.scatter(self.time[:-1][i], self.cell[i], c='r', marker='x')
        plt.savefig(f'{analysis_raw}/Cell {self.c}.pdf')
        if not show_fig: plt.close()
        
        
        return self.indexes, self.heights
    
    
    def Binary(self):
        self.cell_binary = np.zeros(self.rec_duration)
        for i in self.indexes: self.cell_binary[i]+=1
        return self.cell_binary
    
    
    def FindBurst(self, min_event = 4, max_ISI = 5*60):
        max_ISI = int(max_ISI*self.sampling_Hz)
    
        starts = [i for i, j in enumerate(self.cell_binary)
                  if j==1 and np.sum(self.cell_binary[i:(i+max_ISI)])>=min_event]

        split_indexes = [i+1 for i, s in enumerate(starts[:-1]) if (starts[i+1]-s)>max_ISI]

        starts_indexes = [j[0] for j in np.split(starts, split_indexes) if len(j)]
    

        stops_indexes = [i for i, j in enumerate(self.cell_binary)
                         if j==1 and np.sum(self.cell_binary[i:i+max_ISI])<min_event]

        # limits = [[(start, stop[min_event-2]) for stop in
        #            [x for x in stops_indexes if x>start]] for start in starts_indexes]
        limits = []
        for start_index in starts_indexes:
            possible_stop, possible_stop_index = [], (min_event-2)
            for stop_index in stops_indexes:
                if stop_index>start_index: possible_stop.append(stop_index)
            limits.append([start_index , possible_stop[possible_stop_index]])
        return limits
    
    
    def PeakAnalysor(self, show_fig=False):
        self.HM_indexes, borders, fits, rises, decays = [], [], [], [], []
        Mimir = [[[], []] for x in range(5)]
    
        local_mins_indexes = [np.where(self.cell == min(a))[0]
                                   for a in np.split(self.cell, self.indexes)]
        
        LMI = [(local_mins_indexes[0][-1],)*2]
        for lmi in local_mins_indexes[1:-1]:
            if len(lmi) == 1: LMI.append((lmi[0],)*2)
            else: LMI.append((lmi[0], lmi[-1]))
        LMI = [*LMI, (local_mins_indexes[-1][0],)*2]
        LMI = [x for y in LMI for x in y][1:-1]
        LMI = [LMI[i:i+2] for i in range(0, len(LMI), 2)]
        
        for i, index in enumerate(self.indexes):
            self.peak = self.cell[LMI[i][0]:LMI[i][1]]
            ind = index - LMI[i][0]

            Rise, Decay = self.peak[:ind], self.peak[ind:]
            HM = (max(self.peak) - min(self.peak))/2
            RMask, DMask = [x > HM for x in (Rise, Decay)]
            
            try: RM_indx = np.where(RMask == False)[0][-1]
            except IndexError: RM_indx = 0
            
            try: DM_indx = np.where(DMask == False)[0][0]
            except IndexError: DM_indx = len(DMask) - 1

            self.HM_indexes.append((ind-RM_indx, DM_indx))

        plt.figure()
        plt.plot(self.cell)
        ylim = plt.ylim()

        for HMDR, indx in zip(self.HM_indexes, self.indexes):
            HMR, HMD, StartStop, koef = indx-HMDR[0], HMDR[1]+indx, [], []
            x_ = (HMR, indx, HMD)
            y_ = (self.cell[HMR], self.cell[indx], self.cell[HMD])

            [plt.plot(x_[i], y_[i], c) for i, c in zip((0, 1, 2), ('go', 'ko', 'ro'))]
            
            for i in range(2):
                x, y, l = x_[i:i+2], y_[i:i+2], len(self.cell)
                param = np.polyfit(x, y, 1)
                fit = np.poly1d(param)
                x_ax = np.linspace(0, l, l)
                y_ax = fit(x_ax)
                plt.plot(x_ax, y_ax, lw=.8)
                
                decays.append(param[0]) if i else rises.append(param[0])
                if i:
                    try: StartStop.append((np.where(y_ax < 0)[0][0], y_ax))
                    except IndexError: StartStop.append((len(self.cell), y_ax))
                    
                else: StartStop.append((np.where(y_ax*-1 < 0)[0][0], y_ax))
            borders.append([StartStop[j][0] for j in (0, 1)])
            fits.append([StartStop[j][1] for j in (0, 1)])
        plt.ylim(ylim), plt.close()
            
        
        for i, ind in enumerate(self.indexes[:-1]):
            if borders[i][1] > borders[i+1][0]:
                Line1, Line2 = [fits[i+j][abs(j-1)] for j in (0, 1)]
                intersection = np.where(Line1 > Line2)[0][-1]
                borders[i][1], borders[i+1][0] = (intersection,)*2
                
        for i, (ind, border) in enumerate(zip(self.indexes, borders)):
            mini_peak, mini_ind = self.cell[border[0]: border[1]], ind-border[0]
                        
            #Area
            area = np.sum(np.trapz(mini_peak))
            
            #Half Max Full Width
            HMFW = len([x for x in mini_peak if x>(max(mini_peak)/2)])/self.sampling_Hz

            bldr = 0 if ind < self.drug_time else 1
            [Mimir[a][bldr].append(b) for a, b in enumerate([self.heights[i], HMFW, rises[i], decays[i], area])]


        #Separate baseline and drug peaks
        bl, dr = [int(sum(x)) for x in (self.cell_binary[:self.drug_time],
                                        self.cell_binary[self.drug_time:])]

        self.auc_cell = self.cell.clip(0)
        self.auc = (np.trapz(self.auc_cell[:self.drug_time]),
                    np.trapz(self.auc_cell[self.drug_time:])*self.coef)
        
        self.Hz = (bl/self.drug_time/self.sampling_Hz*1000,
                   dr/(self.time[-1]-self.drug_time)/self.sampling_Hz*1000)
        
        return Mimir, self.Hz, self.auc, self.auc_cell, borders
        
    
    
    def Response_finder(self, kernel = 201, show_fig = True):
        
        def Response_seeker(data, dtype):
            roll = np.convolve(data, np.ones(kernel)/kernel, mode='full')
            roll = roll[int(kernel/2):]
            
            mx = max(roll[:int(self.drug_time-kernel/2)])
            if dtype == 'hz': start = np.where(roll > mx)[0][0]
            else:
                try: start = np.where(roll > 1.2 * mx)[0][0]
                except IndexError: start = np.nan
            
            stops = np.where(roll < (mx if mx else .001))
 
            try: stop = [x for x in stops[0] if x > start + kernel][0]
            except IndexError: stop = len(roll)
            
            start, stop = start + kernel/2, stop - kernel/2

            return True, (start, stop)

        
        if (self.Hz[1]>=2*self.Hz[0] and
            self.auc[1]>=2*self.auc[0]) and np.sum(self.cell_binary)>1:
            resp, border1 = Response_seeker(self.auc_cell, 'auc')
            resp, border2 = Response_seeker(self.cell_binary, 'hz')

            border = (min(border1[0], border2[0]), max(border1[1], border2[1]))
            
        elif self.auc[1]>=2*self.auc[0] and np.sum(self.cell_binary)>1:
            resp, border = Response_seeker(self.auc_cell, 'auc')
            
        elif self.Hz[1]>=2*self.Hz[0] and np.sum(self.cell_binary)>1:
            resp, border = Response_seeker(self.cell_binary, 'hz')

        else: resp, border = False, (np.nan,)*2
        
        return resp, (border, (border[1] - border[0])/self.sampling_Hz,
                      (border[0] - self.drug_time)/self.sampling_Hz)