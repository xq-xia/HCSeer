#!/usr/bin/python3

### Made by xiaxingquan
### March 2024



# coding = utf-8
import random

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import LeaveOneOut,GridSearchCV
from scipy.stats import gaussian_kde  # 用于计算概率密度

########################################################################
######     The purpose of this class is to calculate best fit       ####
######     bandwidth and corresponding kernel density function      ####
########################################################################

class KDE:
    def __init__(self, data):
        self.data = data
        self.n = len(data)
        self.defactor = self.density_factor()
        self.hopt = self.h_opt()

    ## get average of all data
    def get_average(self):
        x_average = 0
        for i in self.data:
            x_average += i
        # print(x_average / self.n)
        return x_average/self.n
            
    ## kernel funtion
    def gaussian_kernel(self,x):
        """
        gaussian funtion
        """
        return (1 / np.sqrt(2 * np.pi)) * np.exp(-0.5 * x**2)
    
    ## get density factor
    def density_factor(self):
        j = 0
        for i in self.data:
            j = j + (i - self.get_average()) ** 2
        # print(np.sqrt(j/ self.n))
        return np.sqrt(j / self.n)
    
    ## get the optimal fixed bandwidth
    def h_opt(self):
        defactor = self.defactor
        n = self.n
        return 1.059*np.power(defactor,1/2)*np.power(n,-1/5)
    
    ## get density function 
    def density_function(self,x):
        result = 0
        for i in self.data:
            param = (x - i) / self.hopt
            result += (self.gaussian_kernel(param)) / (self.n * self.hopt)
        # print(result)
        return result

    ## compute bandwidth
    def compute_bandwidth(self,x):
        c = (3 / 8) * np.power(np.pi,-1 / 2) * np.power(self.hopt,-5)
        m = 0.0
        # print(self.h_opt())
        for i in self.data:
            t = (x - i) / self.hopt
            m += (np.power(t,2)*self.gaussian_kernel(t)) / (self.n * self.hopt)
            # print(t)
            # print( self.gaussian_kernel(t))
        # print(m - self.density_function(x))
        g = np.power(np.abs(m - self.density_function(x)),2/5)
        h = np.power(self.hopt,9/5) * np.power(c,1/5) * np.power( self.density_function(x),1/5) / g
        return h
        
      ## compute density function
    # def compute_density(self,x):
    #     result = 0
    #     for i in self.data:
    #         t = (x - i) / self.compute_bandwidth(x)
    #         result += self.gaussian_kernel(t) / (self.n * self.compute_bandwidth(x))
    #     return result
    def compute_density(self, x):
        result = 0
        for i in self.data:
            t = (x - i) / self.hopt
            result += self.gaussian_kernel(t) / (self.n * self.hopt)
        return result


### main function
if __name__=="__main__":  
    np.random.seed(0)
    data = []
    x = np.linspace(np.array(data).min(), np.array(data).max(), 1000) 

#     x = np.linspace(-5, 10, 100)
    kde = KDE(data)
    print('best optaion',kde.hopt)
    bandwidths = []
    band = 0.0

    bandwidths = 10 ** np.linspace(-1, 1, 100)
    grid = GridSearchCV(KernelDensity(kernel='gaussian'),  {'bandwidth': bandwidths}, cv=LeaveOneOut())  
    grid.fit(np.array(data).reshape((-1,1)))
    h_good = grid.best_params_.get('bandwidth')

    #Using optimal bandwidth for kernel density estimation
    fold_kde = KernelDensity(bandwidth=grid.best_params_['bandwidth']).fit((np.array(data)).reshape(-1,1))
    density = gaussian_kde(data) 
     #Calculate the probability density of a column
    y = density(x)  #Calculate the probability density value for each point
    log_dens = fold_kde.score_samples(x[:, None])
    density_estimation = np.array([kde.compute_density(xi) for xi in x])
    line1, = plt.plot(x, density_estimation, color='r', marker='.', linestyle='-', markersize=2, alpha=0.5, linewidth=3)
    line2, = plt.plot(x, np.exp(log_dens), color='g', marker='.', linestyle='-', markersize=2, alpha=0.5, linewidth=3)
    line3, = plt.plot(x, y, color='b', marker='.', linestyle='-', markersize=2, alpha=0.5, linewidth=3)

    plt.legend(['compute-bandwidth','cross validation'])
    plt.xlabel(u'data', fontsize=14, color='r')
    plt.ylabel(u'Nuclear density', fontsize=14, color='b')
    plt.hist(data, bins=40, density=True)
    plt.show()