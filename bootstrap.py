#!/usr/bin/env python3

import os 

import csv
import numpy as np
from time import time
from scipy.stats import norm
import matplotlib.pyplot as plt
from numpy.random import randint

#from saveerror import SaveError

#where to save data
#DATA_ID = "results_1c"
DATA_ID = "results"

def data_path(data_id): 
    return os.path.join(DATA_ID, data_id)

def SaveError(filename, data):
    """
    bootstrap needs to store more per row, so simplest to 
    implement save function locally.
    """
    outfile = open(filename, 'w', newline='')
    csv_writer = csv.writer(outfile, delimiter=";")
    for dat in data:
        csv_writer.writerow(dat)

    #outfile = open(filename, 'w')
    #for dat in data:
    #    outfile.write("%g\n"%dat)

def tsboot(data, statistic, R, l):
    t = np.zeros(R)
    n = len(data)
    k = int(np.ceil(float(n)/l))
    inds = np.arange(n)
    t0 = time()

    #time series bootstrap
    for i in range(R):
        #bootstrap sample from k cunks of data
        #chunksize l
        _data = np.concatenate([data[j:j+l] for j in randint(0, n-1, k)])[0:n]
        t[i] = statistic(_data)

    #analysis
    print(  "="*25  )
    print(  "Runtime: %g sec" % (time()-t0))
    print(  "Bootstrap statistics: ")
    print(  "Original set:      ", statistic(data)  )
    print(  "mean statistics:   ", np.mean(t)   )
    print(  "bias:              ", np.mean(t) - statistic(data) )
    print(  "std. Error:        ", np.std(t)    )
    print(  "="*25  )
    return t

def stat(data):
    return np.mean(data)

def csvError(infile, delim, plot):
    csv_reader = csv.reader(infile, delimiter=delim)
    error = []
    for row in csv_reader:
        row = list(filter(None, row))
        X = np.array(row, dtype=float)
        #t = tsboot(X, stat, int(len(X)*.5), int(len(X)*.2))
        t = tsboot(X, stat, 2**12, 2**10)
        error.append([np.std(t), np.mean(t) - stat(X)])
        #TODO justify numbers used for number of and size of 
        #       bootstraps
        """
        histogram of bootstrapped data:
        norm=1/True is deprecated and threw an error.
        updating with 'density' solved it.
        """
        if plot:
            n, binsboot, patches = plt.hist(t, 50, density=1, facecolor='red', alpha=.75)
            """
            'best fit'
            mpl.normpdf is deprecated/removed. scipy steps up
            """
            y = norm.pdf(binsboot, np.mean(t), np.std(t))
            lt = plt.plot(binsboot, y, 'r--', linewidth=1)
            plt.xlabel('mean energy')
            plt.ylabel('Probability')
            plt.grid(True)
            plt.show()

        #break
        """
        The program goes per row in the energies file. I'm unsure which variables go
        to which line and running over each takes a little while. I've also not found
        a great measure of the number of bootstraps and the amount that ought be sampled
        per bootstrap, but this seems to give a pretty good result so far. 
        """
    return error


#read in data
infile = open(data_path("energies_1c.csv"), "r")
error = csvError(infile, ";", True)
infile.close()
SaveError(data_path("error_bootstrap_1c.dat"), error)

infile = open(data_path("energies_1d.csv"), "r")
error = csvError(infile, ";", False)
infile.close()
SaveError(data_path("error_bootstrap_1d.dat"), error)
"""
cat results/results_1c.txt:
    n Particles;n Dimensions;n Metropolis Steps;Equilibration Fraction;Accepted Steps;Found Energy;Elapsed Time;n Parameters;Parameters (undefinedNumber)
    1;1;10000;0.1;0.73986;0.50353;0.0229235;2;0.46;0.001
"""
