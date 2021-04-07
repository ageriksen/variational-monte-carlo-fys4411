# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 18:22:09 2021

@author: benda
"""

import matplotlib.pyplot as plt
import numpy as np

def make_latex_table(filename, columns):
    
    # assert(len(columns) == len(column_names))
    output = "\\begin{table}[H] \n\\centering \n\\begin{tabular}{|"
    
    for i in range(len(columns)):
        output += "c|"
    output += "} \n\\hline \n"
    
    first = True
    infile = open(filename, 'r')
    for line in infile:
        if (first):
            first = False
        else:
            temp = line.split(";")
            for i in range(len(columns)-1):
                output += temp[columns[i]] + " & "
            output += temp[columns[-1]]#[:-1]
            output += " \\\\ \\hline \n" 
    
    output += "\\end{tabular} \n\\end{table}"
    
    print(output)

make_latex_table("results/results_1h.csv", [0,1,4,5,7,9])

def make_plot():
    
    infile = open("results/results_1b.csv", 'r')
    columns = [3,5]
    first = True
    energies = [] #np.zeros(12*7)
    alphas =[] # np.zeros(12*7)
    for line in infile:
        if (first):
            first = False
        else:
            temp = line.split(";")
            energies.append(float(temp[5]))
            alphas.append(float(temp[8]))
            
    length = int(len(energies)/12)
    for i in range(12):
        temp_en = np.array(energies[length*i:length*(i+1)])
        temp_al = np.array(alphas[length*i:length*(i+1)])
        ind = np.argsort(temp_al)
        plt.plot(temp_al[ind], temp_en[ind]/np.min(temp_en))
        
    plt.grid()
    plt.legend(["1;1", "1;2", "1;3", "10;1", "10;2", "10;3", "100;1", "100;2", "100;3", "500;1", "500;2", "500;3"], loc="best")
    plt.xlabel("α")  
    plt.ylabel("Scaled energy")
    plt.savefig("results/plot1b.pdf")
    plt.show()
    
 

make_plot()
    
            
                

