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

#make_latex_table("results/results_1g.csv", [0,1,4,5,6,7,9])

def make_plot():
    
    infile = open("results/results_1b.csv", 'r')
    columns = [3,5]
    first = True
    energies = np.zeros()
    for line in infile:
        if (first):
            first = False
        else:
            temp = line.split(";")
            
                

