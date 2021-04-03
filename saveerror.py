#!/usr/bin/env python3

"""
bootstrap run results (sans histogram):
=========================
Runtime: 1.00759 sec
Bootstrap statistics:
Original set:       0.4880171004966685
mean statistics:    0.4880226165911547
bias:               5.5160944861798455e-06
std. Error:         0.005676068973146812
=========================
original set already exists, mean of sampled means maybe?, 
bias probably worth saving. std.error definitely should save. 

blocking run results (last few results):
    .
    .
    .
    standard error = 1.23387e-05
**************************************************
Warning: Data size: 89989, not a power of 2.
Truncating data to 65536
standard error = 0.00121674
**************************************************
Warning: Data size: 89989, not a power of 2.
Truncating data to 65536
standard error = 2.62663e-06
**************************************************
Warning: Data size: 899899, not a power of 2.
Truncating data to 524288
standard error = 2.83644e-06
**************************************************
Warning: Data size: 899899, not a power of 2.
Truncating data to 524288
standard error = 5.57943e-07
**************************************************
Warning: Data size: 899899, not a power of 2.
Truncating data to 524288
/home/ms/uni/fys4411/project1/src/blocking.py:47: RuntimeWarning: invalid value encountered in true_divide
  ( (gamma/s)**2 * 2**arange(1, d+1)[::-1] )[::-1] ) )[::-1]
Warning: Use more data
standard error = 0
===
The std.error is the only output here. it should definitely be saved. should the bias and others
be attempted found? 
"""

#first off, save to file filename
def SaveError(filename, data):
    #import csv
    #outfile = open(filename, 'w', newline='')
    #csv_writer = csv.writer(outfile, delimiter=";")
    #for dat in data:
    #    csv_writer.writerow(dat)

    outfile = open(filename, 'w')
    for dat in data:
        outfile.write("%g\n"%dat)


if __name__ == '__main__':
    import numpy as np
    filename = "testfile.txt"
    data = np.array([1, 3, 2, 2, 4, 6, 3, 5, 3, 5])

    SaveError(filename, data)






