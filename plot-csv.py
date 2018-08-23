"""
@author: Seth

for plotting .csv convergence files
"""

import matplotlib.pyplot as plt
import csv
from math import sqrt

def get_csv_data(filename, col):
    data = []
    times = []
    file = open(filename, 'rU')
    reader = csv.reader(file)
    get_data = False
    for row in reader:
        if get_data:
            times.append(float(row[0]))
            data.append(float(row[col]))
        elif len(row) > 0:
            if row[0] == 'time':
                get_data = True
    file.close()
    return [times, data]

def get_csv_data_skip(filename, col, skip):
    data = []
    times = []
    file = open(filename, 'rU')
    reader = csv.reader(file)
    get_data = 0
    for row in reader:
        get_data += 1
        if get_data > skip:
            times.append(float(row[0]))
            data.append(float(row[col]))
    file.close()
    return [times, data]

def plot_data(data, names):
    for i in range(len(data)):
        plt.plot(data[i][0], data[i][1], label=names[i])
    plt.xlabel('time')
    plt.ylabel('Q(t)')
    plt.legend(loc='best')
    return

def plot_data(data, names):
    for i in range(len(data)):
        plt.plot(data[i][0], data[i][1], label=names[i])
    plt.xlabel('time')
    plt.ylabel('Q(t)')
    plt.legend(loc='best')
    return

""" input: full filename (with extension) as string """
def get_i1data(filename, i1=1):
    file = open(filename, 'rU')
    reader = csv.reader(file)
    data = [row[i1] for row in reader]
    file.close()
    return data

def get_i1i2data(filename, i1=1, i2=3):
    file = open(filename, 'rU')
    reader = csv.reader(file)
    data = [[row[i1], row[i2]] for row in reader]
    file.close()
    return data

""" ([[str,str]]) data = [[m(t), m_{r>r0}(t)] for all t] """
def get_mi(data, i0=1):
    return float(data[i0][0]) - float(data[len(data)/2][1])

def get_mf(data):
    return float(data[len(data)-1][0]) - float(data[len(data)-1][1])

"""d0 = val(4h), d1 = val(2h), d2 = val(h)"""
def get_err(d0, d1, d2):
    e01 = (d0 - d1) / 12.0
    e12 = (d1 - d2) / 3.0
    return max(abs(e01), abs(e12))

def make_mfname(fname, resn):
    return "mass-" + str(resn) + "-" + fname + ".csv"

""" give fname as just the outfile variable from the simulation
    and for fout give full file name or '' to prevent writing 
    returns [mass(t(i0)), mass_err(t(i0)), mass(tf), mass_err(tf)]"""
def mass_err(fname, resn0=8, fout='', i0=1):
    filename0 = make_mfname(fname, resn0)
    m0 = get_i1data(filename0, 1)
    
    filename1 = make_mfname(fname, 2*resn0)
    m1 = get_i1i2data(filename1, 1, 3)
        
    filename2 = make_mfname(fname, 4*resn0)
    m2 = get_i1data(filename2, 1)

    mi = get_mi(m1, i0)
    ei = get_err(float(m0[i0]), float(m1[i0][0]), float(m2[i0]))
    mf = get_mf(m1)
    k = len(m1) - 1
    ef = get_err(float(m0[k]), float(m1[k][0]), float(m2[k]))
    mlist = [mi, ei, mf, ef]

    if fout:
        times = get_i1data(filename0, 0)
        fwr = open(fout, 'w')
        writer = csv.writer(fwr)
        writer.writerow(["time","(4h-2h)/12","(2h-h)/3","max(|err|)"])
        for i in range(i0, len(m2)):
            wrow = [float(times[i])]
            wrow.append((float(m0[i]) - float(m1[i][0])) / 12.0)
            wrow.append((float(m1[i][0]) - float(m2[i])) / 3.0)
            wrow.append(max(abs(wrow[1]), abs(wrow[2])))
            writer.writerow(wrow)
        fwr.close()
    
    return mlist

def absorbed(mlist):
    mi = mlist[0]
    ei = mlist[1]
    mf = mlist[2]
    ef = mlist[3]
    absorb = 1 - (mf / mi)
    lob = 1 - ((mf + ef) / (mi - ei))
    elo = absorb - lob
    upb = 1 - ((mf - ef) / (mi + ei))
    eup = upb - absorb
    return [absorb, elo, eup]
    
def merrplot(fnames, dsqs, resn0=8, fout='', i0=1):
    x = [sqrt(float(dsq)) for dsq in dsqs]
    y = []
    errlo = []
    errup = []
    for fname in fnames:
        data = absorbed(mass_err(fname, resn0, fout, i0))
        y.append(100*data[0])
        errlo.append(400*data[1])
        errup.append(400*data[2])
    plt.errorbar(x, y, yerr=[errlo, errup], lw=1, elinewidth=2, capsize=2)
    plt.xlabel('initial pulse width / black hole mass')
    plt.ylabel('mass absorption percent')
    plt.title('scalar field incident on static black hole:\n' \
              'mass absorption vs. pulse width')
    return "enter plt.show() to view plot"
