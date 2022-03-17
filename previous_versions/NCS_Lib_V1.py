import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig
import csv

def openWebplot (files):
  xypoint = [0,0]
  xyvalues = []
  xval = []
  yval = []

  with open(files, 'r') as file:
    reader = csv.reader(file,delimiter=";")
    for row in reader:
      xypoint[0] = float(row[0].replace(",","."))
      xypoint[1] = float(row[1].replace(",","."))
      xyvalues.append(xypoint.copy())   
  
  for i in range(len(xyvalues)):
    already_sorted = True
    for j in range(len(xyvalues)-i-1):
      if xyvalues[j][0] > xyvalues[j+1][0]:
        xyvalues[j], xyvalues[j + 1] = xyvalues[j + 1], xyvalues[j]
        already_sorted = False
    if already_sorted:
      break
  
  for i in range(len(xyvalues)):
    xval.append(xyvalues[i][0])
    yval.append(xyvalues[i][1])

  return xval,yval

def xyResampler (xvalues, yvalues, tstart, tend, samples):
  timesp = np.linspace(tstart,tend,samples)
  y_itp = np.interp(timesp,xvalues,yvalues)
  return timesp,y_itp

def derivateFilter (yvalues, window):
  yfilter = sig.savgol_filter(yvalues,window_length=window,polyorder=3)
  return yfilter

def zeroCross (xvalue, yvalue):
  zeroes = []
  point = []
  yvalue = np.diff(yvalue)
  for points in yvalue:
    if np.sign(yvalue[points]) != np.sign(yvalue[points+1]):
      point[0] = xvalue[points]
      point[1] = yvalue[points]
      zeroes.append(points.copy())
  return zeroes

def maxMin (xvalue, yvalue):
  ypeaks = []
  xpeaks = sig.find_peaks(x=yvalue)
  for points in xpeaks:
    ypeaks.append(yvalue[points])
  return xpeaks, ypeaks

xp,yp = openWebplot(".\Datasets\ncs_cmap_normal_a_elbow_2mvdiv_3msdiv.txt")

xp,yp = xyResampler(xp,yp,0,30e-3,1000)

yp = derivateFilter(yp, 51)

plt.plot()
plt.show()
