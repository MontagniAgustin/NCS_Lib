###############################################################################
###############################################################################
##  
##  Webplot Digitizer data resampler.
##  
##  Montagni Agustín
##
##
###############################################################################
###############################################################################

##--- Imports --------------------------------------------------------------###
import numpy as np
import scipy.signal as sig
import csv

##--- Functions ------------------------------------------------------------###

"""
  Función que recibe un archivo CSV con puntos (x,xxx;y,yyy) correspondientes a
una señal y la devuelve resampleada e interpolada a una frecuencia de muestreo
determinada.

  Formato archivo CSV
      * Datos tipo float con separador decimal ","
      * Dos columnas delimitadas por ";" 
  
      (x,xxx;y,yyy)

  En caso de especificar tiempos de comienzo y fin, así como amplitudes
  máximas y mínimas, la señal se ajustará en esa escala (TODO)  

"""
def getSignal (file, f_sample, t_start=0, t_end=0, a_min=0, a_max=0):
 
  # Puntos X Y del plot
  xval,yval = openWebplot(file)

  # Ajuste de escala
  if not t_start and not t_end:
    t_start = xval[0]
    t_end = xval[-1]
    t_span = t_end - t_start
  else:
    pass  #Agregar ajuste de escala X

  if not a_min and not a_max:
    a_min = np.amin(yval)
    a_max = np.amax(yval)
  else:
    pass  #Agregar ajuste de escala Y

  # Resampleado de la señal
  tval,yval = xyResampler(xval, yval, t_start, t_end, int(f_sample * t_span))

  # Filtrado Savigol para eliminar ruido de derivadas
  window = len(xval) * 0.1   # 10% para tener una ventana proporcional al largo de datos.
  if window%2 == 0:
      window += 1
  yval = derivateFilter(yval, int(window))
  
  return tval, yval

""""
Abre un archivo csv delimitado por ; con separador decimal . con valores xy dispersados
Devuelve dos listas con los valores x e y ordenados de menor a mayor x
"""
def openWebplot (files):
  xypoint = [0,0]
  xyvalues = []
  xval = []
  yval = []

  # Apertura de archivo y carga de vectores
  with open(files, 'r') as file:
    reader = csv.reader(file,delimiter=";")
    for row in reader:
      xypoint[0] = float(row[0].replace(",","."))
      xypoint[1] = float(row[1].replace(",","."))
      xyvalues.append(xypoint.copy())   
  
  # Bubble sort con secuencia de escape para el ordenado de puntos en el eje X
  for i in range(len(xyvalues)):
    already_sorted = True
    for j in range(len(xyvalues)-i-1):
      if xyvalues[j][0] > xyvalues[j+1][0]:
        xyvalues[j], xyvalues[j + 1] = xyvalues[j + 1], xyvalues[j]
        already_sorted = False
    if already_sorted:
      break
  
  # Formato de vectores de salida
  for i in range(len(xyvalues)):
    xval.append(xyvalues[i][0])
    yval.append(xyvalues[i][1])

  return xval,yval

"""
Función que recibe vectores con valores X Y y una frecuencia de resampleo.
Devuelve dos vectores con los valores de entrada resampleados a la frecuencia indicada.
"""
def xyResampler (xvalues, yvalues, tstart, tend, samples):
  timesp = np.linspace(tstart,tend,samples)
  y_itp = np.interp(timesp,xvalues,yvalues)
  return timesp,y_itp

"""
Funcion que recive un vector de valores.
Devuelve un vector de valores filtrados en sus derivadas, eliminando ruido en zonas
de alta pendiente.
"""
def derivateFilter (yvalues, window):
  if len(yvalues) > window:
    if window %2 == 0:
      window = window - 1
    yfilter = sig.savgol_filter(yvalues,window_length=window,polyorder=3)
  else:
    window = len(yvalues)
    if window %2 == 0:
      window = window - 1
      
    yfilter = sig.savgol_filter(yvalues,window_length=window,polyorder=window-1)
  return yfilter
