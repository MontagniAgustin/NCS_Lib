###############################################################################
###############################################################################
##  
##  Digital Filter Definitions
##  
##  Montagni Agustín
##
##
###############################################################################
###############################################################################

##--- Imports --------------------------------------------------------------###
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt

#https://stackoverflow.com/questions/66569244/should-i-use-scipy-signal-firwin-for-amplification

"""
Devuelve los coeficientes de segundo orden para la creación de un filtro Biquad en cascada
"""
def bandpass_filter_init (lowcut, highcut, fs, order):
    nyquist = fs * 0.5
    low = lowcut / nyquist
    high = highcut / nyquist

    sos = sig.butter(N=order, Wn=[low, high], btype="band", output="sos")
    return sos

"""
Construcción del fitro pasabanda
"""
def bandpass_filter (data, lowcut, highcut, fs, order):
    sos = bandpass_filter_init(
                                lowcut=lowcut,
                                highcut=highcut,
                                fs=fs,
                                order=order
                                )
    val = sig.sosfilt(sos=sos, x=data)
    return val


###################################################################
#--- Espacio de prueba del archivo individual
###################################################################
if __name__ == "__main__":

    fs = 200000
    lowcut = 1
    highcut = 1000
    for order in [3, 6, 12]:
        sos = bandpass_filter_init(lowcut, highcut, fs, order)
        w, h = sig.sosfreqz(sos, worN=4096)
        plt.plot((fs*0.5/np.pi)*w, abs(h),
            alpha=(order+1)/13,
            label="order = %d" % order)

    plt.show()