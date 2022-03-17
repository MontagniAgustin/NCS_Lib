###############################################################################
###############################################################################
##  
##  Nerve conduction study main program.
##  
##  Montagni Agustín
##
##
###############################################################################
###############################################################################

##--- Imports --------------------------------------------------------------###
import webplot_extractor as wp
import matplotlib.pyplot as plt
import ncs_study as ncs
import numpy as np

##--- Variables ------------------------------------------------------------###
frecs = 20000       # Frecuencia de muestreo

##--- Files ----------------------------------------------------------------###
# Remover comentario del archivo a analizar

tp,yp = wp.getSignal("Datasets\\ncs_cmap_normal_a_elbow_2mvdiv_3msdiv.txt",f_sample=frecs)
#tp,yp = wp.getSignal("Datasets\\ncs_cmap_normal_b_elbow_2mvdiv_3msdiv.txt",f_sample=frecs)
#tp,yp = wp.getSignal("Datasets\\ncs_cmap_normal_wrist_2mvdiv_3msdiv.txt",f_sample=frecs)
#tp,yp = wp.getSignal("Datasets\\ncs_snap_normal_3rd_digit_10uvdiv_1msdiv.txt",f_sample=frecs)

##--- Calculos -------------------------------------------------------------###

# Seteo de la señal y espacio lineal
ncs_sig1 = ncs.ncs_dataset()
ncs_sig1.set_Dataset(dataset=yp, timespace = tp)

# Seteo de threshold para la segmentación
ncs_sig1.set_SegmentationThreshold(ncs_sig1.dataset.max() * 0.05)

# Computación de parámetros
ncs_sig1.compute()

##--- Plot -----------------------------------------------------------------###

fig, axs = plt.subplots(3, sharex=True, sharey=True)

axs[0].plot(tp,yp)

axs[1].plot(tp,ncs_sig1.x_segm_mask * 0.005)

for i in ncs_sig1.x_segments:
    ti, yi = ncs_sig1.get_Time_Amp(i)
    axs[2].plot(ti, yi)

plt.show()
