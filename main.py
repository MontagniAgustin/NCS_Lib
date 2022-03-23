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
import matplotlib as mp
import ncs_study as ncs
import numpy as np

##--- Variables ------------------------------------------------------------###
frecs = 20000       # Frecuencia de muestreo

##--- Files ----------------------------------------------------------------###
# Remover comentario del archivo a analizar

#tp,yp = wp.getSignal("Datasets\\ncs_cmap_normal_a_elbow_2mvdiv_3msdiv.txt",f_sample=frecs)
#tp,yp = wp.getSignal("Datasets\\ncs_cmap_normal_b_elbow_2mvdiv_3msdiv.txt",f_sample=frecs)
#tp,yp = wp.getSignal("Datasets\\ncs_cmap_normal_wrist_2mvdiv_3msdiv.txt",f_sample=frecs)
tp,yp = wp.getSignal("Datasets\\ncs_snap_normal_3rd_digit_10uvdiv_1msdiv.txt",f_sample=frecs)

##--- Calculos -------------------------------------------------------------###

# Seteo de la señal y espacio lineal
ncs_sig1 = ncs.ncs_dataset()
ncs_sig1.set_Dataset(dataset=yp, timespace = tp)

# Seteo de threshold para la segmentación
ncs_sig1.set_SegmentationThreshold(ncs_sig1.dataset.max() * 0.2,
                                    np.amax(np.diff(ncs_sig1.dataset)) * 0.11)

# Computación de parámetros
ncs_sig1.compute()

##--- Plot -----------------------------------------------------------------###

fig, ax = plt.subplots()

ax.set_title("Estudio de conducción nerviosa")
ax.set_xlabel("Tiempo")
ax.set_ylabel("Volts")
ax.plot(tp,yp, label="Resampleada")

#for i in ncs_sig1.segmentos:
#    ta, ya = ncs_sig1.get_Time_Amp(i)
#
#    ax.plot(ta,ya, label="Resampleada")

ta, ya = ncs_sig1.get_Time_Amp(ncs_sig1.x_segments[1])
ax.plot(ta,ya, label="Resampleada")
print("Latencia")
print(ncs_sig1.t_onset_lat)

#ax.axhline(ncs_sig1.y_vmed)


#ta, ya = ncs_sig1.get_Time_Amp(ncs_sig1.x_mainsegment)
#ax.plot(ta,ya, label="Resampleada")
#ti, yi = ncs_sig1.get_Time_Amp(ncs_sig1.x_stimulus)
#ax.plot(ti,yi)

#ax.plot(tp,yp, label="Resampleada")
#ax.plot(tp,ncs_sig1.x_segm_mask * 0.05, label="Mascara")

"""
fig, axs = plt.subplots(3, sharex=True, sharey=True)

axs[0].plot(tp,yp)

axs[1].plot(tp,ncs_sig1.x_segm_mask * 0.005)

ta, ya = ncs_sig1.get_Time_Amp(ncs_sig1.stimulus(ncs_sig1.segmentos))
axs[2].plot(ta, ya)

print("Latencia Onset")
print(ncs_sig1.t_onset_lat)
#for i in ncs_sig1.segmentos:
#    ti, yi = ncs_sig1.get_Time_Amp(i)
#    axs[1].plot(ti, yi)



#for i in ncs_sig1.x_segments:
#    ti, yi = ncs_sig1.get_Time_Amp(i)
#    axs[2].plot(ti, yi)

#axs[2].axhline(ncs_sig1.y_baselines[1])
"""


plt.show()
