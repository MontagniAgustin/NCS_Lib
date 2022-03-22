###############################################################################
###############################################################################
##  
##  Nerve conduction study feature extraction.
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
import scipy.integrate as integ

##--- Classes --------------------------------------------------------------###

"""
Clase que contiene la información de una señal de estudio de conducción nerviosa
y que realiza el procesamiento necesario para obtener sus parametros.
"""
class ncs_dataset():

    ###################################################################
    #--- Init
    ###################################################################  
    def __init__(self):
        pass

    ###################################################################
    #--- Setters
    ###################################################################
    """
    Setea la base de datos de la señal y su frecuencia de muestreo
    """
    def set_Dataset(self, dataset, f_sample = 0, timespace = 0):
        self.dataset = np.array(dataset)
        self.inv_dataset = np.negative(self.dataset)

        if len(timespace) <= 1:
            self.samplefrec = f_sample
            self.timespace = np.linspace(0,
                                         len(self.dataset)/f_sample,
                                         len(self.dataset)
                                         )
        elif not f_sample:
            self.timespace = timespace
            self.samplefrec = 1 / (timespace[1] - timespace[0])
        else:
            pass # Ver que pasa si no se ingresa ni frecuencia de sampleo ni espacio lineal (TODO)

    """
    Establece un Threshold para el calculo de la segmentación
    """
    def set_SegmentationThreshold(self, y_thld, d_thld):
        self.segthld  = y_thld
        self.dthld = d_thld
    
    ###################################################################
    #--- Getters
    ###################################################################
    
    """
    Devuelve los valores del espacio lineal y la amplitud 
    correspondientes a los indices ingresados
    """
    def get_Time_Amp (self, xpoints):
        try:
            t_val = []
            y_val = []
            for point in xpoints:
                y_val.append(self.dataset[point].copy())
                t_val.append(self.timespace[point].copy())
            return np.array(t_val), np.array(y_val)
        except:
            return np.nan, np.nan

    def get_Amp (self, xpoints):
        y_val = []
        for point in xpoints:
            y_val.append(self.dataset[point].copy())
        return np.array(y_val)
    
    """"
    Devuelve los valores del espacio lineal correspondientes a un valor 
    de amplitud
    """
    def get_TimeLine (self, ypoints):
        t_val = []
        y_val = []

        for point in range(len(self.dataset)):
            t_val.append(self.timespace[point].copy())
            y_val.append(ypoints)
        return t_val, y_val

    """
    Devuelve el tiempo entre dos muestras
    """
    def get_Duration(self,start, end):
        return self.timespace[end] - self.timespace[start]
    
    ###################################################################
    #--- Functions
    ###################################################################
    """
    Genera los cálculos necesarios para la extracción de caracteristicas
    """
    def compute (self):

        # Extracción de parámetros y segmentación
        #self.x_max = self.max_values()
        #self.x_min = self.min_values()
        #self.x_slope = self.slopes()
        self.x_segments, self.x_segm_mask, self.x_segments_div = self.segmentation(self.segthld, self.dthld)

        # 
        self.y_baselines = self.baselines()
        self.t_latencies = self.latencies()
        self.x_mainsegment = self.mainsegment(self.segmentos)
        self.x_stimulus = self.stimulus(self.segmentos)
        self.t_onset_lat = self.onset_latency()
        self.val_mainarea = self.mainarea()
        self.v_area = self.mainarea()
        pass
        #

    """
    Devuelve la duración de cada segmento
    """
    def latencies(self):
        durations = []
        for segment in self.x_segments:
            durations.append(self.get_Duration(segment[0],segment[-1]))
        return durations

    """
    Devuelve segmentos de la señal basado en un threshold sobre la 
    media de la señal
    """
    def segmentation (self, y_thld, d_thld):
        y_mean = np.mean(self.dataset)
        y_mean_band = [y_mean - y_thld, y_mean + y_thld]

        greater_thld = np.argwhere(self.dataset > y_mean)
        greater_band = np.argwhere(self.dataset > y_mean_band[1])
        greater_deriv = np.argwhere(np.absolute(np.diff(self.dataset)) > d_thld)

        greater = np.union1d(np.intersect1d(greater_thld,greater_deriv), greater_band)
        
        lower_thld = np.argwhere(self.dataset < y_mean)
        lower_band = np.argwhere(self.dataset < y_mean_band[0])
        lower_deriv = np.argwhere(np.absolute(np.diff(self.dataset)) > d_thld)

        lower = np.union1d(np.intersect1d(lower_thld, lower_deriv), lower_band)

        segment_mask = np.zeros(len(self.dataset))
        np.put(segment_mask, lower, -1)
        np.put(segment_mask, greater, 1)
        
        segments_active = np.where(np.diff(segment_mask) > 0)[0]
        segments = np.where(np.diff(segment_mask))[0]
        
        segments_divided = np.split(range(len(self.dataset)),segments+1)
        segments_active_divided = np.split(range(len(self.dataset)),segments_active+1)
        
        active = []
        for seg in segments_active_divided:
            if segment_mask[seg[0]:seg[-1]].all() != 0:
                active.append(seg.copy().astype(int))
        
        self.segmentos = active

        return [segments_divided, segment_mask, segments_active_divided]

    """
    Devuelve el segmento principal
    """
    def mainsegment (self, segments):
        seg_len = [len(i) for i in segments]
        ind = np.argmax(seg_len)
        return segments[ind]

    """
    Devuelve el estímulo si existiera
    """
    def stimulus (self, segments):
        if len(segments) == 1:
            return np.nan
        else:
            seg_len = [len(i) for i in segments]
            ind = np.argmin(seg_len)
            return segments[ind]

    """
    Devuelve la latencia entre el estímulo o comienzo de la señal y el
    segmento principal
    """
    def onset_latency (self):
        if not np.isnan(self.x_stimulus):
            t_lat = self.get_Duration(self.x_stimulus[-1], self.x_mainsegment[0])
        else:
            t_lat = self.get_Duration(0, self.x_mainsegment[0])
        return t_lat

    def mainarea (self):
        if not np.isnan(self.x_stimulus):
            amps = self.dataset[self.x_stimulus[-1] : self.x_mainsegment[0]]
            v_med = np.mean(amps)
            area = integ.simps(np.where(amps >= v_med))
        else:
            baseline = self.dataset[0 : self.x_mainsegment[0]]
            v_med = np.mean(baseline)
            amps = self.get_Amp(self.x_mainsegment)

            segment = np.where(amps > v_med)
            values = self.get_Amp(segment)

            v_med_arr = np.full(len(segment), v_med)

            area = integ.simps(values) - integ.simps(v_med_arr)
            
        print("Area positiva")
        print(area)
        return area

        

    """
    Devuelve valores medios de los segmentos que no superan el threshold
    """
    def baselines (self):
        means = []
        index = np.where(np.diff(self.x_segm_mask))[0]
        seg_mask = np.split(self.x_segm_mask, index +1)
        for base, mask in zip(self.x_segments, seg_mask):
            if mask.all() == 0:
                means.append(np.mean(self.get_Time_Amp(base)[1]))
        return means        

    ###########################################################################
    """
    Devuelve los puntos de la señal donde cambia el signo de
    la pendiente (Sin utilizar)
    """
    def slopes (self):
        slope_section = []
        slopex = np.sign(np.diff(self.dataset))
        for point in range(len(slopex)):
            if (slopex[point] * slopex[point-1]) == -1:
                slope_section.append(point)
        return slope_section

    """
    Devuelve los puntos de la señal donde se cruza con el nivel
    ingresado. (Sin utilizar)
    """
    def crossing (self, level):
        x_val = np.where(
                np.diff(
                np.sign([value - level for value in self.dataset])
                ))[0]
        return x_val

    """
    Devuelve los puntos de la señal donde se detecta un máximo (Sin utilizar)
    """
    def max_values(self):
        x_val = sig.find_peaks(self.dataset)[0] 
        return x_val
        
    """
    Devuelve los puntos de la señal donde se detecta un mínimo (Sin utilizar)
    """
    def min_values(self):
        x_val = sig.find_peaks(self.inv_dataset)[0]
        return x_val

    """
    Devuelve los anchos de los distintos máximos (Sin utilizar)
    """
    def peak_widths(self):
        self.peak_width = sig.peak_widths(self.dataset, self.xmaxval, rel_height=1)

###################################################################
#--- Espacio de prueba del archivo individual
###################################################################
if __name__ == "__main__":
    """
    x = np.linspace(0, 6 * np.pi, 1000)
    x = np.sin(x) + 0.6 * np.sin(2.6 * x)

    peaks, _ = sig.find_peaks(x)
    results_half = sig.peak_widths(x, peaks, rel_height=0.5)
    results_half[0]  # widths

    results_full = sig.peak_widths(x, peaks, rel_height=1)
    results_full[0]  # widths
    plt.plot(x)
    plt.plot(peaks, x[peaks], "x")
    plt.hlines(*results_half[1:], color="C2")
    plt.hlines(*results_full[1:], color="C3")
    plt.show()

    """

    


