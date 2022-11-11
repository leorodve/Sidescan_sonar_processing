import numpy as np
from scipy.interpolate import CubicSpline, interp1d
from scipy.signal import medfilt, find_peaks
import seaborn as sns
from matplotlib.colors import ListedColormap
import pandas as pd

'''
SSS functions
'''

def stbd_and_port_data(filepath, port_channel, stbd_channel, resampling,
                       file_number, max_no_pings, pings_per_file, slant_port,
                       slant_stbd, sonar_altitude):
    #matrices for future use
    new_port_channel = np.zeros((resampling, max_no_pings, len(filepath)))
    new_stbd_channel = np.zeros((resampling, max_no_pings, len(filepath)))
    first_port_index = np.zeros((max_no_pings, len(filepath)))
    first_stbd_index = np.zeros((max_no_pings, len(filepath)))
    average_port = np.zeros((1, len(filepath)))
    average_stbd = np.zeros((1, len(filepath)))
    
    port_channel_agc = np.zeros((resampling, max_no_pings, len(filepath)))
    stbd_channel_agc = np.zeros((resampling, max_no_pings, len(filepath)))
    port_channel_median = np.zeros((resampling, max_no_pings, len(filepath)))
    stbd_channel_median = np.zeros((resampling, max_no_pings, len(filepath)))
    
    
    for file_number in range(len(filepath)):
        for i in range(0, max_no_pings):
            print('iteration of max no pings', i)
            if i >= pings_per_file[file_number+1]:
                pass
            else:
                # The peak amplitude should be defined by the user in the future
                port_peaks, _ = find_peaks(port_channel[(3*resampling)//4::, i, file_number], height=1.8)
                first_port_index[i, file_number] = port_peaks[len(port_peaks)-1] + (3*resampling)//4 + 1
                    
                stbd_peaks, _ = find_peaks(stbd_channel[1:resampling//4, i, file_number], height=1.9)
                first_stbd_index[i, file_number] = stbd_peaks[0]
                
                slant_range_port = (resampling - np.arange(0, int(first_port_index[i, file_number])+1)) * slant_port[i, file_number] / resampling
                slant_range_stbd = np.arange(int(first_stbd_index[i, file_number]), resampling) * slant_stbd[i, file_number] / resampling
                horizontal_range_port = np.sqrt(np.square(slant_range_port)-np.square(sonar_altitude[i, file_number]))
                horizontal_range_stbd = np.sqrt(np.square(slant_range_stbd)-np.square(sonar_altitude[i, file_number]))
                
                #Port-side
                try:
                    last_index = np.nonzero(np.nan_to_num(horizontal_range_port, nan=-999) == -999)
                    new_index = resampling - (horizontal_range_port[:last_index[0][0]] // (slant_port[i, file_number]/resampling)).astype(int)
                    if new_index[new_index.size-1] == resampling:
                        new_index[new_index.size-1] = resampling-1
                    new_port_channel[new_index, i, file_number] = port_channel[np.arange(0, last_index[0][0]),i, file_number]
                except IndexError:
                    new_index = resampling - (horizontal_range_port // (slant_port[i, file_number]/resampling)).astype(int)
                    if new_index[new_index.size-1] == resampling:
                        new_index[new_index.size-1] = resampling-1
                    new_port_channel[new_index, i, file_number] = port_channel[np.arange(0, int(first_port_index[i, file_number])+1),i, file_number]
                port_non_zeros = np.nonzero(new_port_channel[:,i, file_number] != 0)
                port_values = new_port_channel[port_non_zeros[0], i, file_number]
                first_index_port = port_non_zeros[0][0]
                last_index_port = port_non_zeros[0][port_non_zeros[0].size-1]+1
                new_port_channel[first_index_port:last_index_port,i, file_number] = interp1d(port_non_zeros[0], port_values)(np.arange(first_index_port, last_index_port))
                
                #Starboard-side
                try:
                    first_index = np.nonzero(np.nan_to_num(horizontal_range_stbd, nan=-999) == -999)
                    new_index = (horizontal_range_stbd[first_index[0][first_index[0].size-1]+1:] // (slant_stbd[i, file_number]/resampling)).astype(int)
                    if new_index[new_index.size-1] == resampling:
                        new_index[new_index.size-1] = resampling-1
                    new_stbd_channel[new_index, i, file_number] = stbd_channel[np.arange(first_index[0][first_index[0].size-1]+1+int(first_stbd_index[i, file_number]), resampling),i, file_number]
                except IndexError:
                    new_index = (horizontal_range_stbd // (slant_stbd[i, file_number]/resampling)).astype(int)
                    if new_index[new_index.size-1] == resampling:
                        new_index[new_index.size-1] = resampling-1
                    new_stbd_channel[new_index, i, file_number] = stbd_channel[np.arange(int(first_stbd_index[i, file_number]), resampling),i, file_number]
                stbd_non_zeros = np.nonzero(new_stbd_channel[:,i, file_number] != 0)
                stbd_values = new_stbd_channel[stbd_non_zeros[0], i, file_number]
                first_index_stbd = stbd_non_zeros[0][0]
                last_index_stbd = stbd_non_zeros[0][stbd_non_zeros[0].size-1]+1
                new_stbd_channel[first_index_stbd:last_index_stbd,i,file_number] = interp1d(stbd_non_zeros[0], stbd_values)(np.arange(first_index_stbd, last_index_stbd))
                
                #Gains
                port_channel_agc[:int(pings_per_file[file_number+1]), i, file_number] = Automatic_Gain_Control(new_port_channel[:int(pings_per_file[file_number+1]),i, file_number], 51, average_port[0, file_number], first_index_port, last_index_port)
                stbd_channel_agc[:int(pings_per_file[file_number+1]), i, file_number] = Automatic_Gain_Control(new_stbd_channel[:int(pings_per_file[file_number+1]),i, file_number], 51, average_stbd[0, file_number], first_index_stbd, last_index_stbd)
                port_channel_median[first_index_port:last_index_port, i, file_number] = sss.Median_filter(new_port_channel[first_index_port:last_index_port, i, file_number], 201, average_port1[0, file_number])
                stbd_channel_median[first_index_stbd:last_index_stbd, i, file_number] = sss.Median_filter(new_stbd_channel[first_index_stbd:last_index_stbd, i, file_number], 201, average_stbd1[0, file_number])
    return new_port_channel, new_stbd_channel, port_channel_agc, stbd_channel_agc, port_channel_median, stbd_channel_median

def bottom_track(swath, mean, N):
    first_value, new_port_channel = np.zeros(len(swath)), np.zeros(len(swath))
    start = len(swath) // 2
    finish = len(swath) - 1
    n = N // 2
    average = mean / 2
    for i in range(start, len(swath)):
        try:
            if swath[i] >= average and swath[i]-swath[i-1]<=0.2 and abs(swath[i]-swath[i+1])>=0.2:
                first_value[i] = 1 # Set row value equal to 1 to identify ocean floor
        except IndexError:
            pass
    while finish >= start:
        if first_value[finish] == 1:
            value = swath[finish] # First Backscatter value of ocean floor
            break
        else:
            finish -= 1
    return first_value, value, finish

def moving_window(data, window_length, window_step):
    global a
    if window_step == 0:
        nrows = 1
    else:
        nrows = ((data.size-window_length)//window_step)+1
    n = data.strides[0]
    return np.lib.stride_tricks.as_strided(data, shape=(nrows,
                                                        window_length),
                                                strides=(window_step*n,n))

def Automatic_Gain_Control(data, window_size, average, first, last):
    correction = np.zeros(len(data))
    step = int(round(abs((data.size / window_size - 1) * window_size)))
    analysis_window = moving_window(data[first:last], window_size, step)
    for i in enumerate(analysis_window):
        local_average = np.mean(i[1])
        correction[first + i[0]*(window_size//2):first + window_size*(i[0]+1)] = data[first + i[0]*(window_size//2):first + window_size*(i[0]+1)] * average / local_average
    last_index = first +  (len(analysis_window)-1) * (len(analysis_window[0])-1)//2
    if last_index < last:
        local_average = np.mean(data[last_index:last+1])
        correction[last_index:last+1] = data[last_index:last+1] * average / local_average
    return correction

def Median_filter(data, window_size, average):
    median = medfilt(data, kernel_size = window_size)
    new_data = np.zeros_like(data)
    factor = average / median
    new_data = data * factor
    new_data = np.nan_to_num(new_data)
    new_data.clip(0, np.max(data), out=new_data)
    return new_data

def Coordinate_calculator(y_value, x_value, Ymin, Xmin, Ymax, Xmax):#,df):
    y_coord = np.round(((Ymax+300)-y_value) / ((Ymax+300-(Ymin-300))/9999),0)
    x_coord = np.round(((Xmax+300)-x_value) / ((Xmax+300-(Xmin-300))/9999),0)
    return y_coord.astype(int), x_coord.astype(int)

def Lat_Long(longitud, latitud, bool):
    #If bool == False, the function transform lat/long to utm
    #If bool == True, the function transform utm to lat/long
    p = Proj(proj='utm',zone=10,ellps='WGS84', preserve_units=False, inverse=False)
    if bool == False:
        utmX, utmY = p(longitud, latitud)
    else:
        utmX, utmY = p(longitud, latitud, inverse=True)
    return utmX, utmY

def Angle(heading):
    factor_y, factor_x = np.zeros(len(heading)), np.zeros(len(heading))
    for i, heading_value in enumerate(heading):
        if heading_value > 90 and heading_value <= 180:
            angle = 180 - heading_value
            factor_y[i] = -np.cos(np.radians(angle))
            factor_x[i] = np.sin(np.radians(angle))
        if heading_value >180 and heading_value <= 270:
            angle = 270 - heading_value
            factor_y[i] = -np.sin(np.radians(angle))
            factor_x[i] = -np.cos(np.radians(angle))
        if heading_value >270 and heading_value <= 360:
            angle = 360 - heading_value
            factor_y[i] = np.cos(np.radians(angle))
            factor_x[i] = -np.sin(np.radians(angle))
        else:
            factor_y[i] = np.cos(np.radians(heading_value))
            factor_x[i] = np.sin(np.radians(heading_value))
    return factor_y, factor_x

def Color_palette():
    mstlbronze = ['#000000','#000300', '#000700', '#000B00', '#2B0F00',
                  '#2B1200', '#451600', '#451900', '#571D00', '#572000',
                  '#652300', '#652600', '#702900', '#702C00', '#7A2F00',
                  '#7A3200', '#833500', '#833700', '#8A3A00', '#8A3C00',
                  '#913F00', '#914100', '#974300', '#974600', '#9C4800',
                  '#9C4A00', '#A14C00', '#A14E00', '#A65000', '#A65200',
                  '#AA5400', '#AA5600', '#AE5800', '#AE5900', '#B25B00',
                  '#B25D00', '#B65E00', '#B66000', '#B96100', '#B96300',
                  '#BC6400', '#BC6600', '#BF6700', '#BF6800', '#C26A00',
                  '#C26B00', '#C56C00', '#C56D00', '#C86F00', '#C87000',
                  '#CA7100', '#CA7200', '#CD7300', '#CD7400', '#CF7500',
                  '#CF7600', '#D17700', '#D17800', '#D47900', '#D47A00',
                  '#D67B00', '#D67C00', '#D87D00', '#D87E00', '#DA8000',
                  '#DA8100', '#DC8200', '#DC8300', '#DE8400', '#DE8500',
                  '#DF8600', '#DF8700', '#E18800', '#E18900', '#E38A00',
                  '#E38B00', '#E58C00', '#E58D00', '#E68E00', '#E68F00',
                  '#E89100', '#E89200', '#E99300', '#E99400', '#EB9500',
                  '#EB9700', '#EC9800', '#EC9900', '#EE9B00', '#EE9C00',
                  '#EF9E00', '#EF9F00', '#F1A100', '#F1A200', '#F2A400',
                  '#F2A600', '#EFA800', '#F0A900', '#F1AB00', '#F2AD00',
                  '#F3B100', '#F4B300', '#F5B500', '#F6B700', '#F7B900',
                  '#F7BC00', '#F8BE00', '#F9C000', '#F9C300', '#FAC500',
                  '#FAC800', '#FBCB00', '#FBCD00', '#FCD000', '#FCD300',
                  '#FCD600', '#FDD900', '#FDDF00' ,'#FEE200', '#FEE600',
                  '#FEE900', '#FEED00', '#FEF000', '#FEF400', '#FEF800',
                  '#FEFC00']
    new_color_palette = ListedColormap(sns.color_palette(mstlbronze))
    return new_color_palette
