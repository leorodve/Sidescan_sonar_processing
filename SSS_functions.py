import numpy as np
from scipy.interpolate import CubicSpline, interp1d
from scipy.signal import medfilt, find_peaks, savgol_filter, resample
from scipy.integrate import cumtrapz
from pyproj import Proj
import seaborn as sns
from matplotlib.colors import ListedColormap
from pyxtf import xtf_read, concatenate_channel, XTFHeaderType
from tqdm import tqdm

'''
SSS functions
'''

def load_data(path, resampling):
    total_pings = 0
    pings_per_file = np.zeros((len(path)+1,1))
    for file_number in range(len(path)):
        (fh, p) = xtf_read(path[file_number])
        '''
        n_channels = fh.channel_count(verbose=True)
        for i in range(0, n_channels):
            actual_chan_info = fh.ChanInfo[i]
            if 'Port' in str(actual_chan_info.ChannelName):
                total_pings += actual_chan_info.Reserved
                pings_per_file[file_number+1] = actual_chan_info.Reserved
                break
            else:
                continue
        '''
        total_pings += len(p[XTFHeaderType.sonar])
        pings_per_file[file_number+1] = len(p[XTFHeaderType.sonar])
    max_no_pings = int(np.max(pings_per_file))
    # Creating matrixes for future use
    port_chan = np.full((resampling, max_no_pings, len(path)), -999,
                         dtype=np.float64)
    stbd_chan = np.full(port_chan.shape, -999, dtype=np.float64)
    sonar_altitude = np.full((max_no_pings, len(path)), -999, dtype=np.float64)
    slant_port = np.full(sonar_altitude.shape, -999, dtype=np.float64)
    slant_stbd = np.full(sonar_altitude.shape, -999, dtype=np.float64)
    sensor_pitch = np.full(port_chan.shape, -999, dtype=np.float64)
    Xposition = np.full(sonar_altitude.shape, -999, dtype=np.float64)
    Yposition = np.full(sonar_altitude.shape, -999, dtype=np.float64)
    sensor_heading = np.full(sonar_altitude.shape, -999, dtype=np.float64)
    sensor_speed = np.full(sonar_altitude.shape, -999, dtype=np.float64)
    ship_gyro = np.full(sonar_altitude.shape, -999, dtype=np.float64)
    layback = np.full(port_chan.shape, -999, dtype=np.float64)
    average_port = np.zeros((1, len(path)))
    average_stbd = np.zeros((1, len(path)))
    
    for file_number in range(len(path)):
        # Read file header and packets
        (fh, p) = xtf_read(path[file_number])
        
        # --------------------------- Get sonar data -------------------------- 
        np_chan1 = concatenate_channel(p[XTFHeaderType.sonar], file_header=fh, channel=0, weighted=True).astype(float)
        np_chan2 = concatenate_channel(p[XTFHeaderType.sonar], file_header=fh, channel=1, weighted=True).astype(float)
        
        upper_limit1 = np.percentile(np_chan1, 99)
        upper_limit2 = np.percentile(np_chan2, 99)
        # Clip to range (max cannot be used due to outliers)
        np_chan1.clip(0, upper_limit1 - 1, out=np_chan1)
        np_chan2.clip(0, upper_limit2 - 1, out=np_chan2)
        
        # Transpose so that the largest axis is horizontal
        # np_chan1 = np_chan1 if np_chan1.shape[0] < np_chan1.shape[1] else np_chan1.T
        # np_chan2 = np_chan2 if np_chan2.shape[0] < np_chan2.shape[1] else np_chan2.T
        np_chan1 = np_chan1.T
        np_chan2 = np_chan2.T
        
        # Resampling data and finding first backscatter value of each ping
        port_chan[:,:int(pings_per_file[file_number+1]),file_number] = abs(resample(np_chan1, resampling, axis = 0, window=2))
        stbd_chan[:,:int(pings_per_file[file_number+1]),file_number] = abs(resample(np_chan2, resampling, axis = 0, window=2))

        average_port[0, file_number] = np.mean(np_chan1)
        average_stbd[0, file_number] = np.mean(np_chan2)
        '''
        Los valores de las variables siguientes se pueden incluir en el loop
        de max_no_pings para evitar crear matrices al comienzo, ya que igual
        se itera entre columnas y files
        '''
        sonar_altitude[max_no_pings - int(pings_per_file[file_number+1]):,file_number] = [ping.SensorPrimaryAltitude for ping in p[XTFHeaderType.sonar]]
        Xposition[max_no_pings - int(pings_per_file[file_number+1]):,file_number] = [ping.ShipXcoordinate for ping in p[XTFHeaderType.sonar]]
        Yposition[max_no_pings - int(pings_per_file[file_number+1]):,file_number] = [ping.ShipYcoordinate for ping in p[XTFHeaderType.sonar]]
        ship_gyro[max_no_pings - int(pings_per_file[file_number+1]):,file_number] = [ping.ShipGyro for ping in p[XTFHeaderType.sonar]]
        sensor_speed[max_no_pings - int(pings_per_file[file_number+1]):,file_number] = [ping.SensorSpeed for ping in p[XTFHeaderType.sonar]]
        sensor_heading[max_no_pings - int(pings_per_file[file_number+1]):,file_number] = [ping.SensorHeading for ping in p[XTFHeaderType.sonar]]
        
        slant_port[:int(pings_per_file[file_number+1]),file_number] = [ping.ping_chan_headers[0].SlantRange for ping in p[XTFHeaderType.sonar]]
        slant_stbd[:int(pings_per_file[file_number+1]),file_number] = [ping.ping_chan_headers[1].SlantRange for ping in p[XTFHeaderType.sonar]]
    sonar_altitude = np.flip(sonar_altitude, axis=0)
    Xposition = np.flip(Xposition, axis=0)
    Yposition = np.flip(Yposition, axis=0)
    ship_gyro = np.flip(ship_gyro, axis=0)
    sensor_speed = np.flip(sensor_speed, axis=0)
    sensor_heading = np.flip(sensor_heading, axis=0)
    return sonar_altitude, Xposition, Yposition, ship_gyro, sensor_speed, sensor_heading, layback, slant_port, slant_stbd, pings_per_file, port_chan, stbd_chan

def stbd_and_port_data(filepath, port_channel, stbd_channel, resampling,
                       file_number, max_no_pings, pings_per_file, slant_port,
                       slant_stbd, sonar_altitude, tracking_method,
                       blanking, threshold, duration):
    #matrices for future use
    new_port_channel = np.zeros((resampling, max_no_pings, len(filepath)))
    new_stbd_channel = np.zeros((resampling, max_no_pings, len(filepath)))
    first_stbd_index = np.zeros((max_no_pings, len(filepath)))
    average_port = np.zeros((1, len(filepath)))
    average_stbd = np.zeros((1, len(filepath)))
    
    port_channel_median = np.zeros((resampling, max_no_pings, len(filepath)))
    stbd_channel_median = np.zeros((resampling, max_no_pings, len(filepath)))
    
    for file_number in range(len(filepath)):
        for i in tqdm(range(max_no_pings)):
            if i >= pings_per_file[file_number+1]:
                pass
            else:
                # The peak amplitude should be defined by the user in the future
                if tracking_method == 1:
                    #Function to automatically find the bottom
                    port_peaks, _ = find_peaks(port_channel[(resampling)//2::, i, file_number], height=500)#1.8)
                    stbd_peaks, _ = find_peaks(stbd_channel[1:resampling//2, i, file_number], height=500)#1.9)
                    first_port_index = port_peaks[-1] + (3*resampling)//4
                else:
                    #Function using user-assigned value for bottom tracking
                    start = int(blanking * resampling / slant_port[i, file_number])
                    interval = int(duration * resampling / slant_port[i, file_number])
                    port_peaks, _ = find_peaks(port_channel[resampling-start-interval:resampling-start, i, file_number], height=threshold)
                    stbd_peaks, _ = find_peaks(stbd_channel[start:start+interval, i, file_number], height=threshold)
                    first_port_index = port_peaks[-1] + resampling-start-interval
                
                first_stbd_index[i, file_number] = stbd_peaks[0]
                
                slant_range_port = (resampling - np.arange(0, int(first_port_index)+1)) * slant_port[i, file_number] / resampling
                slant_range_stbd = np.arange(int(first_stbd_index[i, file_number]), resampling) * slant_stbd[i, file_number] / resampling
                horizontal_range_port = np.sqrt(np.square(slant_range_port)-np.square(sonar_altitude[i, file_number]))
                horizontal_range_stbd = np.sqrt(np.square(slant_range_stbd)-np.square(sonar_altitude[i, file_number]))
                #Port-side
                try:
                    last_index = np.nonzero(np.nan_to_num(horizontal_range_port, nan=-999) == -999)
                    new_index = resampling - (horizontal_range_port[:last_index[0][0]] // (slant_port[i, file_number]/resampling)).astype(int)
                    if new_index[-1] == resampling:
                        new_index[-1] = resampling-1
                    new_port_channel[new_index, i, file_number] = port_channel[np.arange(0, last_index[0][0]),i, file_number]
                except IndexError:
                    new_index = resampling - (horizontal_range_port // (slant_port[i, file_number]/resampling)).astype(int)
                    if new_index[new_index.size-1] == resampling:
                        new_index[new_index.size-1] = resampling-1
                    new_port_channel[new_index, i, file_number] = port_channel[np.arange(0, int(first_port_index)+1),i, file_number]
                port_non_zeros = np.nonzero(new_port_channel[:,i, file_number] != 0)
                port_values = new_port_channel[port_non_zeros[0], i, file_number]
                first_index_port = port_non_zeros[0][0]
                last_index_port = port_non_zeros[0][port_non_zeros[0].size-1]+1
                new_port_channel[first_index_port:last_index_port,i, file_number] = interp1d(port_non_zeros[0], port_values)(np.arange(first_index_port, last_index_port))
                
                #Starboard-side
                try:
                    first_index = np.nonzero(np.nan_to_num(horizontal_range_stbd, nan=-999) == -999)
                    new_index = (horizontal_range_stbd[first_index[0][first_index[0].size-1]+1:] // (slant_stbd[i, file_number]/resampling)).astype(int)
                    if new_index[-1] == resampling:
                        new_index[-1] = resampling-1
                    new_stbd_channel[new_index, i, file_number] = stbd_channel[np.arange(first_index[0][first_index[0].size-1]+1+int(first_stbd_index[i, file_number]), resampling),i, file_number]
                except IndexError:
                    new_index = (horizontal_range_stbd // (slant_stbd[i, file_number]/resampling)).astype(int)
                    if new_index[-1] == resampling:
                        new_index[-1] = resampling-1
                    new_stbd_channel[new_index, i, file_number] = stbd_channel[np.arange(int(first_stbd_index[i, file_number]), resampling),i, file_number]
                stbd_non_zeros = np.nonzero(new_stbd_channel[:,i, file_number] != 0)
                stbd_values = new_stbd_channel[stbd_non_zeros[0], i, file_number]
                first_index_stbd = stbd_non_zeros[0][0]
                last_index_stbd = stbd_non_zeros[0][stbd_non_zeros[0].size-1]+1
                new_stbd_channel[first_index_stbd:last_index_stbd,i,file_number] = interp1d(stbd_non_zeros[0], stbd_values)(np.arange(first_index_stbd, last_index_stbd))
                
                first_value_port = new_port_channel[np.nonzero(new_port_channel[(3*resampling)//4:, i, file_number])[0][-1] + (3*resampling)//4, i, file_number]
                first_value_stbd = new_stbd_channel[np.nonzero(new_stbd_channel[:resampling//4, i, file_number])[0][0], i, file_number]
                
                port_length = resampling - (np.nonzero(new_port_channel[(3*resampling)//4:, i, file_number])[0][-1] + (3*resampling)//4)
                stbd_length = np.nonzero(new_stbd_channel[:resampling//4, i, file_number])[0][0]
                total_length = port_length + stbd_length
                
                new_port_channel[np.nonzero(new_port_channel[(3*resampling)//4:,i,file_number])[0][-1] + (3*resampling)//4 + 1:, i, file_number] = interp1d([0, total_length], [first_value_port, first_value_stbd])(np.arange(1, port_length))
                new_stbd_channel[:np.nonzero(new_stbd_channel[:resampling//4, i, file_number])[0][0], i, file_number] = interp1d([0, total_length], [first_value_port, first_value_stbd])(np.arange(port_length, total_length))
        
    return new_port_channel, new_stbd_channel, port_channel_median, stbd_channel_median

def smooth_data(x_position, y_position, max_no_pings, pings_per_file, step, file_number, gyro):
    global heading
    length = np.arange(0, int(pings_per_file[file_number+1]))
    #window size = step, polynomial order = 3
    smoothX = savgol_filter((x_position[:, file_number], length), step, 3)
    smoothX = smoothX[0,:]
    smoothY = savgol_filter((y_position[:, file_number], length), step, 3)
    smoothY = smoothY[0,:]
    
    heading = savgol_filter((gyro[:, file_number], length), step, 3)
    heading = heading[0,:]
    heading = heading.clip(min=0, max=360, out=heading)
    heading = gyro[:, file_number]
    
    factor_stbd_y = np.cos(np.radians(heading+90))
    factor_stbd_x = np.sin(np.radians(heading+90))
    factor_port_y = np.cos(np.radians(heading-90))
    factor_port_x = np.sin(np.radians(heading-90))

    return smoothX, smoothY, factor_stbd_y, factor_stbd_x, factor_port_y, factor_port_x

def bottom_track(swath, mean, N):
    first_value = np.zeros(len(swath))
    start = len(swath) // 2
    finish = len(swath) - 1
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
  
def Automatic_Gain_Control(port_data, stbd_data, window_size, target_intensity):
    data = np.append(port_data, stbd_data)
    window = np.ones(window_size)
    nonzeros = np.ones(data.shape[0])
    nonzeros[np.argwhere(data == 0)] = 0
    amplitudes = np.convolve(data, window, mode='same')
    number_nonzeros = np.convolve(nonzeros, window, mode='same')
    average = amplitudes / number_nonzeros
    sclr = np.max(data) * target_intensity * 0.01 / average
    data_corrected = sclr * data
    
    return data_corrected

def Median_filter(port_data, stbd_data, window_size):
    data = np.append(port_data, stbd_data)
    average = np.mean(data)
    median = medfilt(data, kernel_size = window_size)
    correction = median - average
    data_corrected = data - correction
    data_corrected.clip(0, np.max(data), out=data_corrected)
    return data_corrected

def Coordinate_calculator(y_value, x_value, Ymin, Xmin, Ymax, Xmax):#,df):
    y_coord = np.round(((Ymax+300)-y_value) / ((Ymax+300-(Ymin-300))/9999),0)
    x_coord = np.round(((Xmax+300)-x_value) / ((Xmax+300-(Xmin-300))/9999),0)
    return y_coord.astype(int), x_coord.astype(int)

def Lat_Long(longitud, latitud, bool):
    #If value == False, the function convert from lat/long to utm
    #If value == True, the function convert from utm to lat/long
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

def SurfRough_attribute(data, window):
    if window % 2 != 1 or window < 1:
        raise TypeError("window_size size must be a positive odd number")
    half_window = (window - 1) // 2
    data_mean = np.asarray([data[:,i] - np.mean(data[:,i]) for i in range(data.shape[1])])
    Ra = cumtrapz(abs(data_mean[:]), axis=1)
    Ra = Ra / (data.shape[0] -1)
    # Add padding at the beggining and end of the data
    firstvals = Ra[0] - np.abs(Ra[1:half_window+1][::-1] - Ra[0])
    lastvals = Ra[-1] + np.abs(Ra[-half_window-1:-1][::-1] - Ra[-1])
    data = np.concatenate((firstvals, Ra, lastvals))
    Ra_new = np.asarray([np.mean(Ra[i-half_window:i+half_window,data.shape[1]-1])
                         for i in range(half_window, data.shape[0]- half_window)])

    for i in range(data_mean.shape[0]-1):
        lower = np.argwhere(data_mean[i, :] < Ra_new[i])#[i,1022])
        data_mean[i, lower] = Ra_new[i]
    
    return data_mean
