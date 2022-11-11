
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from pyxtf import xtf_read, concatenate_channel, XTFHeaderType
from scipy.signal import resample
import tkinter
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import SSS_functions as sss


# Function to upload files selected by the user
def upload():
    global filepath
    filepath = filedialog.askopenfilenames(initialdir="/", title="Select A XTF File", filetypes=(("XTF files", "*.xtf"),("all files", "*.*")))

# ----------- Window widget asking the user for file/parameters --------------
def welcome_window():
        #labels and grid
    ttk.Button(root, text="Please select a XTF file", command=upload).grid(column=2, row=1, columnspan=5)
    ttk.Label(root, text="Select a XTF file").grid(column=0, row=1)
    ttk.Button(root, text="Next", command=display_data_window).grid(column=0, row=3, columnspan=3)
    return root

# Plot the data in waterfall view
def plot_data(window, stbd_min_percent, stbd_max_percent, port_min_percent, port_max_percent):
    my_colors = sss.Color_palette()
    figure = Figure(figsize=(12, 7), dpi=100)
    
    ax1, ax2 = figure.subplots((2, 1), sharex=True)
    
    pos = ax1.imshow(new_port_chan[:,:,0], cmap=my_colors, aspect='auto', vmin=np.percentile(new_port_chan[:,:,0], port_min_percent), vmax=np.percentile(new_port_chan[:,:,0],port_max_percent))
    figure.colorbar(pos, ax=ax1)
    ax1.set_title('Slant-corrected Data - Raw')
    ax1.set_ylabel('Port-side distance [m]')
    axyticks = np.arange(0, np_chan1.shape[0]+1, 256)
    increase = 75/4
    ylimit = np.arange(0, 75+1, increase)
    ylimit1 = np.flip(ylimit)
    ax1.set_yticks(axyticks)
    ax1.set_yticklabels(ylimit)
    
    pos1 = ax2.imshow(new_stbd_chan[:, :, 0], cmap=my_colors, aspect='auto', vmin=np.percentile(new_stbd_chan[:,:,0], stbd_min_percent), vmax= np.percentile(new_stbd_chan[:,:,0], stbd_max_percent))
    figure.colorbar(pos1, ax=ax2)
    ax2.set_xlabel('Along-track distance [ping number]')
    ax2.set_ylabel('Starboard-side distance [m]')
    ax2.set_yticks(axyticks)
    ax2.set_yticklabels(ylimit1)
    
    figure.subplots_adjust(hspace=0)
    
    canvas = FigureCanvasTkAgg(figure, window)
    return canvas

def new_display_window(): #function for live-adjusting gain parameters on the waterfall-view
    # loop for cleaning the GUI window
    for widget in root.winfo_children():
        widget.destroy()
    stbd_min_ent = stbd_min_entry.get() 
    stbd_max_ent = stbd_max_entry.get()
    port_max_ent = port_max_entry.get()
    port_min_ent = port_min_entry.get()
    first_plot = plot_data(root, stbd_min_ent, stbd_max_ent,
                           port_min_ent, port_max_ent)
    first_plot.get_tk_widget().grid(row=0, column=0, columnspan=4, rowspan=8)
    ttk.Label(root, text="% port min").grid(column=5, row=0)
    port_first_entry = ttk.Entry(root, textvariable=port_min_entry)
    port_first_entry.grid(column=5, row=1, columnspan=1)
    ttk.Label(root, text="% port max").grid(column=5, row=2)
    port_second_entry = ttk.Entry(root, textvariable=port_max_entry)
    port_second_entry.grid(column=5, row=3, columnspan=1)
    ttk.Label(root, text="% starboard min").grid(column=5, row=4)
    stbd_first_entry = ttk.Entry(root, textvariable=stbd_min_entry)
    stbd_first_entry.grid(column=5, row=5, columnspan=1)
    ttk.Label(root, text="% starboard max").grid(column=5, row=6)
    stbd_second_entry = ttk.Entry(root, textvariable=stbd_max_entry)
    stbd_second_entry.grid(column=5, row=7, columnspan=1)
    ttk.Button(root, text="Plot again", command=new_display_window).grid(column=5, row=8, columnspan=1)

def display_data_window():
    global np_chan1, port_chan, new_port_chan, resampling, port_chan_agc, port_chan_median, stbd_chan, new_stbd_chan, stbd_chan_agc, stbd_chan_median, utmX, utmY, sonar_data
    # loop for cleaning the GUI window
    for widget in root.winfo_children():
        widget.destroy()
    pings_per_file = np.zeros((len(filepath)+1,1))
    total_pings = 0
    #read the total number of pings per file and save it on a vector for future use
    for file_number in range(len(filepath)):
        (fh, p) = xtf_read(filepath[file_number])
        total_pings += len(p[XTFHeaderType.sonar])
        pings_per_file[file_number+1] = len(p[XTFHeaderType.sonar])  
    '''
    Parameters that should be given by the user
    '''
    resampling = 1024
    
    max_no_pings = int(np.max(pings_per_file))
    # Creating matrixes for future use
    port_chan = np.full((resampling, max_no_pings, len(filepath)), -999,
                         dtype=np.float64)
    stbd_chan = np.full((resampling, max_no_pings, len(filepath)), -999,
                         dtype=np.float64)
    sonar_altitude = np.full((max_no_pings, len(filepath)), -999, dtype=np.float64)
    slant_port = np.full((max_no_pings, len(filepath)), -999, dtype=np.float64)
    slant_stbd = np.full((max_no_pings, len(filepath)), -999, dtype=np.float64)
    sensor_pitch = np.full((resampling, max_no_pings, len(filepath)), -999,
                              dtype=np.float64)
    Xposition = np.full((max_no_pings, len(filepath)), -999, dtype=np.float64)
    Yposition = np.full((max_no_pings, len(filepath)), -999, dtype=np.float64)
    sensor_heading = np.full((max_no_pings, len(filepath)), -999, dtype=np.float64)
    sensor_speed = np.full((max_no_pings, len(filepath)), -999, dtype=np.float64)
    ship_gyro = np.full((max_no_pings, len(filepath)), -999, dtype=np.float64)
    layback = np.full((resampling, max_no_pings, len(filepath)), -999,
                              dtype=np.float64)
    average_port = np.zeros((1, len(filepath)))
    average_stbd = np.zeros((1, len(filepath)))
    new_port_chan = np.zeros((resampling, max_no_pings, len(filepath)))
    new_stbd_chan = np.zeros((resampling, max_no_pings, len(filepath)))
    
    port_chan_agc = np.zeros((resampling, max_no_pings, len(filepath)))
    stbd_chan_agc = np.zeros((resampling, max_no_pings, len(filepath)))
    port_chan_median = np.zeros((resampling, max_no_pings, len(filepath)))
    stbd_chan_median = np.zeros((resampling, max_no_pings, len(filepath)))
    #Retrieve the sonar data from each file
    for file_number in range(len(filepath)):
        # Read file header and packets
        (fh, p) = xtf_read(filepath[file_number])
        
        # --------------------------- Get sonar data --------------------------
        upper_limit = 2 ** 16
        np_chan1 = concatenate_channel(p[XTFHeaderType.sonar], file_header=fh, channel=0, weighted=True)
        np_chan2 = concatenate_channel(p[XTFHeaderType.sonar], file_header=fh, channel=1, weighted=True)
        
        # Clip to range (max cannot be used due to outliers)
        # More robust methods are possible (through histograms / statistical outlier removal)
        np_chan1.clip(0, upper_limit - 1, out=np_chan1)
        np_chan2.clip(0, upper_limit - 1, out=np_chan2)
        
        # The sonar data is logarithmic (dB), add small value to avoid log10(0)
        np_chan1 = np.log10(np_chan1 + 1, dtype=np.float32)
        np_chan2 = np.log10(np_chan2 + 1, dtype=np.float32)
        
        # Transpose so that the largest axis is horizontal
        np_chan1 = np_chan1 if np_chan1.shape[0] < np_chan1.shape[1] else np_chan1.T
        np_chan2 = np_chan2 if np_chan2.shape[0] < np_chan2.shape[1] else np_chan2.T
        
        # Resampling data and finding first backscatter value of each ping
        port_chan[:,:int(pings_per_file[file_number+1]),file_number] = abs(resample(np_chan1, resampling, axis = 0, window=2))
        stbd_chan[:,:int(pings_per_file[file_number+1]),file_number] = abs(resample(np_chan2, resampling, axis = 0, window=2))

        average_port[0, file_number] = np.mean(np_chan1)
        average_stbd[0, file_number] = np.mean(np_chan2)

        max_no_pings - int(pings_per_file[file_number+1])
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
    
    y_min_global, x_min_global, y_max_global, x_max_global = Yposition[0,0], Xposition[0,0], Yposition[0,0], Xposition[0,0]
    #set the minimum/maximum coordinate values
    for i in range(0, len(filepath)):
        y_min = np.min(Yposition[:int(pings_per_file[file_number+1]), file_number])
        x_min = np.max(Xposition[:int(pings_per_file[file_number+1]), file_number])
        y_max = np.max(Yposition[:int(pings_per_file[file_number+1]), file_number])
        x_max = np.min(Xposition[:int(pings_per_file[file_number+1]), file_number])
        if y_min < y_min_global:
            y_min_global = y_min
        if y_max > y_max_global:
            y_max_global = y_max
        if x_min > x_min_global:
            x_min_global = x_min
        if x_max < x_max_global:
            x_max_global = x_max
    #Get the slant-corrected, AGC and Median gain starboard & port data
    new_port_chan, new_stbd_chan, port_chan_agc, stbd_chan_agc, port_chan_median, stbd_chan_median = sss.stbd_and_port_data(filepath, port_chan, stbd_chan, resampling, file_number, max_no_pings, pings_per_file, slant_port, slant_stbd, sonar_altitude)
    abc = np.arange(resampling-1,-1,-1)
    #Get the coordinates (in utm) of every sample value (not only NADIR)
    utmX, utmY, sonar_data = np.zeros((1, 2 * resampling)), np.zeros((1, 2 * resampling)), np.zeros((2 * resampling, 1))
    for file_number in range(len(filepath)):
        print('file number:', file_number)
        x_init, y_init = sss.Lat_Long(Xposition[int(pings_per_file[file_number+1])-1, file_number], Yposition[int(pings_per_file[file_number+1])-1, file_number], False)
        
        factor_y, factor_x = sss.Angle(ship_gyro[:int(pings_per_file[file_number+1]), file_number])
        distance_y = sensor_speed[:int(pings_per_file[file_number+1]),
                                   file_number]*0.51*(0.10000000149011612)*factor_y
        distance_x = sensor_speed[:int(pings_per_file[file_number+1]),
                                   file_number]*0.51*(0.10000000149011612)*factor_x
        lat = y_init + np.flip(np.cumsum(np.flip(distance_y)))
        long = x_init + np.flip(np.cumsum(np.flip(distance_x)))
        heading_index = np.where(sensor_heading[:, file_number] > 270)
        sensor_heading[heading_index, file_number] = sensor_heading[heading_index, file_number] - 360
        factor_stbd_y, factor_stbd_x = sss.Angle(sensor_heading[:int(pings_per_file[file_number+1]), file_number]+90)
        sensor_heading[heading_index, file_number] = sensor_heading[heading_index, file_number] + 360
        heading_index = np.where(sensor_heading[:, file_number] < 90)
        sensor_heading[heading_index, file_number] = sensor_heading[heading_index, file_number] + 360
        factor_port_y, factor_port_x = sss.Angle(sensor_heading[:int(pings_per_file[file_number+1]), file_number]-90)
        # Port and Starboard UTM coordinates
        stbd_utmY = lat[:, None] + np.multiply.outer(factor_stbd_y, abc) * 75 / resampling
        stbd_utmX = long[:, None] + np.multiply.outer(factor_stbd_x, abc) * 75 / resampling
        port_utmY = lat[:, None] + np.multiply.outer(factor_port_y, 1023-abc) * 75 / resampling
        port_utmX = long[:, None] + np.multiply.outer(factor_port_x, 1023-abc) * 75 / resampling
        
        #Matrix with port and starboard coordinates of the ping
        ping_utmX = np.c_[stbd_utmX, port_utmX]
        ping_utmY = np.c_[stbd_utmY, port_utmY]
        #Matrix with port and starboard the sample values of the ping
        ping_sonar_data = np.r_[new_port_chan[:,:int(pings_per_file[file_number+1]), file_number], new_stbd_chan[:,:int(pings_per_file[file_number+1]), file_number]]
        #Matrix with port and starboard coordinates of all pings
        utmX = np.r_[utmX, ping_utmX]
        utmY = np.r_[utmY, ping_utmY]
        #Matrix with port and starboard the sample values of all pings
        sonar_data = np.c_[sonar_data, ping_sonar_data]
        # Min and Max coordinates from lat and long to utm
        x_min_utm, y_max_utm = sss.Lat_Long(x_max, y_max, False)
        x_max_utm, y_min_utm = sss.Lat_Long(x_min, y_min, False)
        
    #GUI
    first_plot = plot_data(root, 3.5, 98, 4, 90)
    first_plot.get_tk_widget().grid(row=0, column=0, columnspan=4, rowspan=8)
    ttk.Label(root, text="% port min").grid(column=5, row=0)
    port_first_entry = ttk.Entry(root, textvariable=port_min_entry)
    port_first_entry.grid(column=5, row=1, columnspan=1)
    ttk.Label(root, text="% port max").grid(column=5, row=2)
    port_second_entry = ttk.Entry(root, textvariable=port_max_entry)
    port_second_entry.grid(column=5, row=3, columnspan=1)
    ttk.Label(root, text="% starboard min").grid(column=5, row=4)
    stbd_first_entry = ttk.Entry(root, textvariable=stbd_min_entry)
    stbd_first_entry.grid(column=5, row=5, columnspan=1)
    ttk.Label(root, text="% starboard max").grid(column=5, row=6)
    stbd_second_entry = ttk.Entry(root, textvariable=stbd_max_entry)
    stbd_second_entry.grid(column=5, row=7, columnspan=1)
    ttk.Button(root, text="Plot again", command=new_display_window).grid(column=5, row=8, columnspan=1)

# ----------- Window widget asking the user for file/parameters --------------
root = Tk()
root.title("Sonar Processing")

port_min_entry = IntVar()
port_min_entry.set("5")
port_max_entry = IntVar()
port_max_entry.set("90")
stbd_min_entry = IntVar()
stbd_min_entry.set("5")
stbd_max_entry = IntVar()
stbd_max_entry.set("90")

welcome_frame = welcome_window()
root.mainloop()
