import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import gridspec
from pyxtf import xtf_read, concatenate_channel, XTFHeaderType
from scipy.signal import resample
import tkinter
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import SSS_functions as sss
import mosaic_to_geotiff as ma

# Function to upload files selected by the user
def upload():
    global filepath
    filepath = filedialog.askopenfilenames(initialdir="/", title="Select A XTF File", filetypes=(("XTF files", "*.xtf"),("all files", "*.*")))

def welcome_window():
        #labels and grid
    ttk.Button(root, text="Please select a XTF file", command=upload).grid(column=1, row=1)
    ttk.Label(root, text="Select a XTF file").grid(column=0, row=1)
    ttk.Label(root, text="Resampling:").grid(column=0, row=2)
    ttk.Radiobutton(root, text='1024', variable=resampling_entry, value=1024).grid(column=1, row=2)
    ttk.Radiobutton(root, text='2048', variable=resampling_entry, value=2048).grid(column=2, row=2)
    ttk.Radiobutton(root, text='4096', variable=resampling_entry, value=4096).grid(column=3, row=2)
    ttk.Label(root, text="Navigation:").grid(column=0, row=3, rowspan = 2)
    ttk.Radiobutton(root, text='Use fish heading', variable=navigation_entry, value=0).grid(column=1, row=3, rowspan = 2)
    ttk.Radiobutton(root, text='Use course made good', variable=navigation_entry, value=1).grid(column=2, row=3, rowspan = 2)
    ttk.Label(root, text='(No. pings for smoothing)').grid(column=3, row=3)
    ttk.Entry(root, textvariable=course_made_good_entry).grid(column=3, row=4)
    ttk.Button(root, text="Next", command=display_data_window).grid(column=0, row=5, columnspan=4)
    return root

# Plot survey lines in mosaic
def plot_data(window):
    my_colors = sss.Color_palette()
    figure = Figure(figsize=(10, 6), dpi=100)

    ax1 = figure.subplots(1, 1)
    ax1.imshow(high_res_grid, cmap=my_colors, aspect='auto')
    ax1.set_xlabel('X distance [m]')
    ax1.set_ylabel('Y distance [m]')
    
    canvas = FigureCanvasTkAgg(figure, window)
    return canvas

# Plot the data in waterfall view
def plot_waterfall(window, stbd_min_percent, stbd_max_percent, port_min_percent, port_max_percent, first_ping, last_ping):

    my_colors = sss.Color_palette()
    file_name_with_ext = os.path.basename(filepath[0])
    file_name = os.path.splitext(file_name_with_ext)[0]
    figure = plt.figure(figsize=(12, 7), dpi=100)
    gs = gridspec.GridSpec(1, 2, wspace=0.0)
    figure.suptitle('Raw data \n File name: %s' %file_name)
    ax1 = figure.add_subplot(gs[0, 0])
    port_min_percent = 0
    port_max_percent = 99.9
    stbd_min_percent = 0
    stbd_max_percent = 99.9
    pos = ax1.imshow(port_chan[:,2400:3600,0].T, cmap=my_colors, aspect='auto', vmin=np.percentile(port_chan[:,:,0], port_min_percent), vmax=np.percentile(port_chan[:,:,0], port_max_percent))
    figure.colorbar(pos, ax=[ax1], location='left', pad=0.2)
    ax1.set_xlabel('Port-side distance [m]')
    ax1.set_ylabel('Along-track distance [ping number]')
    axyticks = np.arange(0, resampling+1, 256)
    increase = 75/4
    ylimit = np.arange(0, 75+1, increase)
    ylimit1 = np.flip(ylimit)
    ax1.set_xticks(axyticks)
    ax1.set_xticklabels(ylimit1)
    
    ax2 = figure.add_subplot(gs[0,1])
    ax2.yaxis.tick_right()
    pos1 = ax2.imshow(stbd_chan[:, 2400:3600,0].T, cmap=my_colors, aspect='auto', vmin=np.percentile(stbd_chan[:,:,0], stbd_min_percent), vmax= np.percentile(stbd_chan[:,:,0], stbd_max_percent))
    figure.colorbar(pos1, ax=[ax2], location='right', pad=0.2)
    ax2.set_xlabel('Starboard-side distance [m]')
    ax2.set_xticks(axyticks)
    ax2.set_xticklabels(ylimit)
    
    canvas = FigureCanvasTkAgg(figure, window)
    return canvas

#Function calling the widget to plot the data in waterfall view
def waterfall_window():
    for widget in root.winfo_children():
        widget.destroy()
    ma.save_mosaic(high_res_grid, Xsmooth_utm, Ysmooth_utm)
    waterfall_plot = plot_waterfall(root, 0, 99.9, 0, 99.9, 2000, 2800)
    waterfall_plot.get_tk_widget().grid(row=0, column=0, columnspan=4, rowspan=14)
    ttk.Label(root, text="port minimum display %").grid(column=5, row=0)
    port_first_entry = ttk.Entry(root, textvariable=port_min_entry)
    port_first_entry.grid(column=5, row=1)
    ttk.Label(root, text="port maximum display %").grid(column=5, row=2)
    port_second_entry = ttk.Entry(root, textvariable=port_max_entry)
    port_second_entry.grid(column=5, row=3)
    ttk.Label(root, text="starboard minimum display %").grid(column=5, row=4)
    stbd_first_entry = ttk.Entry(root, textvariable=stbd_min_entry)
    stbd_first_entry.grid(column=5, row=5)
    ttk.Label(root, text="starboard maximum display %").grid(column=5, row=6)
    stbd_second_entry = ttk.Entry(root, textvariable=stbd_max_entry)
    stbd_second_entry.grid(column=5, row=7)
    ttk.Label(root, text="First ping").grid(column=5, row=8)
    ping_first_entry = ttk.Entry(root, textvariable=first_ping_entry)
    ping_first_entry.grid(column=5, row=9)
    ttk.Label(root, text="Last ping").grid(column=5, row=10)
    ping_second_entry = ttk.Entry(root, textvariable=last_ping_entry)
    ping_second_entry.grid(column=5, row=11)
    ttk.Button(root, text="Plot again", command=display_again).grid(column=5, row=12)
    ttk.Button(root, text="Plot next file", command=display_next_file).grid(column=5, row=13)

# Display the data using user defined parameters
def display_again():
    for widget in root.winfo_children():
        widget.destroy()
    stbd_min_ent = stbd_min_entry.get() 
    stbd_max_ent = stbd_max_entry.get()
    port_max_ent = port_max_entry.get()
    port_min_ent = port_min_entry.get()
    first_ping_ent = first_ping_entry.get()
    last_ping_ent = last_ping_entry.get()
    waterfall_plot = plot_waterfall(root, stbd_min_ent, stbd_max_ent,
                           port_min_ent, port_max_ent, first_ping_ent,
                           last_ping_ent)
    waterfall_plot.get_tk_widget().grid(row=0, column=0, columnspan=4, rowspan=14)
    ttk.Label(root, text="port minimum display %").grid(column=5, row=0)
    port_first_entry = ttk.Entry(root, textvariable=port_min_entry)
    port_first_entry.grid(column=5, row=1)
    ttk.Label(root, text="port maximum display %").grid(column=5, row=2)
    port_second_entry = ttk.Entry(root, textvariable=port_max_entry)
    port_second_entry.grid(column=5, row=3)
    ttk.Label(root, text="starboard minimum display %").grid(column=5, row=4)
    stbd_first_entry = ttk.Entry(root, textvariable=stbd_min_entry)
    stbd_first_entry.grid(column=5, row=5)
    ttk.Label(root, text="starboard maximum display %").grid(column=5, row=6)
    stbd_second_entry = ttk.Entry(root, textvariable=stbd_max_entry)
    stbd_second_entry.grid(column=5, row=7)
    ttk.Label(root, text="First ping").grid(column=5, row=8)
    ping_first_entry = ttk.Entry(root, textvariable=first_ping_entry)
    ping_first_entry.grid(column=5, row=9)
    ttk.Label(root, text="Last ping").grid(column=5, row=10)
    ping_second_entry = ttk.Entry(root, textvariable=last_ping_entry)
    ping_second_entry.grid(column=5, row=11)
    ttk.Button(root, text="Plot again", command=display_again).grid(column=5, row=12)
    ttk.Button(root, text="Plot next file", command=display_next_file).grid(column=5, row=13)
    
# Display next file
def display_next_file():
    for widget in root.winfo_children():
        widget.destroy()

# Main function, for correcting and applying gains to the data, aswell as plotting
def display_data_window():
    global filepath, np_chan1, port_chan, new_port_chan, resampling, port_chan_agc, port_chan_median, stbd_chan, new_stbd_chan, stbd_chan_agc, stbd_chan_median, utmX, utmY, sonar_data, pings_per_file, Xposition, Yposition, Xsmooth, Ysmooth, xy, utmX1, utmY1, factor_stbd_x1, factor_stbd_y1, factor_port_x1, factor_port_y1, Xsmooth_utm, Ysmooth_utm, lat, long, high_res_grid
    # loop for cleaning the GUI window
    for widget in root.winfo_children():
        widget.destroy()
    resampling = resampling_entry.get()
    method_nav = navigation_entry.get()
    pings_smooth = course_made_good_entry.get()
    pings_per_file = np.zeros((len(filepath)+1,1))
    total_pings = 0
    print('Comienzo del prog')
    for file_number in range(len(filepath)):
        (fh, p) = xtf_read(filepath[file_number])
        total_pings += len(p[XTFHeaderType.sonar])
        pings_per_file[file_number+1] = len(p[XTFHeaderType.sonar])
    
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
    
    for file_number in range(len(filepath)):
        # Read file header and packets
        (fh, p) = xtf_read(filepath[file_number])
        
        # --------------------------- Get sonar data --------------------------
        np_chan1 = concatenate_channel(p[XTFHeaderType.sonar], file_header=fh, channel=0, weighted=True)
        np_chan2 = concatenate_channel(p[XTFHeaderType.sonar], file_header=fh, channel=1, weighted=True)
        
        upper_limit1 = np.percentile(np_chan1, 99)
        upper_limit2 = np.percentile(np_chan2, 99)
        # Clip to range (max cannot be used due to outliers)
        # More robust methods are possible (through histograms / statistical outlier removal)
        np_chan1.clip(0, upper_limit1 - 1, out=np_chan1)
        np_chan2.clip(0, upper_limit2 - 1, out=np_chan2)
        
        # Transpose so that the largest axis is horizontal
        np_chan1 = np_chan1 if np_chan1.shape[0] < np_chan1.shape[1] else np_chan1.T
        np_chan2 = np_chan2 if np_chan2.shape[0] < np_chan2.shape[1] else np_chan2.T
        
        # Resampling data and finding first backscatter value of each ping
        port_chan[:,:int(pings_per_file[file_number+1]),file_number] = abs(resample(np_chan1, resampling, axis = 0, window=2))
        stbd_chan[:,:int(pings_per_file[file_number+1]),file_number] = abs(resample(np_chan2, resampling, axis = 0, window=2))

        average_port[0, file_number] = np.mean(np_chan1)
        average_stbd[0, file_number] = np.mean(np_chan2)
                                               
        # Get navitation values from headers
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
    # Resample, slant correct and apply AGC and MedianFilt gains to the data
    new_port_chan, new_stbd_chan, port_chan_agc, stbd_chan_agc, port_chan_median, stbd_chan_median = sss.stbd_and_port_data(filepath, port_chan, stbd_chan, resampling, file_number, max_no_pings, pings_per_file, slant_port, slant_stbd, sonar_altitude)
    
    abc = np.arange(resampling-1,-1,-1)
    utmX, utmY, sonar_data = np.zeros((max_no_pings, 2*resampling, len(filepath))), np.zeros((max_no_pings, 2*resampling, len(filepath))), np.zeros((2*resampling, max_no_pings,len(filepath)))
    Xsmooth = np.zeros((max_no_pings, len(filepath)))
    Ysmooth = np.zeros((max_no_pings, len(filepath)))
    for file_number in range(len(filepath)):
        print('file number:', file_number)
        
        if method_nav:
            
            Xsmooth[:int(pings_per_file[file_number+1]), file_number], Ysmooth[:int(pings_per_file[file_number+1]), file_number], factor_stbd_y, factor_stbd_x, factor_port_y, factor_port_x = sss.smooth_data(Xposition, Yposition, max_no_pings, pings_per_file, filepath, pings_smooth, file_number)
            Xsmooth_utm, Ysmooth_utm = sss.Lat_Long(Xsmooth[:int(pings_per_file[file_number+1]), file_number], Ysmooth[:int(pings_per_file[file_number+1]), file_number], False)
            #Port and starboard UTM coordinates
            stbd_utmY = np.reshape(Ysmooth_utm, (int(pings_per_file[file_number+1]), 1)) + np.multiply.outer(factor_stbd_y, abc) * 75 / resampling
            stbd_utmX = np.reshape(Xsmooth_utm, (int(pings_per_file[file_number+1]), 1)) + np.multiply.outer(factor_stbd_x, abc) * 75 / resampling
            port_utmY = np.reshape(Ysmooth_utm, (int(pings_per_file[file_number+1]), 1)) + np.multiply.outer(factor_port_y, resampling-1-abc) * 75 / resampling
            port_utmX = np.reshape(Xsmooth_utm, (int(pings_per_file[file_number+1]), 1)) + np.multiply.outer(factor_port_x, resampling-1-abc) * 75 / resampling
            
            ping_utmX = np.c_[stbd_utmX, port_utmX]
            ping_utmY = np.c_[stbd_utmY, port_utmY]
            
            #Matrix with Port and Starboard coordinates of all data
            utmX[:int(pings_per_file[file_number+1]), :, file_number] = ping_utmX
            utmY[:int(pings_per_file[file_number+1]), :, file_number] = ping_utmY
        
        else:
            
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
            port_utmY = lat[:, None] + np.multiply.outer(factor_port_y, resampling-1-abc) * 75 / resampling
            port_utmX = long[:, None] + np.multiply.outer(factor_port_x, resampling-1-abc) * 75 / resampling
            
            #Matrix with Port and Starboard ping coordinates
            ping_utmX = np.c_[stbd_utmX, port_utmX]
            ping_utmY = np.c_[stbd_utmY, port_utmY]
            
            #Matrix with Port and Starboard coordinates of all data            
            utmX[:int(pings_per_file[file_number+1]), :, file_number] = ping_utmX
            utmY[:int(pings_per_file[file_number+1]), :, file_number] = ping_utmY

        #Matrix with Port and Starboard ping values
        ping_sonar_data = np.r_[np.flip(new_stbd_chan[:,:int(pings_per_file[file_number+1]), file_number], axis=0), np.flip(new_port_chan[:,:int(pings_per_file[file_number+1]), file_number], axis=0)]
        
        #Matrix with Port and Starboard values of all data
        sonar_data[:, :int(pings_per_file[file_number+1]), file_number] = ping_sonar_data
    #Generate a high resolution grid the the data in mosaic view
    high_res_grid = ma.grid_data(Xsmooth_utm, Ysmooth_utm, utmX, utmY, sonar_data, int(pings_per_file[file_number+1]), resampling)
        
    first_plot = plot_data(root)
    first_plot.get_tk_widget().grid(row=1, column=0)

    ttk.Button(root, text="Waterfall view", command=waterfall_window).grid(column=0, row=0)


# ----------- Window widget asking the user for file/parameters --------------
root = Tk()
root.title("Sonar Processing")

resampling_entry = IntVar()
resampling_entry.set(1024)
navigation_entry = IntVar()
navigation_entry.set(1)
course_made_good_entry = IntVar()
course_made_good_entry.set(300)
port_min_entry = IntVar()
port_min_entry.set("5")
port_max_entry = IntVar()
port_max_entry.set("90")
stbd_min_entry = IntVar()
stbd_min_entry.set("5")
stbd_max_entry = IntVar()
stbd_max_entry.set("90")
first_ping_entry = IntVar()
first_ping_entry.set("0")
last_ping_entry = IntVar()
last_ping_entry.set("600")

welcome_frame = welcome_window()
root.mainloop()
