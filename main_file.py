import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import gridspec
import tkinter
from tkinter import ttk, font, filedialog
from tqdm import tqdm
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
    #Resampling entries
    ttk.Label(root, text="Resampling:").grid(column=0, row=2)
    ttk.Radiobutton(root, text='1024', variable=resampling_entry, value=1024).grid(column=1, row=2)
    ttk.Radiobutton(root, text='2048', variable=resampling_entry, value=2048).grid(column=2, row=2)
    ttk.Radiobutton(root, text='4096', variable=resampling_entry, value=4096).grid(column=3, row=2)
    #Navigation entries
    ttk.Label(root, text="Navigation:").grid(column=0, row=3, rowspan = 2)
    ttk.Radiobutton(root, text='Use fish heading', variable=navigation_entry, value=0).grid(column=1, row=3, rowspan = 2)
    ttk.Radiobutton(root, text='Use course made good', variable=navigation_entry, value=1).grid(column=2, row=3, rowspan = 2)
    ttk.Label(root, text='(No. pings for smoothing)').grid(column=3, row=3)
    ttk.Entry(root, textvariable=course_made_good_entry).grid(column=3, row=4)
    #Bottom tracking entries
    ttk.Label(root, text="Bottom tracking:").grid(column=0, row=5, rowspan = 3)
    ttk.Radiobutton(root, text='Use automatic bottom detection', variable=bottom_tracking_entry, value=1).grid(column=1, row=5, rowspan = 3)
    ttk.Radiobutton(root, text='Define own values:', variable=bottom_tracking_entry, value=0).grid(column=2, row=5, rowspan = 3)
    ttk.Label(root, text='Blanking').grid(column=3, row=5)
    ttk.Entry(root, textvariable=blanking_entry).grid(column=4, row=5)
    ttk.Label(root, text='Duration').grid(column=3, row=6)
    ttk.Entry(root, textvariable=duration_entry).grid(column=4, row=6)
    ttk.Label(root, text='Threshold').grid(column=3, row=7)
    ttk.Entry(root, textvariable=threshold_entry).grid(column=4, row=7)
    #Resolution entries
    ttk.Label(root, text='Mosaic resolution (m)').grid(column=0, row=8)
    ttk.Entry(root, textvariable=resolution_entry).grid(column=1, row=8)
    ttk.Button(root, text="Next", command=display_data_window).grid(column=0, row=9, columnspan=4)
    return root

# Plot survey lines in mosaic
def plot_data(window):
    my_colors = sss.Color_palette()
    figure = Figure(figsize=(10, 6), dpi=100)

    
    ax1 = figure.subplots(1, 1)
    ax1.imshow(high_res_grid, cmap=my_colors, aspect='auto')

    ax1.set_xlabel('X distance [m]')
    ticks_x = np.arange(0, high_res_grid.shape[1], int(high_res_grid.shape[1]/5))
    ticks_labels_x = np.round(np.arange(min(Xutm), max(Xutm), abs((max(Xutm)-min(Xutm))/5)),2)
    if ticks_labels_x.shape[0] < ticks_x.shape[0]:
        ticks_labels_x = np.append(ticks_labels_x, np.round(max(Xutm),2))
    ax1.set_xticks(ticks_x)
    ax1.set_xticklabels(ticks_labels_x)
    ax1.set_ylabel('Y distance [m]')
    ticks_y = np.arange(0, high_res_grid.shape[0], int(high_res_grid.shape[0]/5))
    ticks_labels_y = np.round(np.arange(min(Yutm), max(Yutm), abs((max(Yutm)-min(Yutm))/5)),2)
    if ticks_labels_y.shape[0] < ticks_y.shape[0]:
        ticks_labels_y = np.append(ticks_labels_y, np.round(max(Yutm),2))
    ax1.set_yticks(ticks_y)
    ax1.set_yticklabels(ticks_labels_y)
    
    canvas = FigureCanvasTkAgg(figure, window)
    return canvas

# Plot the data in waterfall view
def plot_waterfall(window, port_data, stbd_data, settings):

    my_colors = sss.Color_palette()
    color_values = np.arange(0, 2**16/2, 2**16/2/255)
    #color_values = np.arange(255)
    file_name_with_ext = os.path.basename(filepath[file_number])
    file_name = os.path.splitext(file_name_with_ext)[0]
    figure = plt.figure(figsize=(12, 7), dpi=100)
    gs = gridspec.GridSpec(1, 2, wspace=0.0)

    # figure.suptitle('Raw data \n File name: %s' %file_name)
    ax1 = figure.add_subplot(gs[0, 0])
    ax2 = figure.add_subplot(gs[0,1])
    
    if settings[7] == "AGC":
        print("Calculating AGC")
        data = np.asarray([sss.Automatic_Gain_Control(port_data[:,i,settings[10]], stbd_data[:,i,settings[10]], settings[8], settings[9]) for i in tqdm(range(int(pings_per_file[settings[10]+1,0])))]).T
        port_data = data[:1024, :]
        stbd_data = data[1024:, :]
        
        pos = ax1.imshow(port_data[:,settings[4]:settings[5]].T, cmap=my_colors, aspect='auto', vmin=np.percentile(color_values, settings[2]), vmax=np.percentile(color_values, settings[3]))
        pos1 = ax2.imshow(stbd_data[:, settings[4]:settings[5]].T, cmap=my_colors, aspect='auto', vmin=np.percentile(color_values, settings[0]), vmax= np.percentile(color_values, settings[1]))
        
    if settings[7] == "MedFilt":
        print("Calculating MedFilt")
        data = np.asarray([sss.Median_filter(port_data[:,i,settings[10]], stbd_data[:,i,settings[10]], settings[8]) for i in tqdm(range(int(pings_per_file[settings[10]+1,0])))]).T
        port_data = data[:1024, :]
        stbd_data = data[1024:, :]
        
        pos = ax1.imshow(port_data[:,settings[4]:settings[5]].T, cmap=my_colors, aspect='auto', vmin=np.percentile(color_values, settings[2]), vmax=np.percentile(color_values, settings[3]))
        pos1 = ax2.imshow(stbd_data[:, settings[4]:settings[5]].T, cmap=my_colors, aspect='auto', vmin=np.percentile(color_values, settings[0]), vmax= np.percentile(color_values, settings[1]))
    
    if settings[7] == None:
        pos = ax1.imshow(port_data[:,settings[4]:settings[5],settings[10]].T, cmap=my_colors, aspect='auto', vmin=np.percentile(color_values, settings[2]), vmax=np.percentile(color_values, settings[3]))
        pos1 = ax2.imshow(stbd_data[:, settings[4]:settings[5],settings[10]].T, cmap=my_colors, aspect='auto', vmin=np.percentile(color_values, settings[0]), vmax= np.percentile(color_values, settings[1]))

    figure.colorbar(pos, ax=[ax1], location='left', pad=0.2)
    # ax1.set_xlabel('Port-side distance [m]')
    ticks_x = np.arange(0, resampling+1, 256)
    increase = 75/4
    ticks_labels_x = np.arange(0, 75+1, increase)
    ticks_labels_x1 = np.flip(ticks_labels_x)
    ax1.set_xticks(ticks_x)
    ax1.set_xticklabels(ticks_labels_x1)
    ax1.spines['right'].set_visible(False)
    
    ax1.set_ylabel('Along-track distance [ping number]')
    ticks_y = np.arange(0, settings[5]-settings[4], int((settings[5]-settings[4])/6))
    ticks_labels_y = np.arange(settings[4], settings[5], int((settings[5]-settings[4])/6))
    ax1.set_yticks(ticks_y)
    ax1.set_yticklabels(ticks_labels_y)
    
    ax2.yaxis.tick_right()
    # pos1 = ax2.imshow(stbd_data[:, settings[4]:settings[5],0].T, cmap=my_colors, aspect='auto', vmin=np.percentile(stbd_data[:,:,0], settings[0]), vmax= np.percentile(stbd_data[:,:,0], settings[1]))
    figure.colorbar(pos1, ax=[ax2], location='right', pad=0.2)
    
    # ax2.set_xlabel('Starboard-side distance [m]')
    
    ax2.set_xticks(ticks_x)
    ax2.set_xticklabels(ticks_labels_x)
    ax2.set_yticks(ticks_y)
    ax2.set_yticklabels(ticks_labels_y)
    ax2.spines['left'].set_visible(False)
    
    if settings[6] == 'no':
        figure.suptitle('Raw data \n File name: %s' %file_name)
        ax1.set_xlabel('Port-side slant distance [m]')
        ax2.set_xlabel('Starboard-side slant distance [m]')
    else:
        if settings[7] == "MedFilt":
            figure.suptitle('MedFilt corrected data \n File name: %s' %file_name)
        if settings[7] == "AGC":
            figure.suptitle('AGC corrected data \n File name: %s' %file_name)
        if settings[7] == None:
            figure.suptitle('Corrected data \n File name: %s' %file_name)
        ax1.set_xlabel('Port-side horizontal distance [m]')
        ax2.set_xlabel('Starboard-side horizontal distance [m]')
    
    canvas = FigureCanvasTkAgg(figure, window)
    return canvas

def labels():
    label_font = font.Font(underline=True)
    ttk.Label(root, text="Display settings", font=label_font).grid(column=5, row=0, columnspan=2)
    ttk.Label(root, text="port minimum display %").grid(column=5, row=1, columnspan=2)
    port_first_entry = ttk.Entry(root, textvariable=port_min_entry)
    port_first_entry.grid(column=5, row=2, columnspan=2)
    ttk.Label(root, text="port maximum display %").grid(column=5, row=3, columnspan=2)
    port_second_entry = ttk.Entry(root, textvariable=port_max_entry)
    port_second_entry.grid(column=5, row=4, columnspan=2)
    ttk.Label(root, text="starboard minimum display %").grid(column=5, row=5, columnspan=2)
    stbd_first_entry = ttk.Entry(root, textvariable=stbd_min_entry)
    stbd_first_entry.grid(column=5, row=6, columnspan=2)
    ttk.Label(root, text="starboard maximum display %").grid(column=5, row=7, columnspan=2)
    stbd_second_entry = ttk.Entry(root, textvariable=stbd_max_entry)
    stbd_second_entry.grid(column=5, row=8, columnspan=2)
    ttk.Label(root, text="First ping").grid(column=5, row=9, columnspan=2)
    ping_first_entry = ttk.Entry(root, textvariable=first_ping_entry)
    ping_first_entry.grid(column=5, row=10, columnspan=2)
    ttk.Label(root, text="Last ping").grid(column=5, row=11, columnspan=2)
    ping_second_entry = ttk.Entry(root, textvariable=last_ping_entry)
    ping_second_entry.grid(column=5, row=12, columnspan=2)

#Function calling the widget to plot the data in waterfall view
def waterfall_window():
    for widget in root.winfo_children():
        widget.destroy()
        stbd_min_ent = stbd_min_entry.get() 
    stbd_max_ent = stbd_max_entry.get()
    port_max_ent = port_max_entry.get()
    port_min_ent = port_min_entry.get()
    first_ping_ent = first_ping_entry.get()
    last_ping_ent = last_ping_entry.get()
    mosaic_resolution = float(resolution_entry.get())
    ma.save_mosaic(high_res_grid, Xsmooth_utm, Ysmooth_utm, mosaic_resolution, filepath)
    display_settings = [stbd_min_ent, stbd_max_ent, port_min_ent, port_max_ent,
                        first_ping_ent, last_ping_ent, 'no', None, 0, 0, file_number]
    waterfall_plot = plot_waterfall(root, port_chan, stbd_chan, display_settings)
    waterfall_plot.get_tk_widget().grid(row=0, column=0, columnspan=4, rowspan=16)
    labels()
    ttk.Button(root, text="Slant-corrected", command=correct_slant).grid(column=5, row=13, columnspan=2)
    ttk.Button(root, text="Plot next file", command=display_next_file).grid(column=5, row=14, columnspan=2)
    ttk.Button(root, text="Save mosaic as geotiff", command=save_geotiff).grid(column=5, row=15, columnspan=2)

#Apply slant correction to the data
def correct_slant():
    for widget in root.winfo_children():
        widget.destroy()
    stbd_min_ent = stbd_min_entry.get() 
    stbd_max_ent = stbd_max_entry.get()
    port_max_ent = port_max_entry.get()
    port_min_ent = port_min_entry.get()
    first_ping_ent = first_ping_entry.get()
    last_ping_ent = last_ping_entry.get()
    display_settings = [stbd_min_ent, stbd_max_ent, port_min_ent, port_max_ent,
                        first_ping_ent, last_ping_ent, 'yes', None, 0, 0, file_number]
    waterfall_plot = plot_waterfall(root, new_port_chan, new_stbd_chan, display_settings)
    waterfall_plot.get_tk_widget().grid(row=0, column=0, columnspan=4, rowspan=22)
    labels()
    label_font = font.Font(underline=True)
    ttk.Button(root, text="Not corrected", command=waterfall_window).grid(column=5, row=13, columnspan=2)
    ttk.Button(root, text="Plot again", command=display_again).grid(column=5, row=14, columnspan=2)
    ttk.Label(root, text="Gain settings", font=label_font).grid(column=5, row=15, columnspan=2)
    ttk.Button(root, text="AGC gain", command=display_AGC).grid(column=5, row=16)
    
    ttk.Label(root, text="Window size").grid(column=5, row=17, columnspan=2)
    window_entry = ttk.Entry(root, textvariable=window_size_entry)
    window_entry.grid(column=5, row=18, columnspan=2)
    
    ttk.Label(root, text="Intensity (only for AGC)").grid(column=5, row=19, columnspan=2)
    AGC_entry = ttk.Entry(root, textvariable=intensity_entry)
    AGC_entry.grid(column=5, row=20, columnspan=2)
    
    ttk.Button(root, text="MedFilt gain", command=display_MedFilt).grid(column=6, row=16)

    ttk.Button(root, text="Plot next file", command=display_next_file).grid(column=5, row=21, columnspan=2)

#Display the data using the AGC parameters defined by the user
def display_AGC():
    for widget in root.winfo_children():
        widget.destroy()
    stbd_min_ent = stbd_min_entry.get() 
    stbd_max_ent = stbd_max_entry.get()
    port_max_ent = port_max_entry.get()
    port_min_ent = port_min_entry.get()
    first_ping_ent = first_ping_entry.get()
    last_ping_ent = last_ping_entry.get()
    window_size_ent = window_size_entry.get()
    intensity_ent = intensity_entry.get()
    display_settings = [stbd_min_ent, stbd_max_ent, port_min_ent, port_max_ent,
                        first_ping_ent, last_ping_ent, 'yes', 'AGC',
                        window_size_ent, intensity_ent, file_number]
    waterfall_plot = plot_waterfall(root, new_port_chan, new_stbd_chan, display_settings)
    waterfall_plot.get_tk_widget().grid(row=0, column=0, columnspan=4, rowspan=15)
    labels()
    ttk.Button(root, text="Without gain", command=correct_slant).grid(column=5, row=13, columnspan=2)
    ttk.Button(root, text="Plot next file", command=display_next_file).grid(column=5, row=14, columnspan=2)

#Display the data using the MedFilt parameters defined by the user
def display_MedFilt():
    for widget in root.winfo_children():
        widget.destroy()
    stbd_min_ent = stbd_min_entry.get() 
    stbd_max_ent = stbd_max_entry.get()
    port_max_ent = port_max_entry.get()
    port_min_ent = port_min_entry.get()
    first_ping_ent = first_ping_entry.get()
    last_ping_ent = last_ping_entry.get()
    window_size_ent = window_size_entry.get()
    display_settings = [stbd_min_ent, stbd_max_ent, port_min_ent, port_max_ent,
                        first_ping_ent, last_ping_ent, 'yes', 'MedFilt',
                        window_size_ent, 0, file_number]
    waterfall_plot = plot_waterfall(root, new_port_chan, new_stbd_chan, display_settings)
    waterfall_plot.get_tk_widget().grid(row=0, column=0, columnspan=4, rowspan=15)
    labels()
    ttk.Button(root, text="Without gain", command=correct_slant).grid(column=5, row=13, columnspan=2)
    ttk.Button(root, text="Plot next file", command=display_next_file).grid(column=5, row=14, columnspan=2)

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
    display_settings = [stbd_min_ent, stbd_max_ent, port_min_ent, port_max_ent,
                        first_ping_ent, last_ping_ent, 'yes', None, 0, 0, file_number]
    waterfall_plot = plot_waterfall(root, new_port_chan, new_stbd_chan, display_settings)
    waterfall_plot.get_tk_widget().grid(row=0, column=0, columnspan=4, rowspan=22)
    labels()
    label_font = font.Font(underline=True)
    ttk.Button(root, text="Not corrected", command=waterfall_window).grid(column=5, row=13, columnspan=2)
    ttk.Button(root, text="Plot again", command=display_again).grid(column=5, row=14, columnspan=2)
    ttk.Label(root, text="Gain settings", font=label_font).grid(column=5, row=15, columnspan=2)
    ttk.Button(root, text="AGC gain", command=display_AGC).grid(column=5, row=16)
    
    ttk.Label(root, text="Window size").grid(column=5, row=17, columnspan=2)
    window_entry = ttk.Entry(root, textvariable=window_size_entry)
    window_entry.grid(column=5, row=18, columnspan=2)
    
    ttk.Label(root, text="Intensity (only for AGC)").grid(column=5, row=19, columnspan=2)
    AGC_entry = ttk.Entry(root, textvariable=intensity_entry)
    AGC_entry.grid(column=5, row=20, columnspan=2)
    
    ttk.Button(root, text="MedFilt gain", command=display_MedFilt).grid(column=6, row=16)

    ttk.Button(root, text="Plot next file", command=display_next_file).grid(column=5, row=21, columnspan=2)
    
# Display next file
def display_next_file():
    global file_number
    for widget in root.winfo_children():
        widget.destroy()
    file_number += 1
    stbd_min_ent = stbd_min_entry.get() 
    stbd_max_ent = stbd_max_entry.get()
    port_max_ent = port_max_entry.get()
    port_min_ent = port_min_entry.get()
    first_ping_ent = first_ping_entry.get()
    last_ping_ent = last_ping_entry.get()
    display_settings = [stbd_min_ent, stbd_max_ent, port_min_ent, port_max_ent,
                        first_ping_ent, last_ping_ent, 'no', None, 0, 0, file_number]
    waterfall_plot = plot_waterfall(root, port_chan, stbd_chan, display_settings)
    waterfall_plot.get_tk_widget().grid(row=0, column=0, columnspan=4, rowspan=16)
    labels()
    ttk.Button(root, text="Slant-corrected", command=correct_slant).grid(column=5, row=13, columnspan=2)
    ttk.Button(root, text="Plot next file", command=display_next_file).grid(column=5, row=14, columnspan=2)
    ttk.Button(root, text="Save mosaic as geotiff", command=save_geotiff).grid(column=5, row=15, columnspan=2)

def save_geotiff():
    mosaic_resolution = float(resolution_entry.get())
    ma.save_mosaic(high_res_grid, Xsmooth_utm, Ysmooth_utm, mosaic_resolution)

# Main function, for correcting and applying gains to the data, as well as plotting
def display_data_window():
    global filepath, np_chan1, port_chan, new_port_chan, resampling, port_chan_agc, port_chan_median, stbd_chan, new_stbd_chan, stbd_chan_agc, stbd_chan_median, utmX, utmY, sonar_data, pings_per_file, Xposition, Yposition, Xsmooth, Ysmooth, xy, utmX1, utmY1, factor_stbd_x1, factor_stbd_y1, factor_port_x1, factor_port_y1, Xsmooth_utm, Ysmooth_utm, lat, long, high_res_grid, file_number
    # loop for cleaning the GUI window
    for widget in root.winfo_children():
        widget.destroy()
    resampling = resampling_entry.get()
    method_nav = navigation_entry.get()
    method_bot_track = bottom_tracking_entry.get()
    pings_smooth = course_made_good_entry.get()
    blanking = blanking_entry.get()
    threshold = threshold_entry.get()
    duration = duration_entry.get()
    mosaic_resolution = float(resolution_entry.get())
    print('Comienzo del prog')
    sonar_altitude, Xposition, Yposition, ship_gyro, sensor_speed, sensor_heading, layback, slant_port, slant_stbd, pings_per_file, port_chan, stbd_chan = sss.load_data(filepath, resampling)
    max_no_pings = int(np.max(pings_per_file))
    
    y_min_global, x_min_global, y_max_global, x_max_global = Yposition[0,0], Xposition[0,0], Yposition[0,0], Xposition[0,0]
    
    for file_number in range(len(filepath)):
        y_min = np.min(Yposition[:int(pings_per_file[file_number+1, 0]), file_number])
        x_min = np.max(Xposition[:int(pings_per_file[file_number+1, 0]), file_number])

        y_max = np.max(Yposition[:int(pings_per_file[file_number+1, 0]), file_number])
        x_max = np.min(Xposition[:int(pings_per_file[file_number+1, 0]), file_number])

        if y_min < y_min_global:
            y_min_global = y_min
        if y_max > y_max_global:
            y_max_global = y_max
        if x_min > x_min_global:
            x_min_global = x_min
        if x_max < x_max_global:
            x_max_global = x_max
    new_port_chan, new_stbd_chan, port_chan_median, stbd_chan_median = sss.stbd_and_port_data(filepath, port_chan, stbd_chan, resampling, file_number, max_no_pings, pings_per_file, slant_port, slant_stbd, sonar_altitude, method_bot_track, blanking, threshold, duration)
    
    abc = np.arange(resampling-1,-1,-1)
    utmX, utmY, sonar_data = np.zeros((max_no_pings, 2*resampling, len(filepath))), np.zeros((max_no_pings, 2*resampling, len(filepath))), np.zeros((2*resampling, max_no_pings,len(filepath)))
    Xsmooth = np.zeros((max_no_pings, len(filepath)))
    Ysmooth = np.zeros((max_no_pings, len(filepath)))
    for file_number in range(len(filepath)):
        print('file number:', file_number)
        
        if method_nav:
            
            Xsmooth[:int(pings_per_file[file_number+1, 0]), file_number], Ysmooth[:int(pings_per_file[file_number+1, 0]), file_number], factor_stbd_y, factor_stbd_x, factor_port_y, factor_port_x = sss.smooth_data(Xposition, Yposition, max_no_pings, pings_per_file, pings_smooth, file_number, ship_gyro)
            Xsmooth_utm, Ysmooth_utm = sss.Lat_Long(Xsmooth[:int(pings_per_file[file_number+1, 0]), file_number], Ysmooth[:int(pings_per_file[file_number+1, 0]), file_number], False)
            #Port and starboard UTM coordinates
            stbd_utmY = np.reshape(Ysmooth_utm, (int(pings_per_file[file_number+1, 0]), 1)) + np.multiply.outer(factor_stbd_y, abc) * 75 / resampling
            stbd_utmX = np.reshape(Xsmooth_utm, (int(pings_per_file[file_number+1, 0]), 1)) + np.multiply.outer(factor_stbd_x, abc) * 75 / resampling
            port_utmY = np.reshape(Ysmooth_utm, (int(pings_per_file[file_number+1, 0]), 1)) + np.multiply.outer(factor_port_y, resampling-1-abc) * 75 / resampling
            port_utmX = np.reshape(Xsmooth_utm, (int(pings_per_file[file_number+1, 0]), 1)) + np.multiply.outer(factor_port_x, resampling-1-abc) * 75 / resampling
            
            ping_utmX = np.c_[stbd_utmX, port_utmX]
            ping_utmY = np.c_[stbd_utmY, port_utmY]
            
            #matriz con coordenadas port y starboard de toda la data
            utmX[:int(pings_per_file[file_number+1, 0]), :, file_number] = ping_utmX
            utmY[:int(pings_per_file[file_number+1, 0]), :, file_number] = ping_utmY
        
        else:
            
            x_init, y_init = sss.Lat_Long(Xposition[int(pings_per_file[file_number+1, 0])-1, file_number], Yposition[int(pings_per_file[file_number+1, 0])-1, file_number], False)
            
            factor_y, factor_x = sss.Angle(ship_gyro[:int(pings_per_file[file_number+1, 0]), file_number])
            distance_y = sensor_speed[:int(pings_per_file[file_number+1, 0]),
                                       file_number]*0.51*(0.10000000149011612)*factor_y
            distance_x = sensor_speed[:int(pings_per_file[file_number+1, 0]),
                                       file_number]*0.51*(0.10000000149011612)*factor_x
            lat = y_init + np.flip(np.cumsum(np.flip(distance_y)))
            long = x_init + np.flip(np.cumsum(np.flip(distance_x)))
            heading_index = np.where(sensor_heading[:, file_number] > 270)
            sensor_heading[heading_index, file_number] = sensor_heading[heading_index, file_number] - 360
            factor_stbd_y, factor_stbd_x = sss.Angle(sensor_heading[:int(pings_per_file[file_number+1, 0]), file_number]+90)
            sensor_heading[heading_index, file_number] = sensor_heading[heading_index, file_number] + 360
            heading_index = np.where(sensor_heading[:, file_number] < 90)
            sensor_heading[heading_index, file_number] = sensor_heading[heading_index, file_number] + 360
            factor_port_y, factor_port_x = sss.Angle(sensor_heading[:int(pings_per_file[file_number+1, 0]), file_number]-90)
            
            #Port and starboard utm coordinates
            stbd_utmY = lat[:, None] + np.multiply.outer(factor_stbd_y, abc) * 75 / resampling
            stbd_utmX = long[:, None] + np.multiply.outer(factor_stbd_x, abc) * 75 / resampling
            port_utmY = lat[:, None] + np.multiply.outer(factor_port_y, resampling-1-abc) * 75 / resampling
            port_utmX = long[:, None] + np.multiply.outer(factor_port_x, resampling-1-abc) * 75 / resampling
            
            #matrix with ping port and starboard coordinates
            ping_utmX = np.c_[stbd_utmX, port_utmX]
            ping_utmY = np.c_[stbd_utmY, port_utmY]
            
            utmX[:int(pings_per_file[file_number+1, 0]), :, file_number] = ping_utmX
            utmY[:int(pings_per_file[file_number+1, 0]), :, file_number] = ping_utmY

        #Matrix with ping port and starboard samples
        ping_sonar_data = np.r_[np.flip(new_stbd_chan[:,:int(pings_per_file[file_number+1, 0]), file_number], axis=0), np.flip(new_port_chan[:,:int(pings_per_file[file_number+1, 0]), file_number], axis=0)]
        
        #Matrix with port and starboard samples of all data
        sonar_data[:, :int(pings_per_file[file_number+1, 0]), file_number] = ping_sonar_data

    #Creating matrix for mosaic
    utm_non_zeros = np.argwhere(utmX != 0)
    X = np.arange(np.min(utmX[utm_non_zeros[:,0], utm_non_zeros[:,1], utm_non_zeros[:,2]])-100, np.max(utmX)+100, mosaic_resolution)
    utm_non_zeros = np.argwhere(utmY != 0)
    Y = np.arange(np.min(utmY[utm_non_zeros[:,0], utm_non_zeros[:,1], utm_non_zeros[:,2]])-100, np.max(utmY)+100, mosaic_resolution)
    high_res_grid = np.zeros((len(Y),len(X)))
    for file_number in range(len(filepath)):
        high_res_grid += ma.grid_data(X, Y, utmX, utmY, sonar_data, int(pings_per_file[file_number+1, 0]), resampling, mosaic_resolution, file_number, filepath)
        high_res_grid.clip(0, sonar_data.max(), out=high_res_grid)
        
    first_plot = plot_data(root, Xsmooth_utm, Ysmooth_utm)
    first_plot.get_tk_widget().grid(row=1, column=0)

    file_number = 0
    ttk.Button(root, text="Waterfall view", command=waterfall_window).grid(column=0, row=0)


# ----------- Window widget asking the user for file/parameters --------------
root = tkinter.Tk()
root.title("Sonar Processing")

resampling_entry = tkinter.IntVar()
resampling_entry.set(1024)
navigation_entry = tkinter.IntVar()
navigation_entry.set(1)
bottom_tracking_entry = tkinter.IntVar()
bottom_tracking_entry.set(1)
course_made_good_entry = tkinter.IntVar()
course_made_good_entry.set(300)
blanking_entry = tkinter.IntVar()
blanking_entry.set(2)
duration_entry = tkinter.IntVar()
duration_entry.set(2)
threshold_entry = tkinter.IntVar()
threshold_entry.set(5)
resolution_entry = tkinter.StringVar()
resolution_entry.set("0.2")
port_min_entry = tkinter.IntVar()
port_min_entry.set("0")
port_max_entry = tkinter.IntVar()
port_max_entry.set("100")
stbd_min_entry = tkinter.IntVar()
stbd_min_entry.set("0")
stbd_max_entry = tkinter.IntVar()
stbd_max_entry.set("100")
first_ping_entry = tkinter.IntVar()
first_ping_entry.set("0")
last_ping_entry = tkinter.IntVar()
last_ping_entry.set("600")
window_size_entry = tkinter.IntVar()
window_size_entry.set("151")
intensity_entry = tkinter.IntVar()
intensity_entry.set("20")

welcome_frame = welcome_window()
root.mainloop()
