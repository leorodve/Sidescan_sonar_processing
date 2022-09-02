
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from pyxtf import xtf_read, concatenate_channel, XTFHeaderType
from scipy.signal import find_peaks, resample
from scipy.interpolate import CubicSpline
import tkinter
from tkinter import *
from tkinter import ttk
from tkinter import filedialog


# Function to upload files selected by the user
def upload():
    global filepath
    filepath = filedialog.askopenfilename(initialdir="/", title="Select A XTF File", filetypes=(("XTF files", "*.xtf"),("all files", "*.*")))

# ----------- Window widget asking the user for file/parameters --------------
def welcome_window():
        #labels and grid
    ttk.Button(root, text="Please select a XTF file", command=upload).grid(column=2, row=1, columnspan=5)
    ttk.Label(root, text="Select a XTF file").grid(column=0, row=1)
    ttk.Button(root, text="Next", command=display_data_window).grid(column=0, row=3, columnspan=3)
    return root

# Plot the data in waterfall view
def plot_data(window):
    figure = Figure(figsize=(12, 5), dpi=100)

    ax1 = figure.add_subplot(1, 2, 1)
    ax1.imshow(port_chan, cmap='gray', aspect='auto', vmin=0, vmax= np.percentile(port_chan, 99))
    ax1.set_title('Raw Port channel data')
    ax1.set_xlabel('Along-track distance [m]')
    ax1.set_ylabel('Across-track distance [m]')
    ax1yticks = np.arange(0, port_chan.shape[0]+1, 1200)
    increase = int(75/2)
    ylimit = np.arange(0, int(slant_range)+1, increase)
    ylimit1 = np.flip(ylimit)
    ax1.set_yticks(ax1yticks)
    ax1.set_yticklabels(ylimit1)
    
    ax2 = figure.add_subplot(1, 2, 2)
    ax2.imshow(new_port_chan, extent=[0, 9557,0,1024], cmap='gray', aspect='auto', vmin=0)
    ax2.set_title('Slant-corrected Port channel data')
    ax2.set_xlabel('Along-track distance [m]')
    ax2.set_ylabel('Across-track distance [m]')
    ax2yticks = np.arange(0, new_port_chan.shape[0]+1, 512)
    ax2.set_yticks(ax2yticks)
    ax2.set_yticklabels(ylimit1)
    
    canvas = FigureCanvasTkAgg(figure, window)
    return canvas

def display_data_window():
    global np_chan1, port_chan, new_port_chan, resampling, upper_limit, slant_range
    # loop for cleaning the GUI window
    for widget in root.winfo_children():
        widget.destroy()
    
    # Read file header and packets
    (fh, p) = xtf_read(filepath)
    n_channels = fh.channel_count(verbose=True)
    actual_chan_info = [fh.ChanInfo[i] for i in range(0, n_channels)]
    print('Number of data channels: {}\n'.format(n_channels))
    
    # Print the first channel out of 4 (the first two are 'high' freq)
    print(actual_chan_info[0])
    
    # Each element in the list is a ping (XTFPingHeader)
    sonar_ch = p[XTFHeaderType.sonar]  # type: List[pyxtf.XTFPingHeader]
    
    # This retrieves the first ping in the file of the sonar type
    sonar_ch_ping1 = sonar_ch[9557]
    
    # The properties in the header defines the attributes common for all subchannels 
    # (e.g sonar often has port/stbd subchannels)
    print(sonar_ch_ping1)
    
    # Each subchannel has a XTFPingChanHeader, 
    # which contains information that can change from ping to ping in each of the subchannels
    sonar_ping1_ch_header0 = sonar_ch_ping1.ping_chan_headers[0]
    print(sonar_ping1_ch_header0)
    
    '''
    Parameters that should be given by the user
    '''
    resampling = 1024
    
    
    # --------------------------- Get sonar data -------------------------------
    upper_limit = 2 ** 16
    np_chan1 = concatenate_channel(p[XTFHeaderType.sonar], file_header=fh, channel=0, weighted=True)
    np_chan2 = concatenate_channel(p[XTFHeaderType.sonar], file_header=fh, channel=1, weighted=True)
    np_chan3 = concatenate_channel(p[XTFHeaderType.sonar], file_header=fh, channel=2, weighted=True)
    # Clip to range (max cannot be used due to outliers)
    # More robust methods are possible (through histograms / statistical outlier removal)
    np_chan1.clip(0, upper_limit - 1, out=np_chan1)
    np_chan2.clip(0, upper_limit - 1, out=np_chan2)
    np_chan3.clip(0, upper_limit - 1, out=np_chan3)

    # The sonar data is logarithmic (dB), add small value to avoid log10(0)
    np_chan1 = np.log10(np_chan1 + 1, dtype=np.float32)
    np_chan2 = np.log10(np_chan2 + 1, dtype=np.float32)
    np_chan3 = np.log10(np_chan3 + 1, dtype=np.float32)

    # Transpose so that the largest axis is horizontal
    np_chan1 = np_chan1 if np_chan1.shape[0] < np_chan1.shape[1] else np_chan1.T
    np_chan2 = np_chan2 if np_chan2.shape[0] < np_chan2.shape[1] else np_chan2.T
    np_chan3 = np_chan3 if np_chan3.shape[0] < np_chan3.shape[1] else np_chan3.T
    
    # Creating matrixes for future use
    port_altitude = np.zeros_like(np_chan1)
    first_port_index = np.zeros(np_chan1.shape[1])
    port_chan, stbd_chan = np.zeros((resampling, np_chan1.shape[1])), np.zeros((resampling, np_chan2.shape[1]))
    new_port_chan = np.zeros((resampling, np_chan1.shape[1]))
    TVG_factor = np.ones_like(np_chan1)

    '''
    Falta mejorar la precision y tiempo de computo de la funcion
    La resolucion across-along track todavia no es 100% precisa
    '''
    # Resampling data and finding first backscatter value of each ping
    for i in range(0, np_chan1.shape[1]):
        port_chan[:,i] = abs(resample(np_chan1[:,i], resampling, window=2))
        stbd_chan[:,i] = abs(resample(np_chan2[:,i], resampling, window=2))
        indexes, _ = find_peaks(port_chan[:, i], height=1.8)
        first_port_index[i] = indexes[len(indexes)-1]
        sonar_altitude = sonar_ch[port_chan.shape[1]-1-i].SensorPrimaryAltitude
        for index in range(0, int(first_port_index[i])+1):
            try:
                slant_range = (port_chan.shape[0] - index) * sonar_ch[port_chan.shape[0]-i].ping_chan_headers[0].SlantRange / port_chan.shape[0]
                if slant_range >= sonar_altitude:
                    horizontal_range = np.sqrt(np.square(slant_range)-np.square(sonar_altitude))
                    new_index = resampling - int(horizontal_range // (sonar_ch[port_chan.shape[0]-i].ping_chan_headers[0].SlantRange/resampling))
                    new_port_chan[new_index, i] = port_chan[index,i]
                    if new_port_chan[new_index-1,i] == 0 and new_port_chan[new_index-2,i] > 0:
                        new_port_chan[new_index-1,i] = (new_port_chan[new_index-2,i] + new_port_chan[new_index,i]) / 2
                else:
                    new_port_chan[resampling, i] = port_chan[index, i]
            except IndexError:
                pass
        for a in range(800,resampling):
            if new_port_chan[a, i] == 0:# and new_port_chan[a-1, i] != 0:
                for b in range(a,resampling):
                    if new_port_chan[b,i] != 0:
                        values = [new_port_chan[a-1, i], new_port_chan[b, i]]
                        distances = [a-1, b]
                        interpolator = CubicSpline(distances, values)
                        for c in range(a, b):
                            new_port_chan[c, i] = interpolator(c)
                        break

    first_plot = plot_data(root)
    first_plot.get_tk_widget().grid(row=0, column=0, columnspan=2)
    #The following plots a waterfall-view of the 100th ping (in the file)
    fig, (ax1, ax2) = plt.subplots(2,1)
    ax1.plot(np.arange(0, np_chan1.shape[0]), np_chan3[:, 196])
    ax2.plot(np.arange(0, np_chan2.shape[0]), np_chan2[:, 196])
    
    plt.show()


root = Tk()
root.title("Sonar Processing")

welcome_frame = welcome_window()
root.mainloop()
