import numpy as np
import matplotlib.pyplot as plt
from bisect import bisect_left
from rasterio.transform import Affine
import rasterio
from scipy.interpolate import interp1d
import SSS_functions as sss
from tqdm import tqdm
import os

def take_closest(myList, myNumber, length):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = [bisect_left(myList, myNumber[i]) for i in range(0,2*length)]
    position = np.zeros_like(pos)
    for i in range(0,2*length):
        if pos[i] == 0:
            position[i] = 0
        if pos[i] == len(myList):
            position[i] = len(myList)-1
        else:
            before = myList[pos[i]-1]
            after = myList[pos[i]]
            if after - myNumber[i] < myNumber[i] - before:
                position[i] = pos[i]
            else:
                position[i] = pos[i]-1
    return position

def right_pos(myList, pos):
    if pos == 0:
        return pos
    if pos == len(myList):
        return pos-1
    before = myList[pos-1]
    after = myList[pos]
    if after - pos < pos - before:
        return pos
    else:
        return pos-1

def grid_data(X_smooth_utm, Y_smooth_utm, X_utm, Y_utm, data, no_pings, resampling, resolution):
    X = np.arange(np.min(X_smooth_utm)-100, np.max(X_smooth_utm)+100, resolution)
    Y = np.arange(np.min(Y_smooth_utm)-100, np.max(Y_smooth_utm)+100, resolution)
    mosaic = np.zeros((len(Y),len(X)))

    # position the data in the correct coordinate of the new grid
    print("Gridding data")
    for a in tqdm(range(no_pings)):
        id_x = take_closest(X, X_utm[a,:,0], resampling)
        id_y = len(Y) - take_closest(Y, Y_utm[a,:,0], resampling)
        mosaic[id_y, id_x] = data[:,a]
        
    for b in tqdm(range(mosaic.shape[1])):
        non_zeros = np.argwhere(mosaic[:,b] != 0)
        try:
            first = non_zeros[0][0]
            last = non_zeros[len(non_zeros)-1][0]+1
            mosaic[first:last,b] = interp1d(np.reshape(non_zeros, len(non_zeros)),
                  np.reshape(mosaic[non_zeros, b],
                              len(mosaic[non_zeros, b])))(np.arange(first, last))
        except:
            IndexError
    mosaic = np.nan_to_num(mosaic)
    return mosaic

def scalar_to_rgb(data):
    #Transform scalar-data to four dimensional RGBA data; the fourth is alpha
    #First we must normalize the data
    norm = plt.Normalize(vmin=data.min(), vmax=data.max())
    #Map the scalar-data using colormap function from matplotlib
    #The default cmap used by matplotlib is 'viridis', but that can be changed
    #using different cmaps in the function below
    # cm = plt.cm.viridis
    #or a cmap defined by the user
    cm = sss.Color_palette()
    #cm = plt.cm.gist_rainbow
    #cm = plt.cm.viridis

    rgb_data = cm(norm(data)) #four dimensional matrix
    """
    It is important to note that once the data is saved using a certain cmap,
    it CAN'T be changed when plotting. 
    """
    #Now we re-escale the data from 0.0-1.0 float scale, to 0-255 integer scale
    rgb_data = (rgb_data[:,:,:3]*255).astype('uint8')
    # rgb_data = rgb_dat

def save_mosaic(mosaic, x_utm, y_utm, resolution, filepath):
    img_rgb = scalar_to_rgb(mosaic)
    X = np.arange(np.min(x_utm)-100, np.max(x_utm)+100, resolution)
    Y = np.arange(np.min(y_utm)-100, np.max(y_utm)+100, resolution)
    transform = Affine.translation(np.min(X) - resolution / 2, np.min(Y) - resolution / 2) * Affine.scale(resolution, resolution)
    #reshape img to save as geotiff
    #goes from [height,width,bands] to [bands, height, width]
    image = np.moveaxis(img_rgb.squeeze(),-1,0)
    file_name_with_ext = os.path.basename(filepath[0])
    file_name = os.path.splitext(file_name_with_ext)[0]
    new_dataset = rasterio.open(
        '%s.tif' %file_name,
        'w',
        driver='GTiff',
        height=mosaic.shape[0],
        width=mosaic.shape[1],
        count=3,
        dtype=img_rgb.dtype,
        crs='EPSG:32610',
        transform=transform,
    )
    new_dataset.write(image)
    new_dataset.close()
