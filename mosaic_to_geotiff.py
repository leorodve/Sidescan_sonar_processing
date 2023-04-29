import numpy as np
import matplotlib.pyplot as plt
from bisect import bisect_left
from rasterio.transform import Affine
import rasterio
from scipy.interpolate import interp1d
import cv2
import SSS_functions as sss

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


def grid_data(X_smooth_utm, Y_smooth_utm, X_utm, Y_utm, data, no_pings, resampling):
    res_x = 0.2
    res_y = 0.2
    X = np.arange(np.min(X_smooth_utm)-100, np.max(X_smooth_utm)+100, res_x)
    Y = np.arange(np.min(Y_smooth_utm)-100, np.max(Y_smooth_utm)+100, res_y)
    mosaic = np.zeros((len(Y),len(X)))

    # position the data in the correct coordinate of the new grid
    for a in range(0,no_pings):
        print(a)
        id_x = take_closest(X, X_utm[a,:,0], resampling)
        id_y = len(Y) - take_closest(Y, Y_utm[a,:,0], resampling)
        mosaic[id_y, id_x] = data[:,a,0]
    
    # interpolate values in empty pixels
    for b in range(0, mosaic.shape[1]):
        print(b)
        non_zeros = np.argwhere(mosaic[:,b] != 0)
        try:
            first = non_zeros[0][0]
            last = non_zeros[len(non_zeros)-1][0]+1
            mosaic[first:last,b] = interp1d(np.reshape(non_zeros, len(non_zeros)),
                  np.reshape(mosaic[non_zeros, b],
                              len(mosaic[non_zeros, b])))(np.arange(first, last))
        except:
            IndexError
    return mosaic

def save_mosaic(mosaic, x_utm, y_utm):
    #save it first as png to create the rgb matrix
    plt.imsave('mosaic_png.png', mosaic, cmap=sss.Color_palette())
    
    img = cv2.imread('mosaic_png.png')  #cv2 read the images as bgr
    img_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

    res_x = 0.2
    res_y = 0.2
    X = np.arange(np.min(x_utm)-100, np.max(x_utm)+100, res_x)
    Y = np.arange(np.min(y_utm)-100, np.max(y_utm)+100, res_y)
    transform = Affine.translation(np.min(X) - res_x / 2, np.min(Y) - res_y / 2) * Affine.scale(res_x, res_y)
    #reshape img to save as geotiff
    #goes from [height,width,bands] to [bands, height, width]
    image = np.moveaxis(img_rgb.squeeze(),-1,0)
    new_dataset = rasterio.open(
        'new_mosaic.tif',
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