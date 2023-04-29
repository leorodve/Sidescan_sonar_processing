import cv2
import os
import numpy as np

# Read an image. User must define image path below
mosaic = cv2.imread('line3_trimmed_res02mt_nadirfilt_egn_STARBOARD.tif')

# User must define tile size below (number of pixels)
tile_size = 100

tile_row = 'row_'
tile_column = 'column_'

# User must define path to save tiles below
save_path = os.path.join('G:\\', 'My Drive', 'data', 'stbd', 'all_line', '4x4')
for a in range(0, mosaic.shape[0]-tile_size, int(tile_size/5)): #iterates through all the rows, with a 5th of the data length as overlapping
    for b in range(0, mosaic.shape[1]-tile_size, int(tile_size/5)):#int(tile_size/2)):#int(tile_size/4)): #iterates through all columns, with a 5th of the data length as overlapping
        mosaic2 = mosaic[a:a+tile_size, b:b+tile_size]
        zero_pixels = np.nonzero(mosaic2[:,:,2] == 0)
        if len(zero_pixels[0]) < (tile_size*tile_size * 0.025):
            file_name = '%s%i_%s%i.png' % (tile_row, a, tile_column, b)
            outfile = os.path.join(save_path, file_name)
            cv2.imwrite(outfile, mosaic2)
