import numpy as np
import matplotlib.pyplot as plt
from rasterio.transform import Affine
import rasterio
import SSS_functions as sss
import os

def grid_data(X, Y, X_utm, Y_utm, data, no_pings, resampling, resolution, file_number, file,
              progress_callback=None, stage_offset=0, max_gap_factor=1.0):
    """
    Vectorized gridding with ADAPTIVE along-track interpolation.
    """
    mosaic = np.zeros((len(Y), len(X)))
    file_name = os.path.splitext(os.path.basename(file[file_number]))[0]

    # Stage 1: Slice data
    X_sub = X_utm[:no_pings, :, file_number]
    Y_sub = Y_utm[:no_pings, :, file_number]
    data_sub = data[:, :no_pings, file_number].T
    if progress_callback:
        progress_callback(stage_offset + 1, f"[{file_name}] Preparing array slices...")

    # Stage 2: Determine how many interpolated points each ping-to-ping
    # interval needs. We look at every column (bin) across the swath and use
    # the worst-case (largest) consecutive-ping distance per interval
    if no_pings > 1:
        dx = np.diff(X_sub, axis=0)
        dy = np.diff(Y_sub, axis=0)
        gap_dist = np.sqrt(dx ** 2 + dy ** 2).max(axis=1)  # (no_pings - 1,)

        target_spacing = max(resolution * max_gap_factor, 1e-6)
        subdivisions = np.maximum(1, np.ceil(gap_dist / target_spacing).astype(int))

        # Build fractional "virtual ping index" query positions, fully
        # vectorized (no per-ping Python loop): interval i contributes
        # subdivisions[i] points at i, i + 1/n, i + 2/n, ...
        interval_start_idx = np.repeat(np.cumsum(subdivisions) - subdivisions, subdivisions)
        offset_in_interval = np.arange(subdivisions.sum()) - interval_start_idx
        query = np.repeat(np.arange(no_pings - 1), subdivisions) + offset_in_interval / np.repeat(subdivisions, subdivisions)
        query = np.append(query, no_pings - 1)
    else:
        query = np.array([0.0])
    if progress_callback:
        progress_callback(stage_offset + 2, f"[{file_name}] Computing adaptive sample density...")

    # Stage 3: Linear interpolation of UTM X/Y and sonar data at the query
    # positions (same query positions reused for all three so they stay
    # aligned -- equivalent to the previous order=1 ndimage.zoom, just
    # non-uniformly spaced).
    floor_idx = np.minimum(np.floor(query).astype(int), max(no_pings - 2, 0))
    ceil_idx = np.minimum(floor_idx + 1, no_pings - 1)
    frac = (query - floor_idx)[:, None]

    X_dense = X_sub[floor_idx] * (1 - frac) + X_sub[ceil_idx] * frac
    Y_dense = Y_sub[floor_idx] * (1 - frac) + Y_sub[ceil_idx] * frac
    data_dense = data_sub[floor_idx] * (1 - frac) + data_sub[ceil_idx] * frac
    if progress_callback:
        progress_callback(stage_offset + 3, f"[{file_name}] Interpolating UTM coordinates and sonar data...")

    # Stage 4: Calculate grid indices
    id_x = np.clip(np.round((X_dense - X[0]) / resolution).astype(int), 0, len(X) - 1)
    id_y = np.clip((len(Y) - 1) - np.round((Y_dense - Y[0]) / resolution).astype(int), 0, len(Y) - 1)
    if progress_callback:
        progress_callback(stage_offset + 4, f"[{file_name}] Mapping to raster matrix...")

    # Stage 5: Assign and cleanup (in-place, avoids an extra full-size copy)
    mosaic[id_y, id_x] = data_dense
    mosaic = np.nan_to_num(mosaic, copy=False)
    if progress_callback:
        progress_callback(stage_offset + 5, f"[{file_name}] Finalizing raster grid...")

    return mosaic

def scalar_to_rgb(data, color_min, color_max):
    norm = plt.Normalize(vmin=color_min, vmax=color_max)
    cm = sss.Color_palette()

    rgb_data = cm(norm(data))
    rgb_data = (rgb_data[:, :, :3] * 255).astype('uint8')
    return rgb_data

def save_mosaic(mosaic, x_utm, y_utm, resolution, filepath, color_min, color_max):
    img_rgb = scalar_to_rgb(mosaic, color_min, color_max)

    X = np.arange(np.min(x_utm) - 100, np.max(x_utm) + 100, resolution)
    Y = np.arange(np.min(y_utm) - 100, np.max(y_utm) + 100, resolution)
    transform = Affine.translation(np.min(X) - resolution / 2, np.min(Y) - resolution / 2) * Affine.scale(resolution, resolution)

    image = np.moveaxis(img_rgb.squeeze(), -1, 0)
    file_name_with_ext = os.path.basename(filepath[0])
    file_name = os.path.splitext(file_name_with_ext)[0]
    
    new_dataset = rasterio.open(
        f'{file_name}.tif',
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
