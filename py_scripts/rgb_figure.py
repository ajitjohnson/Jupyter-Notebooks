import os

import numpy as np
import matplotlib.pyplot as plt

from skimage.external import tifffile
from matplotlib.colors import ListedColormap

from PIL import ImageOps

if __name__ == '__main__':
    # path
    image_filename = 'BaselTMA_SP41_23.475kx17.66ky_10000x5000_12_20170905_89_254_X11Y6_78_a0_full.tiff'
    image_folderpath = '/mnt/d/SCPL_breast_cancer/decompressed_data/OMEnMasks/ome/'

    # channel
    # format: {channel_name: [channel_index, color_name]}
    # note: index here starts from one, not zero
    channel_dict = {
            'CD163': [9, 'cyan'],
            'DNA': [0, 'blue'],
            'CD3D': [10, 'green'],
            'CD45': [15, 'yellow'],
            'KI67': [21, 'red'],
            'CD20': [27, 'magenta'],
            }

    # color code
    # format: {color_name: RGB code}
    # note: range here is [0, 1]
    color_dict = {
            'cyan': [0., 1., 1.],
            'blue': [0., 0., 1.],
            'green': [0., 1., 0.],
            'yellow': [1., 1., 0.],
            'red': [1., 0., 0.],
            'magenta': [1., 0., 1.],
            }

    # load data
    image_filepath = os.path.join(image_folderpath, image_filename)

    with tifffile.TiffFile(image_path) as infile:
        image = infile.series[0].asarray()

    # prepare rgb image
    # note: image shape format is [n_channel, x, y]
    image_rgb = np.zeros((image.shape[1], image.shape[2], 3))
    


    # plot
    fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(12,8))
    axes = axes.flatten()
    for channel_name, ax in zip(channel_dict, axes):
        # unpack input
        channel_index, color_name = channel_dict[channel_name]
        rgb_code = color_dict[color_name]
        layer = image[channel_index-1, ...]
        # normalize layer
        if np.issubdtype(image.dtype, np.integer):
            dtype_max = np.iinfo(image.dtype).max
        elif np.issubdtype(image.dtype, np.floating):
            dtype_max = np.finfo(image.dtype).max
        else:
            raise TypeError('Image dtype not recognized: {}'.format(image.dtype))
        layer_adj = layer.astype(float) / dtype_max
        #layer_adj = auto_contrast.auto_contrast(layer_adj)

        # colormap
        N = 256
        vals = np.ones((N, 4))
        for i, fraction in enumerate(rgb_code):
            vals[:, i] = np.linspace(0, fraction, N)
        custom_cmap = ListedColormap(vals)

        # imshow
        ax.imshow(layer_adj, cmap=custom_cmap)
        ax.set_title(channel_name)
        ax.set_xticks([])
        ax.set_yticks([])

        # add to rgb image
        for i, fraction in enumerate(rgb_code):
            print (fraction)
            image_rgb[..., i] = np.maximum(image_rgb[..., i], layer_adj * fraction)

    fig.tight_layout()
    #plt.savefig('./fig2a_separated.png', dpi=600)
    plt.show()
    plt.close()

    fig = plt.figure()
    plt.imshow(image_rgb)
    plt.xticks([])
    plt.yticks([])
    plt.title('combined')
    #plt.savefig('./fig2a_combined.png', dpi=600)
    plt.show()
    plt.close()
