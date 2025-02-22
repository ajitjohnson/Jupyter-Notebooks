#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat May 28 17:31:18 2022
# @author: Ajit Johnson Nirmal
# Animation in matplotlib

# libs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import matplotlib.colors as colors
import seaborn as sns
import matplotlib.patches as mpatches


# function
def animate (adata, color=None,
             palette=None,
             embedding='umap', 
             x_coordinate='X_centroid', 
             y_coordinate='Y_centroid',
             imageid='imageid', subset=None,
             use_layer=None, use_raw=False, log=False,
             subsample=None,random_state=0,
             n_frames=50, interval=50,forward=True,
             s=None, alpha=1,  cmap='vlag',
             tight_layout=True,plot_legend=False,
             figsize=(5,5),
             save_animation=None,**kwargs):
    """
Parameters:
    ----------
    adata : AnnData Object  
        
    color : list, optional
        Keys for annotations of observations in `adata.obs.columns` or genes in `adata.var.index`. 
        e.g. `color = ['CD3D']` or `color = ['phenotype']`. Please note only one value can be passed at a time.
        The default is None.
        
    palette : dict, optional  
        Colors to use for plotting categorical annotation groups. 
        It accepts a `dict` mapping categories to colors. 
        e.g. `palette = {'T cells': '#000000', 'B cells': '#FFF675'}`.
        Auto color will be generated for categories that are not specified. The default is None.

    embedding : string, optional  
        The `label key` used when running `sm.tl.umap()`. The default is 'umap'.
        
    x_coordinate : string, optional  
        Column that contains the `x_coordinates`. The default is 'X_centroid'.
        
    y_coordinate : string, optional  
        Column that contains the `y_coordinates`. The default is 'Y_centroid'.
        
    imageid : string, optional  
        Name of the column that contains the unique imageid. The default is 'imageid'.
        
    subset : list, optional  
        Unique imageid of a image to be subsetted for plotting. Please note as the coordinate
        system for each images would be unique, only a single image should be passed at a time. 
        Please use this parameter in conjuction with `imageid` to subset a single 
        image. The Function automatically subsets the `UMAP` coordinates. The default is None.
        
    use_layer : string, optional  
        Pass name of any `Layer` in AnnData. The default is `None` and `adata.X` is used.
        
    use_raw : bool, optional  
        If set to `True`, values in `adata.raw.X` will be used to color the plot. The default is False.
        
    log : bool, optional  
        If set to `True`, the data will natural log transformed using `np.log1p()` for coloring. The default is False.
        
    subsample : float, optional  
        Accepts a value between 0-1; Randomly subsamples the data if needed for large images. The default is None.
        
    random_state : int, optional  
        Seed for random number generator. The default is 0.
        
    n_frames : int, optional  
        Number of frames inbetween the UMAP coordinates and the physical coordinates. 
        Higher numbers create a smoother animation. The default is 50.
        
    interval : int, optional  
        interval between frames in milliseconds. The default is 50.
        
    forward : bool, optional  
        If `True` animation will be from `UMAP -> Physical`, if `False` animation 
        will be from `Physical -> UMAP`. The default is True.
        
    s : int, optional  
        The marker size in points. The default is None.
        
    alpha : float, optional  
        blending value, between 0 (transparent) and 1 (opaque). The default is 1.
        
    cmap : string, optional  
        Color map to use for continous variables. Can be a name or a Colormap 
        instance (e.g. "magma”, "viridis"). The default is 'vlag'.
        
    tight_layout : bool, optional  
        Adjust the padding between and around subplots. If True it will ensure that
        the legends are visible. The default is True.
        
    plot_legend : bool, optional  
        Plots the legend. The default is False.
        
    figsize : tuple, optional  
        Width, height in inches. The default is (10, 10).
        
    save_animation : string, optional  
        Pass path to saving animation. Please note depending on the computer specs the live 
        view may not be optimal and hence saving the animation is recommended. 
        e.g `\path\to\directory\figure` The default is None.
        
    **kwargs : Other `matplotlib` parameters.   

Returns:

    Animation
        Can be saved as `gif` using save_animation parameter.

Example:
```python

# Run UMAP
adata = sm.tl.umap(adata)

# Run animation and color it by the identified cell-types
sm.hl.animation (adata, color='phenotype')

```
    """
    
    
    # intrapolation function between co-ordinate sytems
    def tween(e1, e2, n_frames):
        for i in range(5):
            yield e1
        for i in range(n_frames):
            alpha = i / float(n_frames - 1)
            yield (1 - alpha) * e1 + alpha * e2
        for i in range(5):
            yield(e2)
        return
    
    # check if umap tool has been run
    try:
        adata.obsm[embedding]
    except KeyError:
        raise KeyError("Please run `sm.tl.umap(adata)` first")
    
    # identify the coordinates
    umap_coordinates = pd.DataFrame(adata.obsm[embedding],index=adata.obs.index, columns=['umap-1','umap-2'])
    real_coordinates = adata.obs[[x_coordinate,y_coordinate]]
    
    # other data that the user requests
    if color is not None:
        if isinstance(color, str):
            color = [color]
            
        # identify if all elemets of color are available        
        if len(color) > 1:
            raise ValueError("Only a single value in `color` is supported")
            
        # identify if all elemets of color are available        
        if set(color).issubset(list(adata.var.index) + list(adata.obs.columns)) is False:
            raise ValueError("Element passed to `color` is not found in adata, please check!")
        
        # organise the data
        if any(item in color for item in list(adata.obs.columns)):
            adataobs = adata.obs.loc[:, adata.obs.columns.isin(color)]
        else:
            adataobs = None
            
        if any(item in color for item in list(adata.var.index)):
            # find the index of the marker
            marker_index = np.where(np.isin(list(adata.var.index), color))[0]
            if use_layer is not None:
                adatavar = adata.layers[use_layer][:, np.r_[marker_index]]
            elif use_raw is True:
                adatavar = adata.raw.X[:, np.r_[marker_index]]
            else:
                adatavar = adata.X[:, np.r_[marker_index]]
            adatavar = pd.DataFrame(adatavar, index=adata.obs.index, columns = list(adata.var.index[marker_index]))
        else:
            adatavar = None

        # combine all color data
        if adataobs is not None and adatavar is not None:
            color_data = pd.concat ([adataobs, adatavar], axis=1)
        elif adataobs is not None and adatavar is None:
            color_data = adataobs
        elif adataobs is None and adatavar is not None:
            color_data = adatavar        
    else:
        color_data = None
    
    # combine color data with umap coordinates
    if color_data is not None:
        final_data = pd.concat([umap_coordinates, real_coordinates, color_data], axis=1)
    else:
        final_data = umap_coordinates
    
    # subset the final data if nedded
    if subset is not None:
        if isinstance(subset, str):
            subset = [subset]
        cell_to_keep = adata[adata.obs[imageid].isin(subset)].obs.index
        final_data = final_data.loc[cell_to_keep]
    
    # subsample the data if user requests
    if subsample is not None:
        final_data = final_data.sample(frac=subsample, replace=False, random_state=random_state)
    
    # extract the spaces
    e1 = final_data[['umap-1', 'umap-2']].values
    e2 = final_data[[x_coordinate,y_coordinate]].values


    # rescale to same co-ordinates system
    e1[:, 0] -= (max(e1[:, 0]) + min(e1[:, 0])) / 2
    e1[:, 1] -= (max(e1[:, 1]) + min(e1[:, 1])) / 2
    # scale
    scale = max(max(e1[:, 0]) - min(e1[:, 0]), max(e1[:, 1]) - min(e1[:, 1]))
    e1[:, 0] /= scale
    e1[:, 1] /= scale
    # Translate
    e1[:, 0] += 0.5
    e1[:, 1] += 0.5
    
    # rescale co-ordinates
    e2[:, 0] -= (max(e2[:, 0]) + min(e2[:, 0])) / 2
    e2[:, 1] -= (max(e2[:, 1]) + min(e2[:, 1])) / 2
    # scale
    scale = max(max(e2[:, 0]) - min(e2[:, 0]), max(e2[:, 1]) - min(e2[:, 1]))
    e2[:, 0] /= scale
    e2[:, 1] /= scale
    # Translate
    e2[:, 0] += 0.5
    e2[:, 1] += 0.5
    
    # run the interpolation
    if forward is True:
        interpolation = list(tween(e1, e2, n_frames=n_frames))
    else:
        interpolation = list(tween(e2, e1, n_frames=n_frames))
    
    # generate colors
    if s is None:
        s = 100000 / final_data.shape[0]
    
    # if there are categorical data then assign colors to them
    if final_data.select_dtypes(exclude=["number","bool_","object_"]).shape[1] > 0:
        # find all categories in the dataframe
        cat_data = final_data.select_dtypes(exclude=["number","bool_","object_"])
        # find all categories
        all_cat = []
        for i in cat_data.columns:
            all_cat.append(list(cat_data[i].cat.categories))
        
        # generate colormapping for all categories
        less_9 = [colors.rgb2hex(x) for x in sns.color_palette('Set1')]
        nineto20 = [colors.rgb2hex(x) for x in sns.color_palette('tab20')]
        greater20 = [colors.rgb2hex(x) for x in sns.color_palette('gist_ncar', max([len(i) for i in all_cat]))]
        
        all_cat_colormap = dict()
        for i in range(len(all_cat)):
            if len(all_cat[i]) <= 9:
                dict1 = dict(zip(all_cat[i] , less_9[ : len(all_cat[i]) ]   ))
            elif len(all_cat[i]) > 9 and len(all_cat[i]) <= 20:
                dict1 = dict(zip(all_cat[i] , nineto20[ : len(all_cat[i]) ]   ))
            else:
                dict1 = dict(zip(all_cat[i] , greater20[ : len(all_cat[i]) ]   ))
            all_cat_colormap.update(dict1)
        
        # if user has passed in custom colours update the colors
        if palette is not None:
            all_cat_colormap.update(palette)
    else:
        all_cat_colormap = None
        
    # number of plots
    nplots = len(final_data.columns) - 4 # total number of plots
    if nplots > 0:
        column_to_plot = [e for e in list(final_data.columns) if e not in ('umap-1', 'umap-2',x_coordinate,y_coordinate)][0]
        if all_cat_colormap is not None:
            custom_color = list(final_data[column_to_plot].map(all_cat_colormap).values)


    # plot
    fig, ax = plt.subplots(figsize=figsize)
    plt.rcdefaults()
    ax.set(xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
    
    
    if nplots == 0:
        scat = ax.scatter(x = interpolation[0][:, 0], y = interpolation[0][:, 1], s=s, cmap=cmap, alpha=alpha, **kwargs)
        plt.tick_params(right= False,top= False,left= False, bottom= False)
        ax.get_xaxis().set_ticks([]); ax.get_yaxis().set_ticks([])
        if tight_layout is True:
            plt.tight_layout()
    
    if nplots > 0:
        if all_cat_colormap is None:
            scat = ax.scatter(x = interpolation[0][:, 0], y = interpolation[0][:, 1], s=s, 
                           c=final_data[column_to_plot],
                           cmap=cmap, alpha=alpha, **kwargs)
            if plot_legend is True:
                plt.colorbar(scat, ax=ax)
        else:
            scat = ax.scatter(x = interpolation[0][:, 0], y = interpolation[0][:, 1], s=s, 
                           c=custom_color,
                           cmap=cmap, alpha=alpha, **kwargs)
            # create legend
            if plot_legend is True:
                patchList = []
                for key in list(final_data[column_to_plot].unique()):
                    data_key = mpatches.Patch(color=all_cat_colormap[key], label=key)
                    patchList.append(data_key)    
                    ax.legend(handles=patchList,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
        #plt.title(column_to_plot)
        plt.tick_params(right= False,top= False,left= False, bottom= False)
        ax.set(xticklabels = ([])); ax.set(yticklabels = ([]))
        if tight_layout is True:
            plt.tight_layout()
    
 
    def animate(i):
        scat.set_offsets(interpolation[i])
        
    anim = FuncAnimation(
        fig, animate, interval=interval, frames=len(interpolation)-1)
     
    if save_animation is not None:
        anim.save( str(save_animation) + '_scimap.gif', writer='imagemagick', fps=24)

    # save animation
    #anim.save('/Users/aj/Downloads/filename.mp4')
    
