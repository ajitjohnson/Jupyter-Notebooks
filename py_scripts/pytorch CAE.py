#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 21:52:22 2023

@author: aj
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torchvision.utils as vutils


import tifffile as tiff
import numpy as np
import matplotlib.pyplot as plt
import os
import umap
import seaborn as sns


# load images
img_tensor = load_images_to_tensor('/Users/aj/Desktop/cae')
model = ConvAutoencoder(in_channels=1)
trained_model = train_cae(model, 
                          img_tensor, 
                          epochs=10, 
                          learning_rate=0.0001, 
                          batch_size=20, 
                          verbose=True)

# view the reconstruction
img = tiff.imread('/Users/aj/Desktop/cae/8683.tif')
img = tiff.imread('/Users/aj/Desktop/cae/4108_23.tif')
# plot
img = torch.from_numpy(img).unsqueeze(0).unsqueeze(0).float()
reconstructed_img = trained_model(img)
# viz
reconstructed_img_np = reconstructed_img.detach().numpy().squeeze()
plt.imshow(reconstructed_img_np, cmap='gray')
plt.imshow(img, cmap='gray')

# Define the encoder model
encoder = nn.Sequential(*list(trained_model.children())[:2])

# Pass your data through the encoder
latent_representation = encoder(img_tensor)

# umap on latesnt spcae
embedding = perform_umap(latent_representation)
# plot it
plot_umap(embedding)

# heatmap of latent space
visualize_latent_representation(latent_representation)

# heatmap on the image to see which parts have been used for learning
img = tiff.imread('/Users/aj/Desktop/cae/8683.tif')
encoder = nn.Sequential(*list(trained_model.children())[:2])
visualize_activations(image=img, encoder=encoder)





# CAE

def ConvAutoencoder(in_channels=1, out_channels=10, kernel_size=3, stride=1, padding=1, activation='relu', use_batch_norm=True):
    """
    A basic implementation of a Convolutional Autoencoder (CAE) in PyTorch.

    Parameters:
        in_channels (int): number of input channels
        out_channels (int): number of output channels
        kernel_size (int or tuple): size of the convolutional kernel
        stride (int or tuple): stride of the convolutional kernel
        padding (int or tuple): zero-padding added to both sides of the input
        activation (str): activation function to use ('relu', 'sigmoid', etc.)
        use_batch_norm (bool): flag for using batch normalization after the convolution

    Returns:
        nn.Module: a PyTorch module representing the CAE
    """
    model = nn.Sequential()
    # Encoder
    model.add_module('conv1', nn.Conv2d(in_channels, out_channels, kernel_size, stride, padding))
    if use_batch_norm:
        model.add_module('bn1', nn.BatchNorm2d(out_channels))
    if activation == 'relu':
        model.add_module('relu1', nn.ReLU())
    elif activation == 'sigmoid':
        model.add_module('sigmoid1', nn.Sigmoid())
    # Decoder
    model.add_module('t_conv1', nn.ConvTranspose2d(out_channels, in_channels, kernel_size, stride, padding))
    if use_batch_norm:
        model.add_module('t_bn1', nn.BatchNorm2d(in_channels))
    if activation == 'relu':
        model.add_module('t_relu1', nn.ReLU())
    elif activation == 'sigmoid':
        model.add_module('t_sigmoid1', nn.Sigmoid())
    return model



# Function to process a directory full of images

def load_images_to_tensor(dir_path):
    # Create a list to store all the image tensors
    img_tensors = []
    
    # files
    filenames = [filename for filename in os.listdir(dir_path) if filename.endswith('.tif')]
    
    # Loop through all the files in the directory
    for filename in filenames:
        # Open the image using PIL
        img = tiff.imread(os.path.join(dir_path, filename))
        
        # Convert the numpy array to a tensor
        img_tensor = torch.from_numpy(img).unsqueeze(0).unsqueeze(0).float()
        
        # Append the tensor to the list
        img_tensors.append(img_tensor)
        
    # Concatenate all the image tensors along the first dimension (batch dimension)
    img_tensor = torch.cat(img_tensors, dim=0)
    
    return img_tensor


# Train the model
def train_cae(model, img_tensor, epochs=10, learning_rate=0.001, batch_size=20, verbose=True):
    """
    Trains a Convolutional Autoencoder (CAE) on an image tensor.

    Parameters:
        model (nn.Module): a PyTorch CAE model
        img_tensor (torch.Tensor): image tensor to train the CAE on
        epochs (int): number of training epochs
        learning_rate (float): learning rate for the optimizer
        batch_size (int): batch size to use during training
        verbose (bool): flag for printing the loss during training

    Returns:
        nn.Module: a trained PyTorch CAE model
    """
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

    for epoch in range(epochs):
        for i in range(0, img_tensor.shape[0], batch_size):
            optimizer.zero_grad()
            output = model(img_tensor[i:i + batch_size])
            loss = criterion(output, img_tensor[i:i + batch_size])
            loss.backward()
            optimizer.step()

        if verbose:
            print(f'Epoch {epoch + 1}/{epochs} - Loss: {loss.item()}')

    return model


# UMAP
def perform_umap(latent_representation):
    
    latent_representation = latent_representation.detach().numpy().squeeze()
    
    # Flatten the tensor to 2D
    flattened = latent_representation.reshape(latent_representation.shape[0], -1)
    
    # Perform UMAP
    embedding = umap.UMAP().fit_transform(flattened)
    
    return embedding


# plot UMPA
def plot_umap(embedding):
    # Plot the UMAP embedding
    plt.scatter(embedding[:, 0], embedding[:, 1])
    plt.show()


# heatmap of latent spacce
def visualize_latent_representation(latent_representation):
    # convert to numpy
    latent_representation = latent_representation.detach().numpy().squeeze()
    # Flatten the tensor to 2D
    flattened = latent_representation.reshape(latent_representation.shape[0], -1)
    
    # Visualize the latent representation using a heat map
    sns.heatmap(flattened, cmap='hot')
    plt.show()


# heatmap on the image to see which regions have been used for learning
def visualize_activations(image, encoder):
    # Convert the image to a tensor
    image_tensor = torch.from_numpy(image).unsqueeze(0).unsqueeze(0).float()

    # Register the hook function
    activations = []
    def hook_function(module, input, output):
        activations.append(output)
    
    hooks = []
    # Loop through the layers in the encoder model
    for i, layer in enumerate(encoder):
        if isinstance(layer, nn.Conv2d):
            hook = layer.register_forward_hook(hook_function)
            hooks.append(hook)
    
    # Pass the image through the encoder model
    _ = encoder(image_tensor)

    # Unregister the hook function
    for hook in hooks:
        hook.remove()

    # Plot the activations with contors
    #plt.figure(figsize=(10, 10))
    #plt.axis("off")
    #plt.imshow(activations[0][0][0].detach().numpy(), cmap='gray')
    #plt.show()
    
    # Plot the activations side by side
    #fig, ax = plt.subplots(1, 2, figsize=(10, 10))
    #ax[0].axis('off')
    #ax[0].imshow(image, cmap='gray')
    #ax[1].axis('off')
    #ax[1].imshow(activations[0][0][0].detach().numpy(), cmap='vlag_r')
    #plt.show()
    
    # plot with images overlayed
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.axis('off')
    ax.imshow(image, cmap='gray')
    ax.imshow(activations[0][0][0].detach().numpy(), cmap='vlag_r', alpha=0.5)
    plt.show()







