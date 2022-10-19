# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 17:20:02 2022

@author: ajn16
"""
import napari
from oiffile import imread
from mrc import DVFile, imread
import os

image = imread('Y:/lsp-data/cycif-techdev/confocal/20220311/Sample20um_FoV3_CD11C_MHC1_SOX10_TD.oib')
image = imread("Y:\lsp-analysis\cycif-techdev\ThickSamplesDeltavision\melanomacycle1\1-520b03_R3D_D3D.dv")
image = imread("Y:/lsp-analysis/cycif-techdev/ThickSamplesDeltavision/melanomacycle1/1-520b03_R3D_D3D.dv")


viewer = napari.view_image(image)
viewer = napari.view_image(image, channel_axis=0, name= ['Something','SOX10','MHC1','CD11C'])



import tifffile


image_path = "//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/PTCL_jackson_temporary/ometif/BatchA/melanoma.ome.tif"
image_path = "//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/PTCL_jackson_temporary/ometif/BatchB/12.ome.tif"
image_path = "//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/PTCL_jackson_temporary/ometif/BatchB/15.ome.tif"
image_path = "//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/PTCL_jackson_temporary/ometif/BatchB/20.ome.tif"
image_path = "//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/PTCL_jackson_temporary/ometif/BatchA/crc.ome.tif"


os.chdir("D:/tif_files/batchA")
os.chdir("D:/tif_files/batchB/12")
os.chdir("D:/tif_files/batchB/15")
os.chdir("D:/tif_files/batchB/20")
os.chdir("D:/tif_files/crc")


tiff = tifffile.TiffFile(image_path)
n_channels = tiff.series[0].shape[0]

for i in range(0, n_channels, 5):
    img = tiff.series[0][i].asarray()
    tifffile.imsave(f"output-channel_{i:02}.tif", img, tile=(1024, 1024), bigtiff=True)

tiff.close()


# running Jeremey scriptt

python figure_registration_flow.py batchA/output-channel_00.tif batchA/output-channel_05.tif output_1_2.tif

python figure_registration_flow.py batchB/12/output-channel_00.tif batchB/12/output-channel_05.tif batchB_1_2.tif

python figure_registration_flow.py crc/output-channel_00.tif crc/output-channel_5.tif crc_1_2.tif --intensity-threshold 2000