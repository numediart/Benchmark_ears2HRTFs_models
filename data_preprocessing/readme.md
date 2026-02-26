# Preprocessing of data #

## Setup ##

You will need for the 1D processing to download the ear landmark estimation model from: [link](https://github.com/kbulutozler/ear-landmark-detection-with-CNN/blob/master/my_model.h5)
and then to install tensorflow to run it.

## Data processing ##
### 1D and 2D data ###
(before running those, you need to verify that your input data are correctly located ).
You need to run first the 2D preprocessing codes (in order to obtain cropped ear images that will be used to estimates the 1D datas)

You can then run the preprocessing code concerning 1D data.

### 3D data ###
You need for those to run first 3D-mesh-preprocessing.py (this code separates the ears from the head and save the ears)

You can then run the 3D-mesh-to-voxel-auto.py (which then convert the stl mesh into voxels saved as numpy arrays).

## Visualisation ##
The data-stats-plot.py plot the distribution of ears distances accross different modalities to verify the representativity of our data

this can be seen in the following figure:
![fig stats 1D](https://github.com/numediart/Benchmark_ears2HRTFs_models/blob/main/data_preprocessing/Stats%201D.png)

further visualisation for 2D and 3D datas will be added later.
