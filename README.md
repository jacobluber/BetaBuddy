BetaBuddy
=========

Your very own automated pipeline for Beta-cell calcium imaging segmentation, registration, tracking, and analysis!

Written by Anne Alsup and Kelli Fowlds through The University of Texas at Arlington.

BetaBuddy is run through the Jupyter Notebook BetaBuddy.ipynb. Follow the steps bellow to install all necessary packages and 

How to Use BetaBuddy
====

Our pipeline has been tested on ND2 and TIFF Beta-cell image sequences. Other cell types have not been tested on this pipeline. Your image set will need to be in the same directory as BetaBuddy.ipynb. 

### Cell 1

If you need to convert ND2 images into TIFF, change the filenames after `!./bftools/bfconvert-overwrite ./`. The first filename will call on your current ND2 image. Due to the Linux OS you cannot have any spaces in the name. The second filename is how you want to save your TIFF file. Change the names for your DAPI and calcium dye images.

You will now have two TIFF files in your current working directory.

*Note: You can comment out the first line if you don't have a DAPI image*

### Cell 2

A directory for the merged DAPI and calcium dye images is cre 

Installation
====

### System Requirements
This pipeline must be run on a Linux OS preferably with GPU access and at least 16GB for image stacks. 


### Directory Setup

The directory that hold your BetaBuddy.ipynb must also have the following:
* Calcium image stack and/or DAPI image stack as TIFF or ND2 file(s)
* DAPIMerge.ijm
* MaskMerge.ijm
* LabelLapLinuxTerminalJup.py
* bftools unzipped directory
* Beta_Buddy_Stats.R
* Fiji.app unzipped directory 

## Package Dependencies
To use the widgets within the Jupyter Notebook:
~~~sh
pip3 install ipywidgets
jupyter nbextension enable --py widgetsnbextension
~~~

~~~sh
pip3 install skitimage
~~~

### Cellpose Dependencies 
CUDA and Pytorch will be necessary to run Cellpose. The following command is for CUDA 11.0, but you can find the correct CUDA installation command for you [here](https://pytorch.org/get-started/locally/). 

~~~sh
conda install pytorch torchvision torchaudio cudatoolkit=11.0 -c pytorch
~~~

You can pip or conda install Cellpose:

~~~sh
conda install cellpose
or
pip3 install cellpose
~~~

Spefici OpenCV installation for Cellpose:

~~~sh
pip3 install "opencv-python-headless<4.3"
~~~


## Bio-Format Tools Installation

The Bio-Formats software tool will allow us to convert our ND2 files into TIFF files. The "bftools" directory will need to be inside the folder with your BetaBuddy.ipynb

~~~sh
wget https://downloads.openmicroscopy.org/bio-formats/6.9.1/artifacts/bftools.zip
unzip bftools
~~~

## Fiji Installation

You can install Fiji through the terminal with the following commands. Once you unzip the folder, you will need to move the 'Fiji.app' folder out of 'fiji-linux64.zip' and directly into the directory with BetaBuddy.ipynb.

~~~sh
wget https://downloads.imagej.net/fiji/latest/fiji-linux64.zip
unzip fiji-linux64.zip
~~~





