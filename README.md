BetaBuddy
=========

Your very own automated pipeline for Beta-cell calcium imaging segmentation, registration, tracking, and analysis!

Written by Anne Alsup and Kelli Fowlds through The University of Texas at Arlington.

BetaBuddy is run through the Jupyter Notebook BetaBuddy.ipynb. Follow the steps bellow to run the pipeline and install all necessary packages. 

How to Use BetaBuddy & Step Breakdown
====

Our pipeline has been tested on ND2 and TIFF Beta-cell image sequences. Other cell types have not been tested on this pipeline. Your image set will need to be in the same directory as BetaBuddy.ipynb. 

### Cell 1

If you need to convert ND2 images into TIFF, change the filenames after `!./bftools/bfconvert-overwrite ./`. The first filename will call on your current ND2 image. Due to the Linux OS you cannot have any spaces in the name. The second filename is how you want to save your TIFF file. Change the names for your DAPI and calcium dye images.

You will now have two TIFF files in your current working directory.

*Note: You can comment out the first line if you don't have a DAPI image*

### Cell 2

A directory for the merged DAPI and calcium dye images is created. *If you do not have a DAPI channel, move your calcium dye TIF image(s) into the new CellPoseImg directory and comment out the* `./Fiji.app` *command.* The DAPI and calcium dye images are merged with an ImageJ macro and saved in the CellPoseImg directory.

### Cells 3-6

Configuration and parameterization of Cellpose. Cells 4 and 6 need to be run to have the widget options appear in the output section. Following the instructions before these cells is very important to run Cellpose correctly. Cell 5 will show you which channels of your image contains the calcium dye and DAPI.

### Cell 7
Test Cellpose on only one image to check on your parameters. You can always go back and adjust your parameters if the sample output looks incorrect.

### Cell 8

Run Cellpose on all images and save the masks in the Mask directory in CellPoseImg. This step may take around 5-15 minutes depending on the size of your image stack. 

###  Cell 9

The first command will merge your new mask images with your original Beta cell images with an ImageJ macro. The second command will run this merged image stack through the Fiji plugin TrackMate. 

An XML titled *MaskMerged.xml* will be saved in the current directory along with the *MaskMerged.tif*. You can load this into Fiji through *Plugins -> Tracking -> Load TrackMate File* and see the calculated tracks. You can also open the *Track* tab and click on TrackIDs to see their location. 

A CSV will be saved with fluorescent intensities from every cell grouped together by their TrackID over time. This CSV will be saved in the current directory as *experimentjup.csv*.

*Note: The TrackIDs in the CSV files will be off by one*

### Cell 10

Background subtraction is performed by sampling 100 pixels outside of the ROIs created by Cellpose. These pixel intensities are saved as a CSV.

### Cell 11

Last step!! The R script Beta_Buddy_Stats.R is run from the terminal to conduct background subtraction, normalization, preliminary statistical analyses, and data visualization. Three plots and the final CSV file has only tracks that span the full image stack and identifies tracks in the top 20% of signal variance over all frames.


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

### Cellpose Dependencies 
CUDA and Pytorch will be necessary to run Cellpose. The following command is for CUDA 11.0, but you can find the correct CUDA installation command for you [here](https://pytorch.org/get-started/locally/). 

~~~sh
conda install pytorch torchvision torchaudio cudatoolkit=11.0 -c pytorch
~~~

You can pip or conda install Cellpose:

~~~sh
python -m pip install cellpose
~~~

Spefici OpenCV installation for Cellpose:

~~~sh
pip3 install "opencv-python-headless<4.3"
~~~

~~~sh
pip3 install scikit-image
~~~

~~~sh
pip3 install matplotlib
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





