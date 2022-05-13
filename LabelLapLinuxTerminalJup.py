# Author: Anne Alsup, anne.alsup@uta.edu
# Large contributer to script Gayathri Nadar, SCF MPI-CBG, nadar@mpi-cbg.de
# https://github.com/jandoG/fiji_jython_useful_functions/blob/e2285c4f3a8e5da160a3d0735408ca902df97da3/code/track_droplets_trackmate_stardist_v0.0.py

import os,math,sys
from datetime import datetime as dt
import csv


from ij import IJ, ImagePlus
from ij import WindowManager
from ij.io import FileSaver 
from ij.measure import Measurements, ResultsTable



from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.util import LogRecorder;
from fiji.plugin.trackmate.tracking.sparselap import SparseLAPTrackerFactory
from fiji.plugin.trackmate.tracking import LAPUtils
from fiji.plugin.trackmate.util import TMUtils
from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate.detection import LabeImageDetectorFactory 
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettings
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.action import CaptureOverlayAction
import java.util.ArrayList as ArrayList

# We have to do the following to avoid errors with UTF8 chars generated in 
# TrackMate that will mess with our Fiji Jython.
reload(sys)
sys.setdefaultencoding('utf-8')

  
# ----------------------------------------------------------------
# 	EDIT PATHS BELOW - MACRO MERGES MASK AND ORIGINAL IMAGES
# ----------------------------------------------------------------

pathpath = os.path.dirname(os.path.realpath(__file__))


# ----------------------------------------------------------------
# 	BEGIN TRACKMATE PROCESS
# ----------------------------------------------------------------
print('Beginning TrackMate process...\n')

# Open compostie image
imp = IJ.openImage(pathpath+'/MaskMerged.tif') 

#imp.show()
print('Compisite image is loaded')


#----------------------------
# Set TrackMate Parameters
#----------------------------

"""
Creates a TrackMate instance configured to operated on the specified
ImagePlus imp. 

params: imp
returns: trackmate object with all the settings 
"""

model = Model()
model.setLogger(Logger.IJ_LOGGER)

settings = Settings(imp)
setup = settings.toStringImageInfo()
        
# Configure Label Detector and LAP Tracker
settings.detectorFactory = LabeImageDetectorFactory() ##Anne 2.9 to 53
settings.detectorSettings = {
    'SIMPLIFY_CONTOURS' : False,
    'TARGET_CHANNEL' : 2
}
print('Detector Settings Configured')

# LAP Tracker Settings 
settings.trackerFactory = SparseLAPTrackerFactory()
settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap() 
settings.trackerSettings['LINKING_MAX_DISTANCE'] = 10.0 #Bigger for larger distance in cell movement (number set for 900x900 image)
settings.trackerSettings['ALLOW_GAP_CLOSING'] = False
settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = False
settings.trackerSettings['ALLOW_TRACK_MERGING'] = False

print('Tracker settings configured')

# analyzers 
settings.addAllAnalyzers()


# trackmate   
trackmate = TrackMate(model, settings)
trackmate.computeSpotFeatures(True)
trackmate.computeTrackFeatures(True)

print("Initiating TrackMate")
#-------------------------
# Initate TrackMate
#------------------------- 

ok = trackmate.checkInput()
if not ok:
    print( str( trackmate.getErrorMessage() ) )
    
ok = trackmate.process()
if not ok:
    print( str( trackmate.getErrorMessage() ) )
    


# Get track model
model = trackmate.getModel()
tm = model.getTrackModel()
fm = model.getFeatureModel()
trackIDs = tm.trackIDs(True)
table = ResultsTable() 

model = trackmate.getModel()
selectionModel = SelectionModel(model)

model.getLogger().log(str(model))
logger = model.getLogger()
settings = trackmate.getSettings()

saveFile = TMUtils.proposeTrackMateSaveFile(settings, logger)
writer = TmXmlWriter(saveFile,logger)
writer.appendLog(logger.toString())
writer.appendModel(trackmate.getModel())
writer.appendSettings(trackmate.getSettings())
writer.writeToFile();

print(str(settings))
sm = SelectionModel( model )

print("Saving XML File")
print("XML file saved to: " + pathpath)

#----------------
# CREATE CSV
#----------------

print('Writing CSV File...')
# The feature model, that stores edge and track features.
fm = model.getFeatureModel()
outpath = (pathpath+'/experimentjup.csv')
with open(outpath,'w') as resultFile:
    csvWriter = csv.writer(resultFile, delimiter=',', lineterminator='\n')
    csvWriter.writerow(['TRACK_ID','POSITION_X','POSITION_Y','FRAME','MEAN_INTENSITY_CH1','MEDIAN_INTENSITY_CH1','MAX_INTENSITY_CH1','MIN_INTENSITY_CH1','AREA'])
# Iterate over all the tracks that are visible.
    for id in model.getTrackModel().trackIDs(True):
        # Get all the spots of the current track.
        track = model.getTrackModel().trackSpots(id)
        for spot in track:
            sid = spot.ID()
            # Fetch spot features directly from spot.
            # Note that for spots the feature values are not stored in the FeatureModel
            # object, but in the Spot object directly. This is an exception; for tracks
            # and edges, you have to query the feature model.
            x=spot.getFeature('POSITION_X')
            y=spot.getFeature('POSITION_Y')
            t=spot.getFeature('FRAME')
            q=spot.getFeature('QUALITY')
            snr=spot.getFeature('SNR_CH1')
            mean=spot.getFeature('MEAN_INTENSITY_CH1')
            median = spot.getFeature('MEDIAN_INTENSITY_CH1')
            maxin = spot.getFeature('MAX_INTENSITY_CH1')
            minin = spot.getFeature('MIN_INTENSITY_CH1')
            area1 = spot.getFeature('AREA')
            model.getLogger().log('\tspot ID = ' + str(sid) + ': x='+str(x)+', y='+str(y)+', t='+str(t)+', q='+str(q) + ', snr='+str(snr) + ', mean = ' + str(mean))
            csvWriter.writerow([str(id+1),str(x),str(y),str(t),str(mean),str(median),str(maxin),str(minin),str(area1)])

print("Finished :)")

exit()
