inputs = split(getArgument());
print("Directory: " + inputs[0] +"");

dir = inputs[0];
maskdir = dir + "CellPoseImg/Masks/";
scaleddir = dir+ "Scaled_Beta_Cells.tif";
savename = dir + "MaskMerged.tif";
masksave = dir + "MaskStack.tif";

File.openSequence(maskdir);
saveAs("Tiff",masksave);
maskname = getTitle();
open(scaleddir);
betaname = getTitle();
run("Merge Channels...", "c2="+betaname+" c4= "+maskname+" create");

saveAs("Tiff", savename);
print("Merged Mask/Beta Cell Image Saved");
