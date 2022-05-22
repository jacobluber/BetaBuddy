// ImageJ Macro that will load beta cell and DAPI images
// One DAPI frame will be merged with every Beta cell image 
// Each frame will be saved on its own in the CellPoseImg/ folder

xy = split(getArgument());
print("Directory:" +xy[0]+"");
print("Beta-cell image:"+ xy[1]+"");
print("DAPI image:"+xy[2]+"");

dir = xy[0];
FluoImage = xy[1];
DAPI = xy[2];
arg = xy[3];

Beta = dir + FluoImage;
depi = dir + DAPI;
savedapi = dir +"Scaled_DAPI.tif";
scaleog = dir + "Scaled_Beta_Cells.tif";
cellposeimg = dir + "CellPoseImg/DAPIMerge";

open(Beta);
run("Scale...", "x=- y=- z=1.0 width=900 height=900 depth=25 interpolation=Bilinear average process create");
run("Properties...", "channels=1 slices=1 frames=25 pixel_width=0 pixel_height=0 voxel_depth=1.0000000");
saveAs("Tiff",scaleog);
close("*");
print("Beta Image Scaled and Saved :)");

if (arg == "DAPI"){
    open(depi);
    run("Scale...", "x=- y=- width=900 height=900 interpolation=Bilinear average create");
    run("Properties...", "channels=1 slices=1 frames=1 pixel_width=0 pixel_height=0 voxel_depth=1.0000000");
    saveAs("Tiff",savedapi);
    close("*");


    open(scaleog);
    Stack.getDimensions(width, height, channels,slices,frames);
    print(frames);
    run("Stack to Images");
    list1 = getList("image.titles");
    print(list1[9]);

    for (i=0;i<frames;i++){
        open(savedapi);
        dpi = getTitle();
        print("Merging Image: " +list1[i]+"");
        run("Merge Channels...", "c2="+list1[i]+" c3="+dpi+" create");
        run("RGB Color");
        savename = cellposeimg+i;
        saveAs("Tiff", savename);
    }

    close("*");
    print("All images merged :)");
}
