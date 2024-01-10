inputs = split(getArgument());
print("Directory: " + inputs[0] +"");
close("*");

//Comment out if you use the script within Fiji
dir = inputs[0];
arg = inputs[1];
Bpre = inputs[2];
Dpre = inputs[3];
Bpost = inputs[4];
Dpost = inputs[5];

//Uncomment if you use the script within Fiji
//dir = "/home/BetaBuddy/";
//arg = "NoStitch";
//Bpre = "1VPre";
//Dpre = "1VDAPIPre";
//Bpost = "1VPost";
//Dpost = "1VDAPIPost";

//Defining directories
maskdir = dir + "BBB_Cellpose/Masks/";

//pre ca2+ dir
predir = dir + Bpre + ".tif";

//pre mask dir
m_predir = maskdir + "MASK_BBB_" + Dpre + ".tif";

//pre merge stack savename
pre_savename = "BBB_" + Bpre + "_Merge.tif"
pre_save = dir + pre_savename;

//total mask merge savename
mergesave = dir + "BBB_MaskMerge.tif";

//pre mask stack savename
interimpre = "BBB_" + Bpre + "_Mask.tif"
premasksave = dir + interimpre;

open(predir);
Stack.getDimensions(width, height, channels,slices,frames);
print(frames);

for (i=0;i<frames;i++){
	open(m_predir);
}
run("Images to Stack", "name=Masks title=[] use");
saveAs("Tiff", premasksave);

run("Merge Channels...", "c1="+Bpre+".tif c4="+interimpre+" create");
run("Properties...", "channels=2 slices=1 pixel_width=0 pixel_height=0 voxel_depth=1.0000000");
saveAs("Tiff", pre_save);


if (arg != "NoStitch"){
    
    print("Stitching of 2 Images Sets");
    
    //post ca2+ dir
    postdir = dir + Bpost + ".tif";

    //post mask dir
    m_postdir = maskdir + "MASK_BBB_" + Dpost + ".tif"; 

    //post merge stack savename
    post_savename = "BBB_" + Bpost + "_Merge.tif";
    post_save = dir + post_savename;
    
    //post mask stack saename
    interimpost = "BBB_" + Bpost + "_Mask.tif";
    pomasksave = dir + interimpost;
    
    
    open(postdir);
    Stack.getDimensions(width, height, channels,slices,frames);
    print(frames);


    for (i=0;i<frames;i++){
        open(m_postdir);
    }
    
    run("Images to Stack", "name=PostMasks title=[] use");
    saveAs("Tiff", pomasksave);

    run("Merge Channels...", "c1="+Bpost+".tif c4="+interimpost+" create");
    saveAs("Tiff",post_save);

    run("Concatenate...", "  title=MaskMerge keep image1="+pre_savename+" image2="+post_savename+"");

    run("Properties...", "channels=2 slices=1 pixel_width=0 pixel_height=0 voxel_depth=1.0000000");


    saveAs("TIFF",mergesave);
  
}else{

    saveAs("TIFF", mergesave);
}



print("Mask Merged TIFF Saved!");
print("Tracking will now start :)");







