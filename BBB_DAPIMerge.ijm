close("*");

// ImageJ Macro that will load beta cell and DAPI images
// One DAPI frame will be merged with every Beta cell image 
// Each frame will be saved on its own in the CellPoseImg/ folder

xy = split(getArgument());
print("Directory: " +xy[0]+"");
//print("Pre Beta-cell image: "+ xy[2]+"");
print("Pre DAPI image: "+ xy[2]+"");

//xy = "NoStitch";

if (xy[1] == "NoStitch"){
    print("No Stitching of Control & Experimental Images");
    dir = xy[0];
    PreDAPI = xy[2] + ".tif";
    
    DAPre = dir + PreDAPI;
    PreBeta = dir + "BBB_Pre.png";
    PreSave = dir + "BBB_Cellpose/BBB_"+PreDAPI;
    
    //Merging Pre Ca & DAPI

    open(DAPre);
    open(PreBeta);
    run("Merge Channels...",  "c1=BBB_Pre.png c2="+PreDAPI +" create");
    saveAs("Tiff", PreSave);
    close("*");    
    
    print("DAPI images merged with max images :)");
    print("You can now segment with Cellpose");    
    
    
    
}else{
    dir = xy[0];
    PreDAPI = xy[2] + ".tif";
    PoDAPI = xy[3] + ".tif";
    
    DAPre = dir + PreDAPI;
    DAPost = dir + PoDAPI;

    PreBeta = dir + "BBB_Pre.png";
    PostBeta = dir + "BBB_Post.png";

    PreSave = dir + "BBB_Cellpose/BBB_"+PreDAPI;
    PostSave = dir + "BBB_Cellpose/BBB_"+PoDAPI;

    //Merging Pre Ca & DAPI

    open(DAPre);
    open(PreBeta);
    run("Merge Channels...",  "c1=BBB_Pre.png c2="+PreDAPI +" create");
    saveAs("Tiff", PreSave);
    close("*");

    //Merging Post Ca & DAPI

    open(DAPost);
    open(PostBeta);
    run("Merge Channels...",  "c1=BBB_Post.png c2="+PoDAPI +" create");
    saveAs("Tiff", PostSave);
    close("*");

    print("DAPI images merged with max images :)");
    print("You can now segment with Cellpose");
}
