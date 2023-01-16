%% Load FOV image files

codeDir = 'C:\Users\f0046\OneDrive\Desktop\MotherMachineAnalysis';
cd(codeDir)

% Directory with all FOVs in tif stacks
parFile = '\\dartfs-hpc.dartmouth.edu\rc\lab\S\SchultzLab\Ritz\Pseudomonas data\mmachine\845.7 mmachine\845.7_analyze';

% Directory for flat-field correction tif images
correctFile = '\\dartfs-hpc.dartmouth.edu\rc\lab\S\SchultzLab\Ritz\Pseudomonas data\mmachine\MotherMachineAnalysis\corrections';

% Load files from directories.. loop through all FOVs
fov = 1;

% rawVid(y,x,frame,flu channel)
[rawVid, outName] = loadMother(parFile,fov);

%% Get dark scope images for flat-field correction
cd(codeDir)
if fov == 1
    [AVG_Dark,AVG_TotalGREEN,AVG_TotalBLUE] = loadCorrectMother(correctFile);
end
%% Flat-field correction on the segmentation flu channel, before segmentation

% For vid, only put in 3D with segmentation flu channel selected (CFP = 1, GFP = 2)
% Can also correct other flu channels, but unneeded
[correctSeg,correctFlu] = correctMother(AVG_TotalGREEN,rawVid(:,:,:,2));

%% Segment mother channels (1) automatically, then with (2) user-inputted corrections

% Where mask tif files will be outputted. Each FOV goes into one folder,
% which is named by datetime, OG file name
maskLoc = 'C:\Users\f0046\OneDrive\Desktop\MotherMachineAnalysis\matlab seg\masks\';

% mask will be saved at each iteration into variable and outputted to
% maskLoc..
[mask] = segmentMother(correctSeg,correctFlu,maskLoc,outName);

%% Analyze cells






