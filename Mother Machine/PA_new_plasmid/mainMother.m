%% Eventually uncomment this
%blobsGlobal = 0;
global maskEdit1;
%% Load FOV image files

codeDir = 'C:\Users\f0046\OneDrive\Desktop\MotherMachineAnalysis\matlab seg';

% Directory with all FOVs in tif stacks

%parFile = '\\dartfs-hpc.dartmouth.edu\rc\lab\S\SchultzLab\Ritz\Pseudomonas data\mmachine\845.7 mmachine\845.7_analyze';
parFile = 'D:\all pa14\komexy';
%parFile = '\\dartfs-hpc.dartmouth.edu\rc\lab\S\SchultzLab\Ritz\Pseudomonas data\mmachine\MotherMachineAnalysis\WT_PA14_vid';

% Directory for flat-field correction tif images

%correctFile = '\\dartfs-hpc.dartmouth.edu\rc\lab\S\SchultzLab\Ritz\Pseudomonas data\mmachine\MotherMachineAnalysis\corrections';
correctFile = 'D:\background';

% cd into directory and order files by alphabetical order
cd(parFile)
parDir = dir;
[~,ind]=sort({parDir.name});
parDir = parDir(ind);

% Remove any blank files from dir variable
blank = 0;
for fov=1:size(parDir,1)
    fov = fov-blank;
    if parDir(fov).bytes == 0
        parDir(fov) = [];
        blank=blank+1;
    end
end

fovTot = size(parDir,1)/2;
drugAdd = 243;
strt = 123;
fnsh = strt+200;

%% Load files from directories.. loop through all FOVs
for fov = 5:fovTot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load raw files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cd(codeDir)
    % rawVid(y,x,frame,flu channel); GFP = 1, RFP = 2
    [rawVid, outName] = loadMother(parFile,fov);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get dark scope images for flat-field correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cd(codeDir)
    if fov == 1
        [AVG_TotalGREEN,AVG_TotalRED] = loadCorrectMother(correctFile);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Flat-field correction on the segmentation flu channel, before segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cd(codeDir)
    % For vid, only put in 3D with segmentation flu channel selected (GFP = 1, RFP = 2)
    % Can also correct other flu channels, but unneeded
    
    [correctSeg,correctFlu] = correctMother(AVG_TotalGREEN,rawVid(:,:,:,1));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Segment mother channels (1) automatically, then with (2) user-inputted corrections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cd(codeDir)

    % Where mask tif files will be outputted. Each FOV goes into one folder,
    % which is named by datetime, OG file name
    maskLoc = 'C:\Users\f0046\OneDrive\Desktop\MotherMachineAnalysis\masks\KOmexY';
    
    % mask will be saved at each iteration into variable and outputted to
    % maskLoc..
    correctSeg = correctSeg(:,:,strt:fnsh);
    correctFlu = correctFlu(:,:,strt:fnsh);
    
    mask = segmentMother(correctSeg,correctFlu,maskLoc,outName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyze cells in FOV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    poo
    cd(codeDir)
    % add to remove blank frames
    %mask = mask(:,:,123:323);
    rawVid = rawVid(:,:,strt:fnsh,:);

    if max(max(max(mask(:,:,1)))) ~= 0
        num = sprintf('%03d',fov);
        blobsName = ['blobsGlobal_KOmexY_',num,'.mat'];
        % analyze cells and append cells to blobsGlobal
        if fov == 1
            [blobsGlobal] = analyzeMother(mask,rawVid);
        else
            [blobsAppend] = analyzeMother(mask,rawVid);
            blobsGlobal = [blobsGlobal;blobsAppend];
        end
        save(blobsName,'blobsGlobal')
    end

end
%% Analyze all cells in blobsGlobal













