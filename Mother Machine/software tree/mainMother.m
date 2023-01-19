%% Eventually uncomment this
%blobsGlobal = 0;
%% Load FOV image files

codeDir = 'C:\Users\f0046\OneDrive\Desktop\MotherMachineAnalysis';

% Directory with all FOVs in tif stacks
parFile = '\\dartfs-hpc.dartmouth.edu\rc\lab\S\SchultzLab\Ritz\Pseudomonas data\mmachine\845.7 mmachine\845.7_analyze';

% Directory for flat-field correction tif images
correctFile = '\\dartfs-hpc.dartmouth.edu\rc\lab\S\SchultzLab\Ritz\Pseudomonas data\mmachine\MotherMachineAnalysis\corrections';

% 
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


% Load files from directories.. loop through all FOVs
%fov = 1:fovTot;
for fov = 2:fovTot
    %% Load raw files
    cd(codeDir)
    % rawVid(y,x,frame,flu channel); CFP = 1, GFP = 2
    [rawVid, outName] = loadMother(parFile,fov);
    
    %% Get dark scope images for flat-field correction
    cd(codeDir)
    if fov == 1
        [AVG_Dark,AVG_TotalGREEN,AVG_TotalBLUE] = loadCorrectMother(correctFile);
    end
    %% Flat-field correction on the segmentation flu channel, before segmentation
    cd(codeDir)
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
    
    %% Analyze cells in FOV
    cd(codeDir)
    if max(max(max(mask))) ~= 0
        % analyze cells and append cells to blobsGlobal
        if fov == 1
            [blobsGlobal] = analyzeMother(mask,rawVid);
            %save('blobsGlobal.mat','blobsGlobal')
        else
            [blobsAppend] = analyzeMother(mask,rawVid);
            blobsGlobal = [blobsGlobal;blobsAppend];
            % blobsAppend CFP flu values look weird.. need to assign
            % correct flu to correct cells as well.
            poo
            %save(filename,variables)
        end
    end
end
%% Analyze all cells in blobsGlobal





