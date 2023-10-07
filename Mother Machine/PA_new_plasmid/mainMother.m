%% Eventually uncomment this
%blobsGlobal = 0;
global maskEdit1;

%% User inputs
strt = 110;
frstFrame = 1;
stp = 170+strt;
drugAdd = 167-strt;
drugT = drugAdd;
%% Load FOV image files

codeDir = 'C:\Users\f0046\OneDrive\Desktop\MotherMachineAnalysis\matlab seg';

% Directory with all FOVs in tif stacks
parFile = 'D:\all pa14\komexZ';

% Directory for flat-field correction tif images
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

fnsh = strt+300;
%drugAdd = 120-strt;

%% Load files from directories.. loop through all FOVs
for fov = 1:fovTot
%% Preliminary filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load raw files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fov=10;
    cd(codeDir)
    % rawVid(y,x,frame,flu channel); GFP = 1, RFP = 2
    [rawVid, outName] = loadMother(parFile,fov);
    rawVid = rawVid(:,:,strt:fnsh,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get dark scope images for flat-field correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cd(codeDir)
    % fov == 1
    if fov == fov
        [AVG_TotalGREEN,AVG_TotalRED] = loadCorrectMother(correctFile);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Flat-field correction on the segmentation flu channel, before segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cd(codeDir)
    % For vid, only put in 3D with segmentation flu channel selected (GFP = 1, RFP = 2)
    % Can also correct other flu channels, but unneeded
    
    [correctSeg,correctFlu] = correctMother(AVG_TotalGREEN,rawVid(:,:,:,1),drugAdd);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Segment mother channels (1) automatically, then with (2) user-inputted corrections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cd(codeDir)

    % Where mask tif files will be outputted. Each FOV goes into one folder,
    % which is named by datetime, OG file name
    maskLoc = 'C:\Users\f0046\OneDrive\Desktop\MotherMachineAnalysis\masks\wt';
    
%%  Segmentation

    % mask will be saved at each iteration into variable and outputted to
    % maskLoc..
    mask = segmentMother(correctSeg,correctFlu,maskLoc,outName,rawVid,frstFrame,stp,drugAdd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyze cells in FOV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis    
    cd(codeDir)

    if max(max(max(mask(:,:,1)))) ~= 0
        num = sprintf('%03d',fov);
        blobsName = ['blobsGlobal_wt_',num,'.mat'];
        % analyze cells and append cells to blobsGlobal
        if fov == 1
            [blobsGlobal] = analyzeMother(mask,rawVid);
        else
            % 188 to 164
            [blobsAppend] = analyzeMother(mask(:,:,140:end),rawVid(:,:,140:end,:));
            blobsGlobal = [blobsGlobal;blobsAppend];
        end
        save(blobsName,'blobsGlobal')
        % Can toss blobsGlobal #32 KOmexY.. noisy growth
        % Toss #14 KOmexY? Isn't growing.. very bright
        % KOmexZ: start at 110, drug at 177, finish at 310
        % KOmexY: start at 123, drug at 243, finish at 323
        % KOmexZ missing fov6, 2nd cell (hard to capture bc on the edge and
        % curved channel).. inputted as cell #29
        % KOmexZ: not using first cell on FOV9.. too dim for segmentation
        % at later frames
        % wt: 60 strt, 360 fnsh, drug add at 120-60
        % wt cell 3: mother cell lyses at frame 30, new mother cell used
        % for segmentation
        % wt cell 10: grew into mother channel just before drug
        % wt surviving cells: 4,7,11,13,14,16,22,23,24,25,26,28
        % KOmexZ: cell 10 survives, but gets pushed out of FOV by
        % non-glowing mother cell after frame 139
        % KOmexZ: no fov5 in blobsGlobal?
        % KOmexZ: cell 9 dead
        % KOmexZ: cell 19 lyses
        % KOmexZ: cell 30 survives? not sure where this one is
    end

end
%% Analyze all cells in blobsGlobal













