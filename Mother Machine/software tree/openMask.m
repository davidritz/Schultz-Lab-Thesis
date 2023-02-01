% cd into directory and order files by alphabetical order
parFile = 'C:\Users\f0046\OneDrive\Desktop\MotherMachineAnalysis\matlab seg\masks\24-Jan-2023 18-25-55';
cd(parFile)
parDir = dir;
[~,ind]=sort({parDir.name});
parDir = parDir(ind);

% Remove any blank files from dir variable
blank = 0;
for f=1:size(parDir,1)
    f = f-blank;
    if parDir(f).bytes == 0
        parDir(f) = [];
        blank=blank+1;
    end
end

% vid(y,x,frame,flu channel)

frame = size(parDir,1);
%%
for f=1:frame
    fname = [parFile,'\',parDir(f).name];
    info = imfinfo(fname);

    % Get tif data for single frame
    if f == 1
        p={info.Width};
        p=p{1};
        q={info.Height};
        q=q{1};
        frames = size(info,1);
        mask = false(q,p,frame);
    end
    currentImg = imread(fname,'Info', info);
    mask(:,:,f) = currentImg;
end



