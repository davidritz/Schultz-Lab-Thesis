% cd into directory and order files by alphabetical order
parFile = 'C:\Users\f0046\OneDrive\Desktop\MotherMachineAnalysis\masks\KOmexY14-Jul-2023 08-33-13';
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
    currentImg = im2bw(imread(fname,'Info', info));
    mask(:,:,f) = currentImg;
end
% se = strel('disk',1);
% ermask = imerode(mask,se);
ermask = mask;
perim = bwperim(ermask,4);
se=strel('disk',1,0);
dilMask3=imdilate(perim,se);
se=strel('disk',6,0);
dilMask3=imclose(dilMask3,se);


