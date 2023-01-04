%% Input mask, GFP tiff exported from python
disp('Input mask tiff: ')
[fname,mfp]= uigetfile('*.tif');
cd (mfp)

disp('Input GFP tiff: ')
[fname2,mfp2]= uigetfile('*.tif');

%% Turn off findpeaks warning, so it doesn't flood the cmd window later
findpeaks([1 1 1 1 ]);
[msg id] = lastwarn;
warning('off',id)
close all

%% Get file info
info = imfinfo(fname);
info2 = imfinfo(fname2);
numberOfImages = length(info);
p={info.Width};
p=p{1};
q={info.Height};
q=q{1};
images = uint16(zeros(q,p,numberOfImages));
imagesRaw = uint16(zeros(q,p,numberOfImages));

%% Create variable with all pics

for k = 1:numberOfImages
    
    % extract plate image for each time point
    currentImg = imread(fname, k, 'Info', info);
    currentImg2 = imread(fname2, k, 'Info', info2);
    
    % put img into variable
    images(:,:,k) = currentImg;
    imagesRaw(:,:,k) = currentImg2;
end

%% Auto tracking
mask=false(size(images));
for k=1:size(images,3)
    % Don't look at media channel
    B = images(:,:,k);
    A = images(1:1000,:,k);
    
    % Find x-location of channels
    xArr=sum(double(A),1);
    [bla,xCh]=findpeaks(xArr,'MinPeakHeight',500);
    
    
    % Search in python mask
    C = mat2gray(B);
    % Scan through mother channels
    for i=1:size(xCh,2)
        x=xCh(i);
        peakCheck = smoothdata(double(sum(B(:,x-30:x+30),2)),'gaussian',[5,5]);
        thresh1 = 50;
        thresh2 = 80;
        [troughs,trlocs] = findpeaks(-peakCheck,'MinPeakProminence',thresh2);
        [pks, pklocs] = findpeaks(peakCheck,'MinPeakDistance',10,'MinPeakHeight',thresh1,'MinPeakProminence',thresh2);
        
        % If only one peak (no trough detected)
        if size(troughs,1) == 0
            trlocs = find(peakCheck(pklocs:end)==0,1,'first')+pklocs;
        end


        %figure(1)
        %findpeaks(peakCheck,'MinPeakDistance',10,'MinPeakHeight',thresh1,'MinPeakProminence',thresh2)
        %figure(2)
        %findpeaks(-peakCheck,'MinPeakProminence',thresh2);
        
        % Identify start and stop of first peak
        y1 = find(-peakCheck,~0,'first');
        y2 = trlocs(1);
        % start at first trough in y, go backwards until hit end of peak 
        % from trough
    
        % scan through x-direction around xCh mid location.. put
        % white mark near troughs at each y
        for y=y1:y2
    
            % Find cell flu peak at each y..
            cellCheck = smoothdata(double(B(y,x-30:x+30)),'gaussian',[6 6]);
            [pks, pklocs] = findpeaks(cellCheck,'MinPeakHeight',20,'MinPeakProminence',20);
            
            if size(pks,2) ~= 0
    
                % Find first and last location from peak where value is less
                % than thresh (start and end tails of peak)

                cellL = find(cellCheck>5,1,'first')-31+xCh(i)+1;
                cellR = (61-find(flip(cellCheck,2)>5,1,'first'))-31+xCh(i)-1;
                mask(y,cellL:cellR,k) = 1;
                mask(y,cellR,k) = 1;

                % overlay on B.. see how it looks

                C(y,cellL) = 1;
                C(y,cellR) = 1;
            end
        end
    end
    %figure(1)
    %imshow(C)

    %figure(2)
    %imshow(mask(:,:,k))
end
%% Mask correction
% assign blobs... click nearest a blob... can subtract all pixels in
% y-direction, up or down.. save a copy of OG, draft1 (all edited ones, including one
% working on), draft2 (all edited, not including current one)


perim = bwperim(mask,4);
se=strel('disk',1,0);
dilMask=imdilate(perim,se);

% show mask starting at frame=1, prompt user if edit needed, then run below
for k=1:size(mask,3)

    % Display raw image frame overlaid with segmentation predictions
    img = mat2gray(imagesRaw(:,:,k));
    ovrImg = dilMask(:,:,k);
    fused = imfuse(img,ovrImg,'falsecolor');
    imshow(fused)

    %imshow(labeloverlay(img,ovrImg,'Transparency',0))
    % Display UI. Can hit stop if happy with segmentation
    zoom on
    title('Press E+Enter if cell on frame needs to be edited. If none, hit only Enter') 
    ButtonHandle = uicontrol('Style', 'PushButton', ...
     'String', 'Stop Button', ...
     'Callback', 'delete(gcbf)');
    
    % Prompt user for edit
    uPrompt = input(['Frame ',num2str(k),':  '],'s');

    if ~ishandle(ButtonHandle)
        % Stop if cancel button was pressed
        close all;
        delete(findall(0));
        break;
    end

    waitforbuttonpress()
    
    try
        [mossPoint(1,1),mossPoint(1,2)] = ginput(1);
        close;
    catch
        close;
        disp('Canceled...')
        disp('  ')
        break
    end

end


% dilate the mask slightly to fill in any gaps
se=strel('disk',2);
dilMask=imdilate(mask,se);

% Find connected regions and their centroids
conn = bwconncomp(dilMask(:,:,100),4);
props = regionprops(conn,'centroid','PixelList','PixelIdxList');

% Mark the blobs
blobs = bwlabel(dilMask(:,:,100),4);



bYX = [bY bX];






















