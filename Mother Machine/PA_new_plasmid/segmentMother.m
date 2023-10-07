function [mask] = segmentMother(correctSeg,correctFlu,maskLoc,outName,rawVid,frstFrame,stp,drugT)
    

    %% Turn off findpeaks warning, so it doesn't flood the cmd window later
    findpeaks([1 1 1 1 ]);
    [msg,id] = lastwarn;
    warning('off',id)
    warning('off','images:bwfilt:tie')
    close all
    disp('Starting auto tracking...')
    %% Auto tracking
    maxW = 10;
    mask=false(size(correctSeg));
    % C=false(size(mask(:,:,1)));
    for k=1:size(correctSeg,3)

        % Don't look at media channel
        B = correctSeg(:,:,k);
        A = correctSeg(1:800,:,k);
        
        % Find x-location of channels
        xArr=sum(double(A),1);
        [bla,xCh]=findpeaks(xArr,'MinPeakDistance',40,'MinPeakHeight',100000,'MinPeakProminence',50000);
        
       
        
        for i = 1:size(xCh,2)
            % Find first bright top pixel.. find x-value of this top pixel..
            % Replace in xCh
    
            % This is in the case of a curved channel
            % Also check in case cell is on edge of image and search area is
            % out of bounds
    
            if xCh(i)-40 < 1
                leftX = 1;
            else
                leftX = xCh(i)-40;
            end
    
            if xCh(i)+40 > size(correctSeg,2)
                rightX = size(correctSeg,2);
            else
                rightX = xCh(i)+40;
            end

            for y=200:800
                findPk = smoothdata(double(B(y,leftX:rightX)),'gaussian',[4,4]);
                [pks, pklocs] = findpeaks(findPk,'MinPeakDistance',60,'MinPeakHeight',5000);
                if size(pks,2) > 0
                    xCh(i) = pklocs - round((rightX-leftX)/2) + xCh(i);
                    break;
                end
            end
        end

        
                    

        % Scan through mother channels IF there are any mother channels
        
        if size(bla,2) > 0
            for i=1:size(xCh,2)

                x=xCh(i);
                
                % Check in case cell is on edge of image and search area is
                % out of bounds
                wid1 = 8;
                if x-wid1 < 1
                    wid1New = x-1;
                
                elseif x+wid1 > size(correctSeg,2)
                    wid1New = size(correctSeg,2) - x;
                else
                    wid1New = wid1;
                end
                wid1 = wid1New;

                peakCheck = smoothdata(double(mean(B(:,x-wid1:x+wid1),2)),'gaussian',[4,4]);
                thresh1 = 5000;
                thresh2 = 3000;
                
                [troughs,trlocs] = findpeaks(-peakCheck,'MinPeakProminence',thresh2,'MinPeakDistance',30);
                [pks, pklocs] = findpeaks(peakCheck,'MinPeakDistance',30,'MinPeakHeight',thresh1,'MinPeakProminence',thresh2);

                % If only one peak (no trough detected)
                if size(troughs,1) == 0 && size(pks,1) == 1
                    y2 = find(flip(peakCheck(pklocs(1):end))>pks(1)/5,1)+pklocs(1);
                elseif size(troughs,1) > 0
                    y2 = trlocs(1) - find(flip(peakCheck(1:trlocs(1)))>-troughs(1)*1.2,1) +1;
                end

                if size(y2,1) == 0 && size(trlocs,1) > 0
                    y2 = trlocs(1);
                end
        
                % Identify start and stop of first peak
                if size(pks,1) > 0
                    y1 = find(-peakCheck<-pks(1)/5,1);
                end

                % If no peak/trough detected, this frame will have same y1,
                % y2 for cell as last frame. This happens if, e.g., cell is too
                % dim for detection.
                
                % scan through x-direction around xCh mid location.. put
                % pixel=1 mark near troughs at each y
                flagPk = 0;
                idx=0;

                for y=y1:y2
                    
                    % Use peak closest to last center..
                    
                    % Find cell flu peak at each y..
                    % 16000 height, 500 prominence

                    wid = 30;

                    % Check in case cell is on edge of image and search area is
                    % out of bounds
                    if xCh(i)-wid < 1
                        leftX = 1;
                        widNew = xCh(i)-1;
                        rightX = xCh(i)+wid;
                    
                    elseif xCh(i)+wid > size(correctSeg,2)
                        widNew = size(correctSeg,2) - xCh(i);
                        rightX = size(correctSeg,2);
                        leftX = xCh(i)-wid;
                    
                    else
                        leftX = xCh(i)-wid;
                        rightX = xCh(i)+wid;
                        widNew = wid;
                    end

                    wid = widNew;
                    

                    cellCheck = smoothdata(double(B(y,leftX:rightX)),'gaussian',[3 3]);
                    [pksCheck, pklocsCheck] = findpeaks(cellCheck,'MinPeakHeight',5000,'MinPeakProminence',500,'MinPeakDistance',5);
                    
                    % Check the amount of cells on this y
                    if size(pksCheck,2) > 1
                        [pks, pklocs] = findpeaks(cellCheck,'MinPeakHeight',16000,'MinPeakProminence',500,'MinPeakDistance',5);
                    else
                        [pks, pklocs] = findpeaks(cellCheck,'MinPeakHeight',5000,'MinPeakProminence',500,'MinPeakDistance',5);
                    end
                    [troughs,trlocs] = findpeaks(-cellCheck,'MinPeakProminence',200);

                    
                    pkThresh = 8000;
                    if size(pks,2) ~= 0

                        % First time in loop, grab first peak
                        if flagPk == 0
                            pkIntense = pks(1);
                            flagPk = 1;
                        else

                            % If two peaks (2 cells) for a certain y-value, use
                            % the peak closest to the last peak on the last
                            % y-value (keep cell choice consistent). If target cell is on the left, x-search
                            % should be from left-most detection until first
                            % flu trough. If target cell is on the right,
                            % x-search should be from flu trough until
                            % right-most detection.
                            if exist('lastPk','var')
                                [~,~,idx]=unique(round(abs(pklocs-lastPk)),'stable');
                            else
                                [~,~,idx]=unique(round(abs(pklocs)),'stable');
                            end
                            pkIntense = pks(idx==1);
                            if size(pkIntense,1) > 1
                                pkIntense = pkIntense(1);
                            end
                        end

                        % If more than 1 cell, keep choice constant and
                        % change segmentation depending on which side
                        % the cell is on.
                        if size(idx,1) > 1
                            xEnd = trlocs(1)+xCh(i)-(wid+1);
                        end
                        
                        
                        % Find first and last location from peak where value is less
                        % than thresh (start and end tails of peak)
                        threshSides = pkIntense/1.5;
                        cellLChk = find(cellCheck>threshSides,1,'first')-(wid+1)+xCh(i);
                        if size(idx,1) == 1 && pkIntense > pkThresh
                            
                            cellL = find(cellCheck>threshSides,1,'first')-(wid+1)+xCh(i);
                            cellR = ((wid*2+1)-find(flip(cellCheck,2)>threshSides,1,'first'))-(wid+1)+xCh(i);
                            
                            mask(y,cellL:cellR,k) = 1;

                        % Find between tracked cell (on left) and daughter cell
                        % (on right). Make sure the x-values detected
                        % aren't too wide
                        elseif size(idx,1) > 1 && cellLChk < xEnd && pkIntense > pkThresh
                            if size(cellLChk:xEnd-5,2) < maxW
                                mask(y,cellLChk:xEnd-5,k) = 1;
                            else
                                mask(y,xEnd-5-maxW:xEnd-5,k) = 1;
                            end

                        % Find between tracked cell (on right) and daughter cell
                        % (on left). Make sure the x-values detected aren't
                        % too wide
                        elseif size(idx,1) > 1 && cellLChk > xEnd && pkIntense > pkThresh 
                            if size(xEnd+5:cellLChk,2) < maxW
                                mask(y,xEnd+5:cellLChk,k) = 1;
                            else
                                mask(y,xEnd+5:xEnd+5+maxW,k) = 1;
                            end
                            
                        end
%                         lastPk = round((cellL+cellR)/2);
%                         maxW = size(cellL:cellR,2)+3;
                        if exist('cellL','var') && size(cellL,2) > 0
                            lastPk = round((cellL+cellR)/2);
                            maxW = size(cellL:cellR,2)+3;
                        end
                        % overlay on img.. see how segmentation looks
        
%                         C(y,cellL) = 1;
%                         C(y,cellR) = 1;
                    end
                end
                % Check for total blobs in mother column, choose largest and discard others

                BW = mask(:,leftX:rightX,k);
                mask(:,leftX:rightX,k) = bwareafilt(BW,1);
            end
%             figure(1)
%             imshow(C)
            
%             figure(1)
%             imshow(correctSeg(:,:,k))
%         
%             figure(2)
%             imshow(mask(:,:,k))
        end 
    end
    %% Mask user correction
    global maskEdit1;
    disp('Finished auto tracking...')
    % If using openMask
%     se=strel('disk',6);
%     dilMask=imclose(dilMask3,se);
%     perim = bwperim(dilMask,4);
%     se=strel('disk',1,0);
%     dilMask=imdilate(perim,se);
%     
%     mask = dilMask;

% Otherwise:
% Dilate mask slightly, so user can see segmentation easier
    se=strel('disk',1);
    mask = imdilate(mask,se);
    se=strel('disk',6);
    dilMask=imclose(mask,se);
    perim = bwperim(dilMask,4);
    se=strel('disk',1,0);
    dilMask=imdilate(perim,se);

    % Update every action
    maskEdit2 = mask;
    dilMask2 = dilMask;
    
    % Update every frame
    maskEdit1 = mask;
    dilMask1 = dilMask;
    %%
    % show mask starting at frame=1, prompt user if edit needed, then run below
    k=1;
    %k=strt;
    flag = 0;
    flag2 = 0;
    flag3 = 0;
    
    while flag == 0
        
        % Reset flags
        flag2 = 0;
        flag3 = 0;
    
        % Circle back to beginning if jump ahead too many frames
        if k > size(mask,3)
            k = k-size(mask,3);
        end
    
        % Display raw image frame overlaid with segmentation predictions
        img = mat2gray(correctFlu(:,:,k));
        ovrImg = dilMask2(:,:,k);
        fused = imfuse(img,ovrImg,'falsecolor');
        figure(1)
        imshow(fused)
        
        % Display UI. Can hit stop if happy with segmentation or type 's'
        
        zoom on
        title('Type S+Enter if happy with segmentation') 
        ButtonHandle = uicontrol('Style', 'PushButton', ...
         'String', 'Stop Button', ...
         'Callback', 'delete(gcbf)');
    
        if ~ishandle(ButtonHandle)
            % Stop if cancel button was pressed
            flag = 1;
            close all;
            %delete(findall(0));
            break;
        end
        
        % Prompt user for edit or go to next frame
    
        if flag == 0
            uPrompt = input(['Frame ',num2str(k),':  '],'s');
            if size(uPrompt,1) == 0
                uPrompt = 'N';
            end
    
            if lower(uPrompt) == 's'
                flag = 1;
                close all
                break;
            end
            if isnan(str2double(uPrompt))

                % If miss-typed, just go to next frame
                if size(uPrompt,2) ~= 1
                    uPrompt = 'N';
                end

                if lower(uPrompt) == 'x' || lower(uPrompt) == 'p' || lower(uPrompt) == 'l' || lower(uPrompt) == 'y'
                    while flag2 == 0
                        
                        if flag3 == 1
    
                            % Update the figure for the user
                            ovrImg = dilMask2(:,:,k);
                            fused = imfuse(img,ovrImg,'falsecolor');
                            imshow(fused)
                            
                            % Ask user what they want to do with that cell
                            uPrompt = input('Following action:  ','s');
                            
                            % Change to 'N' if empty to avoid error
                            if size(uPrompt,1) == 0
                                uPrompt = 'N';
                            end
                        end
    
                        % This if checked the second time and onwards.. don't do
                        % pre-processing if user wants to exit loop
    
                        if uPrompt ~='N' && lower(uPrompt) ~='r'
                            % Click the cell you want to edit
                            % Using ginputWhite instead of ginput to make crosshairs white
                            % instead of black
                            [click(1,1),click(1,2)] = ginputWhite(1);
                            delX = round(click(1));
                            delY = round(click(2));
    
                            % Turn on flag for another user-prompt on this FOV
                            flag3 = 1;
                        
                        % Reset to before edits on this frame
                        elseif lower(uPrompt) == 'r'
                            maskEdit2 = maskEdit1;
                            dilMask2 = dilMask1;
                        end
    
                        % If user puts in X, remove that cell from this frame onwards
                        if lower(uPrompt) == 'x'

                            wid1 = 30;
                            wid2 = 30;

                            % Check in case cell is on edge of image and search area is
                            % out of bounds
                            if delX-wid1 < 1
                                leftX = 1;
                            else
                                leftX = delX-wid1;
                            end
                    
                            if delX+wid2 > size(correctSeg,2)
                                rightX = size(correctSeg,2);
                            else
                                rightX = delX+wid2;
                            end

                            maskEdit2(:,leftX:rightX,k:end) = 0;
                            dilMask2(:,leftX:rightX,k:end) = 0;
                        
                        % If user puts in Y, remove that cell from this
                        % frame until frame 140
                        elseif lower(uPrompt) == 'y'

                            wid1 = 30;
                            wid2 = 30;

                            % Check in case cell is on edge of image and search area is
                            % out of bounds
                            if delX-wid1 < 1
                                leftX = 1;
                            else
                                leftX = delX-wid1;
                            end
                    
                            if delX+wid2 > size(correctSeg,2)
                                rightX = size(correctSeg,2);
                            else
                                rightX = delX+wid2;
                            end
                            % change to 163 to 139
                            maskEdit2(:,leftX:rightX,1:139) = 0;
                            dilMask2(:,leftX:rightX,1:139) = 0;
                        
                        % If user puts in L, remove all pixels lower than
                        % clicked spot

                        elseif lower(uPrompt) == 'l'
                            
                            wid1 = 30;
                            wid2 = 30;

                            % Check in case cell is on edge of image and search area is
                            % out of bounds
                            if delX-wid1 < 1
                                leftX = 1;
                            else
                                leftX = delX-wid1;
                            end
                    
                            if delX+wid2 > size(correctSeg,2)
                                rightX = size(correctSeg,2);
                            else
                                rightX = delX+wid2;
                            end
                            
                            % Remove cell segmentation from figure for user
                            maskEdit2(:,leftX:rightX,k) = 0;
                            dilMask2(:,leftX:rightX,k) = 0;
                
                            fused = imfuse(img,ovrImg,'falsecolor');
                            
                            % Zoom in on cell
                            xlim([delX-20 delX+20]);
                            ylim([delY-30 delY+120]);
                            
                            % Get y-value of user-click
                            [remClick(1,1),remClick(1,2)] = ginputWhite(1);
                            remY = round(remClick(2));
                            
                            % Update mask with everything in channel above click
                            maskEdit2(1:remY+1,leftX:rightX,k) = maskEdit1(1:remY+1,leftX:rightX,k);

                            % Update dilMask2
                            se=strel('disk',6);
                            dilMask2(:,leftX:rightX,k)=imclose(maskEdit2(:,leftX:rightX,k),se);
                            perim = bwperim(dilMask2(:,leftX:rightX,k),4);
                            se=strel('disk',1,0);
                            dilMask2(:,leftX:rightX,k)=imdilate(perim,se);

                        % If user puts in P, erase the original segmentation and replace
                        % with user-added polygon
                        elseif lower(uPrompt) == 'p'
    
                            wid1 = 30;
                            wid2 = 30;

                            % Check in case cell is on edge of image and search area is
                            % out of bounds
                            if delX-wid1 < 1
                                leftX = 1;
                            else
                                leftX = delX-wid1;
                            end
                    
                            if delX+wid2 > size(correctSeg,2)
                                rightX = size(correctSeg,2);
                            else
                                rightX = delX+wid2;
                            end
                            
                            % Remove cell segmentation from figure for user
                            maskEdit2(:,leftX:rightX,k) = 0;
                            dilMask2(:,leftX:rightX,k) = 0;
                
                            fused = imfuse(img,ovrImg,'falsecolor');
                            
                            % Zoom in on cell
                            xlim([delX-30 delX+30]);
                            ylim([delY-50 delY+150]);
                
                            % User draws roi polygon
                            roi = drawpolygon;
    
                            % Select all points within polygon
                            roiX = round(roi.Position(:,1));
                            roiY = round(roi.Position(:,2));
                            
                            % Create roiMask from user ROI
                            roiMask = poly2mask(roiX,roiY,size(mask,1),size(mask,2));
                            roiPerim = bwperim(roiMask,4);
                            se=strel('disk',1,0);
                            dilRoi=imdilate(roiPerim,se);
                            
                            % Update mask with user ROI
                            maskEdit2(:,:,k) = maskEdit2(:,:,k) + roiMask;
                            dilMask2(:,:,k) = dilMask2(:,:,k)+dilRoi;
                        
                        % Else go to next frame, save changes
                        else
                            flag2 = 1;
                            maskEdit1 = maskEdit2;
                            dilMask1 = dilMask2;
                            if lower(uPrompt) ~= 'r'
                                k=k+1;
                            end
                        end
                    end
        
                % Go back one frame
                elseif lower(uPrompt) == 'b'
                    k=k-1;
                    if k == 0
                        k = size(mask,3);
                    end
                
                % display graph of all data over mask
                elseif lower(uPrompt) == 'g'
                    % Analyze mask so far
                    se=strel('disk',6);
                    maskCheck = imclose(dilMask1,se);
                    blobsCheck = analyzeMother(maskCheck,rawVid);
                    % strt, stop, drug
                    graphMother(blobsCheck,frstFrame,stp,drugT)
                    clear maskCheck
                    clear blobsCheck

                % Else exit the while loop if nothing entered and go to next frame
                else
                    flag2 = 1;
                    k=k+1;
                end
    
            % Go to user-inputted frame
            elseif flag ~= 1
                k=str2double(uPrompt);
            end
        end
    end
   
    se=strel('disk',6);
    mask = imclose(dilMask1,se);
    maskEdit1 = mask;
    %% Save as tiff

    % Create a folder from datetime and put into masks directory
    dateString = datestr(datetime);
    dateString = strrep(dateString,':','-');
    fileString = [maskLoc,dateString];
    mkdir(fileString)
    cd(fileString)
    filePre = ['segMother_',outName(1:end-7)];
    
    % Output mask as tifs
    for k = 1:size(maskEdit1,3)
        fileNum = sprintf('%03d',k);
        fileName = [filePre,fileNum,'.tif'];
        imwrite(maskEdit1(:,:,k),fileName)
    end
    mask = maskEdit1;
end
















