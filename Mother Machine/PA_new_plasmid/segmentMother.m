function [mask] = segmentMother(correctSeg,correctFlu,maskLoc,outName)


    %% Turn off findpeaks warning, so it doesn't flood the cmd window later
    findpeaks([1 1 1 1 ]);
    [msg,id] = lastwarn;
    warning('off',id)
    warning('off','images:bwfilt:tie')
    close all
    
    %% Auto tracking
    maxW = 10;
    mask=false(size(correctSeg));
    % C=false(size(mask(:,:,1)));
    %for k=1:size(correctSeg,3)
    % KOy -> start = 123 ; end = 123+200
    strt = 1;
    fnsh = 201;
    for k=strt:fnsh

        % Don't look at media channel
        B = correctSeg(:,:,k);
        A = correctSeg(1:800,:,k);
        
        % Find x-location of channels
        xArr=sum(double(A),1);
        [bla,xCh]=findpeaks(xArr,'MinPeakDistance',40,'MinPeakHeight',100000,'MinPeakProminence',50000);

        % Scan through mother channels IF there are any mother channels

        if size(bla,2) > 0
            for i=1:size(xCh,2)
                x=xCh(i);
                peakCheck = smoothdata(double(mean(B(:,x-8:x+8),2)),'gaussian',[4,4]);
                thresh1 = 5000;
                thresh2 = 500;
                [troughs,trlocs] = findpeaks(-peakCheck,'MinPeakProminence',thresh2,'MinPeakDistance',30);
                [pks, pklocs] = findpeaks(peakCheck,'MinPeakDistance',30,'MinPeakHeight',thresh1,'MinPeakProminence',thresh2);

                % If only one peak (no trough detected)
                if size(troughs,1) == 0
                    trlocs = find(peakCheck(pklocs:end)==0,1,'first')+pklocs;
                end
        
                % Identify start and stop of first peak
                y1 = find(-peakCheck,~0,'first');
                if size(trlocs,1) == 0
                    trlocs = y1;
                end

                % expand search in y-direction slightly
                y2 = trlocs(1)+1;

                % scan through x-direction around xCh mid location.. put
                % pixel=1 mark near troughs at each y
                flagPk = 0;
                idx=0;
                for y=y1:y2
                    
                    % Use peak closest to last center..
                    
                    % Find cell flu peak at each y..
                    % 16000 height, 500 prominence
                    wid = 30;
                    cellCheck = smoothdata(double(B(y,x-wid:x+wid)),'gaussian',[3 3]);
                    [pksCheck, pklocsCheck] = findpeaks(cellCheck,'MinPeakHeight',5000,'MinPeakProminence',500,'MinPeakDistance',5);
                    
                    % Check the amount of cells on this y
                    if size(pksCheck,2) > 1
                        [pks, pklocs] = findpeaks(cellCheck,'MinPeakHeight',16000,'MinPeakProminence',500,'MinPeakDistance',5);
                    else
                        [pks, pklocs] = findpeaks(cellCheck,'MinPeakHeight',3000,'MinPeakProminence',500,'MinPeakDistance',5);
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
                        end

                        % If more than 1 cell, keep choice constant and
                        % change segmentation depending on which side
                        % the cell is on.
                        if size(idx,1) > 1
                            xEnd = trlocs(1)+xCh(i)-(wid+1);
                        end
                        
                        % Find first and last location from peak where value is less
                        % than thresh (start and end tails of peak)
                        cellLChk = find(cellCheck>pkIntense/1.5,1,'first')-(wid+1)+xCh(i);
                        if size(idx,1) == 1 && pkIntense > pkThresh
                            cellL = find(cellCheck>pkIntense/1.5,1,'first')-(wid+1)+xCh(i);
                            cellR = (61-find(flip(cellCheck,2)>pkIntense/1.5,1,'first'))-(wid+1)+xCh(i);
                            
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
                        if exist('cellL','var')
                            lastPk = round((cellL+cellR)/2);
                            maxW = size(cellL:cellR,2)+3;
                        end
                        % overlay on img.. see how segmentation looks
        
%                         C(y,cellL) = 1;
%                         C(y,cellR) = 1;
                    end
                end
                % Check for total blobs in mother column, choose largest and discard others

                BW = mask(:,xCh(i)-30:xCh(i)+30,k);
                mask(:,xCh(i)-30:xCh(i)+30,k) = bwareafilt(BW,1);
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

                if lower(uPrompt) == 'x' || lower(uPrompt) == 'p' || lower(uPrompt) == 'l'
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
                            maskEdit2(:,delX-30:delX+30,k:end) = 0;
                            dilMask2(:,delX-30:delX+30,k:end) = 0;
                        
                        % If user puts in L, remove all pixels lower than
                        % clicked spot

                        elseif lower(uPrompt) == 'l'
    
                            % Remove cell segmentation from figure for user
                            maskEdit2(:,delX-30:delX+30,k) = 0;
                            dilMask2(:,delX-30:delX+30,k) = 0;
                
                            fused = imfuse(img,ovrImg,'falsecolor');
                            
                            % Zoom in on cell
                            xlim([delX-30 delX+30]);
                            ylim([delY-50 delY+150]);
                            
                            % Get y-value of user-click
                            [remClick(1,1),remClick(1,2)] = ginputWhite(1);
                            remY = round(click(2));
                            
                            % Update mask with everything in channel above click
                            maskEdit2(1:remY+1,delX-30:delX+30,k) = maskEdit1(1:remY+1,delX-30:delX+30,k);

                            % Update dilMask2
                            se=strel('disk',6);
                            dilMask2(:,delX-30:delX+30,k)=imclose(maskEdit2(:,delX-30:delX+30,k),se);
                            perim = bwperim(dilMask2(:,delX-30:delX+30,k),4);
                            se=strel('disk',1,0);
                            dilMask2(:,delX-30:delX+30,k)=imdilate(perim,se);

                        % If user puts in P, erase the original segmentation and replace
                        % with user-added polygon
                        elseif lower(uPrompt) == 'p'
    
                            % Remove cell segmentation from figure for user
                            maskEdit2(:,delX-30:delX+30,k) = 0;
                            dilMask2(:,delX-30:delX+30,k) = 0;
                
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
    for k = 1:size(mask,3)
        fileNum = sprintf('%03d',k);
        fileName = [filePre,fileNum,'.tif'];
        imwrite(maskEdit1(:,:,k),fileName)
    end
    mask = maskEdit1;
end
















