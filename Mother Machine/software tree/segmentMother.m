function [mask] = segmentMother(correctSeg,correctFlu,maskLoc,outName)

    %% Turn off findpeaks warning, so it doesn't flood the cmd window later
    findpeaks([1 1 1 1 ]);
    [msg,id] = lastwarn;
    warning('off',id)
    close all
    
    %% Auto tracking
    mask=false(size(correctSeg));
    for k=1:size(correctSeg,3)
        % Don't look at media channel
        B = correctSeg(:,:,k);
        A = correctSeg(1:800,:,k);
        
        % Find x-location of channels
        xArr=sum(double(A),1);
        [bla,xCh]=findpeaks(xArr,'MinPeakDistance',40,'MinPeakHeight',100000,'MinPeakProminence',50000);
    
        % Search in python mask
        % C = mat2gray(B);
        % Scan through mother channels
        for i=1:size(xCh,2)
            x=xCh(i);
            peakCheck = smoothdata(double(sum(B(:,x-30:x+30),2)),'gaussian',[6,6]);
            thresh1 = 5000;
            thresh2 = 5000;
            [troughs,trlocs] = findpeaks(-peakCheck,'MinPeakProminence',thresh2);
            [pks, pklocs] = findpeaks(peakCheck,'MinPeakDistance',10,'MinPeakHeight',thresh1,'MinPeakProminence',thresh2);
            
            % If only one peak (no trough detected)
            if size(troughs,1) == 0
                trlocs = find(peakCheck(pklocs:end)==0,1,'first')+pklocs;
            end
    
            % Identify start and stop of first peak
            y1 = find(-peakCheck,~0,'first');
            y2 = trlocs(1);
        
            % scan through x-direction around xCh mid location.. put
            % pixel=1 mark near troughs at each y
            for y=y1:y2
                
                % Find cell flu peak at each y..
                cellCheck = smoothdata(double(B(y,x-30:x+30)),'gaussian',[5 5]);
                [pks, pklocs] = findpeaks(cellCheck,'MinPeakHeight',2000,'MinPeakProminence',1000,'MinPeakDistance',20);
    
                if size(pks,2) ~= 0
                    
                    % Find first and last location from peak where value is less
                    % than thresh (start and end tails of peak)
    
                    cellL = find(cellCheck>pks/1.5,1,'first')-31+xCh(i);
                    cellR = (61-find(flip(cellCheck,2)>pks/1.5,1,'first'))-31+xCh(i);
                    mask(y,cellL:cellR,k) = 1;
                    mask(y,cellR,k) = 1;
    
                    % overlay on img.. see how segmentation looks
    
%                     C(y,cellL) = 1;
%                     C(y,cellR) = 1;
                end
            end
        end
%         figure(1)
%         imshow(C)
%     
%         figure(2)
%         imshow(mask(:,:,k))
    end
    %% Mask user correction
    
    % Dilate mask slightly, so user can see segmentation easier
    se=strel('disk',2);
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
    
    % show mask starting at frame=1, prompt user if edit needed, then run below
    k=1;
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
                if lower(uPrompt) == 'x' || lower(uPrompt) == 'p'
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
        fileName = filePre+fileNum+'.tif';
        imwrite(maskEdit1(:,:,k),fileName)
    end
    mask = maskEdit1;
end
















