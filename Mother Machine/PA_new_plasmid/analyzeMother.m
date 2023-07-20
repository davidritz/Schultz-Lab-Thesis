function [blobsGlobal] = analyzeMother(mask,rawVid)

    %% Get each mother cell's location and length data for a FOV

    % Only a max of 8 mother cells can be in one FOV at 100x
    % blobsGlobal (y,x,instant_length,doubling_length,doubling_counter,CFP,GFP,normal GFP,instant_growth(um/min)) for each frame
    blobsGlobal = NaN(8,9,size(mask,3));
    momCount = 0;
    
    % Remove small pixel blobs from mask
    maskClose = mask;
    
    for k = 1:size(mask,3)
        connMask = bwlabel(mask(:,:,k),4);
        numCells = max(max(connMask));
    
        % Loop through number of cells
        for c = 1:numCells
            [y,x] = find(connMask == c);
            totPix = size(x,1);
    
            % Zero small things
            if totPix < 60
                for pix = 1:totPix
                    maskClose(y(pix),x(pix),k) = 0;
                end
            end
        end
    end
    
    mask = maskClose;
    %% Get mother cell characteristics for each frame
    for k = 1:size(mask,3)
    
        % Blobs assigned L to R (column major) by regionprops
        % IF blob is within 50 x-values of an already assigned blob location,
        % it is the same blob. If not, add a new blob entry
    
    
        blobsLocal = regionprops(mask(:,:,k),'Centroid','PixelList','MajorAxisLength');
        
        if size(blobsLocal,1) > 0
            for b = 1:size(blobsLocal,1)
                x = blobsLocal(b).PixelList(:,1);
                y = blobsLocal(b).PixelList(:,2);

                totPix = size(x,1);

                if k == 1
    
                    % y = 1st column
                    blobsGlobal(b,1,k) = round(blobsLocal(b).Centroid(2));
                    % x = 2nd column
                    blobsGlobal(b,2,k) = round(blobsLocal(b).Centroid(1));
                    
                    % Get length info
                    %[feretOut, LM] = bwferet(mask(:,:,k),'MaxFeretProperties');
                    len = blobsLocal(b).MajorAxisLength(1);

                    % 0.06 microns per pixel
                    %blobsGlobal(b,3,k) = feretOut.MaxDiameter(b).*0.06;
                    blobsGlobal(b,3,k) = len*0.06;
                    blobsGlobal(b,4,k) = blobsGlobal(b,3,k);
                    blobsGlobal(b,5,k) = 0;
                    blobsGlobal(b,9,k) = 0;
                    
                    % Fluorescence 

                    % Loop through each pixel that's in the blob, b
                    for pix = 1:totPix
                        
                        % Add all fluorescence values together for each blob
                        if pix == 1
                            % rawVid(y,x,frame,flu channel); GFP = 1, RFP = 2
                
                            % RFP
                            blobsGlobal(b,6,k) = rawVid(y(pix),x(pix),k,2);
                
                            % GFP
                            blobsGlobal(b,7,k) = rawVid(y(pix),x(pix),k,1);
                        else
                            % RFP
                            blobsGlobal(b,6,k) = double(rawVid(y(pix),x(pix),k,2))+double(blobsGlobal(b,6,k));
                
                            % GFP
                            blobsGlobal(b,7,k) = double(rawVid(y(pix),x(pix),k,1))+double(blobsGlobal(b,7,k));
                        end
                    end
                    
                    % Avg RFP
                    blobsGlobal(b,6,k) = blobsGlobal(b,6,k)/totPix;
                    % Avg GFP
                    blobsGlobal(b,7,k) = blobsGlobal(b,7,k)/totPix;
                    % Normalized GFP
                    blobsGlobal(b,8,k) = double(blobsGlobal(b,7,k))/double(blobsGlobal(b,6,k));

                    % Count total # of mothers in FOV on frame 1
                    momCount = momCount+1;
                else
    
                    % Check to see if this blob already exists.. Mother cells
                    % should be locked in x-direction, so can use this to check
                    [minDistance, indexOfMin] = min(min(abs(blobsGlobal(:,2,1:k-1)-blobsLocal(b).Centroid(1)),[],3));
    
                    if minDistance > 50
                        % y = 1st column
                        blobsGlobal(momCount,1,k) = round(blobsLocal(b).Centroid(2));
                        % x = 2nd column
                        blobsGlobal(momCount,2,k) = round(blobsLocal(b).Centroid(1));
                        
                        % Get length info
                        %[feretOut, LM] = bwferet(mask(:,:,k),'MaxFeretProperties');
                        len = blobsLocal(b).MajorAxisLength(1);
        
                        % 0.06 microns per pixel
                        %blobsGlobal(momCount,3,k) = feretOut.MaxDiameter(b).*0.06;
                        blobsGlobal(momCount,3,k) = len*0.06;
                        blobsGlobal(momCount,4,k) = blobsGlobal(momCount,3,k);
                        blobsGlobal(momCount,5,k) = 0;
                        blobsGlobal(momCount,9,k) = 0;
    

                        % Fluorescence 

                        % Loop through each pixel that's in the blob, b
                        for pix = 1:totPix
                            
                            % Add all fluorescence values together for each blob
                            if pix == 1
                                % rawVid(y,x,frame,flu channel); GFP = 1, RFP = 2
                    
                                % RFP
                                blobsGlobal(momCount,6,k) = rawVid(y(pix),x(pix),k,2);
                    
                                % GFP
                                blobsGlobal(momCount,7,k) = rawVid(y(pix),x(pix),k,1);
                            else
                                % RFP
                                blobsGlobal(momCount,6,k) = double(rawVid(y(pix),x(pix),k,2))+double(blobsGlobal(momCount,6,k));
                    
                                % GFP
                                blobsGlobal(momCount,7,k) = double(rawVid(y(pix),x(pix),k,1))+double(blobsGlobal(momCount,7,k));
                            end
                        end
                        
                        % Avg RFP
                        blobsGlobal(momCount,6,k) = blobsGlobal(momCount,6,k)/totPix;
                        % Avg GFP
                        blobsGlobal(momCount,7,k) = blobsGlobal(momCount,7,k)/totPix;
                        % Normalized GFP
                        blobsGlobal(momCount,8,k) = double(blobsGlobal(momCount,7,k))/double(blobsGlobal(momCount,6,k));

                        % Add to total # of mothers in FOV
                        momCount = momCount+1;
                    else
                        % y = 1st column
                        blobsGlobal(indexOfMin,1,k) = round(blobsLocal(b).Centroid(2));
                        % x = 2nd column
                        blobsGlobal(indexOfMin,2,k) = round(blobsLocal(b).Centroid(1));
                        
                        % Get length info
                        %[feretOut, LM] = bwferet(mask(:,:,k),'MaxFeretProperties');
                        len = blobsLocal(b).MajorAxisLength(1);
        
                        % 0.06 microns per pixel
                        %blobsGlobal(indexOfMin,3,k) = feretOut.MaxDiameter(b).*0.06;
                        blobsGlobal(indexOfMin,3,k) = len*0.06;
                        
                        % Determine if cell doubled from previous frame
                        [blobsGlobal,doubling] = doubleCheckMother(blobsGlobal,indexOfMin,k);
                        
                        % Update total length of cell if no doubling occurred
                        % between frame and frame-1
                        if doubling == 0
                            blobsGlobal(indexOfMin,4,k) = blobsGlobal(indexOfMin,4,k-1) + (blobsGlobal(indexOfMin,3,k) - blobsGlobal(indexOfMin,3,k-1));
                        end
                        
                        % Running count of doublings
                        blobsGlobal(indexOfMin,5,k) = sum(blobsGlobal(indexOfMin,5,k-1))+doubling;

                        % Instant growth per min
                        blobsGlobal(indexOfMin,9,k) = (blobsGlobal(indexOfMin,4,k) - blobsGlobal(indexOfMin,4,k-1))/10;

                        % Fluorescence 

                        % Loop through each pixel that's in the blob, b
                        for pix = 1:totPix
                            
                            % Add all fluorescence values together for each blob
                            if pix == 1
                                % rawVid(y,x,frame,flu channel); GFP = 1, RFP = 2
                    
                                % RFP
                                blobsGlobal(indexOfMin,6,k) = rawVid(y(pix),x(pix),k,2);
                    
                                % GFP
                                blobsGlobal(indexOfMin,7,k) = rawVid(y(pix),x(pix),k,1);
                            else
                                % RFP
                                blobsGlobal(indexOfMin,6,k) = double(rawVid(y(pix),x(pix),k,2))+double(blobsGlobal(indexOfMin,6,k));
                    
                                % GFP
                                blobsGlobal(indexOfMin,7,k) = double(rawVid(y(pix),x(pix),k,1))+double(blobsGlobal(indexOfMin,7,k));
                            end
                        end
                        
                        % Avg RFP
                        blobsGlobal(indexOfMin,6,k) = blobsGlobal(indexOfMin,6,k)/totPix;
                        % Avg GFP
                        blobsGlobal(indexOfMin,7,k) = blobsGlobal(indexOfMin,7,k)/totPix;
                        % Normalized GFP
                        blobsGlobal(indexOfMin,8,k) = double(blobsGlobal(indexOfMin,7,k))/double(blobsGlobal(indexOfMin,6,k));

                    end
                end
            end
        end
    end
    
    % Remove all rows with NaN from all frames
    idxToRemove = all(all(isnan(blobsGlobal),3),2);
    blobsGlobal(idxToRemove,:,:) = [];

end
































