function logBlobs = conv2Log(blobsGlobal)

    % blobsGlobal (y,x,instant_length,doubling_length,doubling_counter,RFP,GFP,normal GFP,instant_growth(um/min)) for each frame
    logBlobs = blobsGlobal(:,:,1:size(blobsGlobal,3));
    
    % Change instantaneous growth to log2 scale
    logBlobs(:,3,:) = log2(logBlobs(:,3,:));
    
    % Copy instant length to doubling_length column
    logBlobs(:,4,:) = logBlobs(:,3,:);
    
    doublingCount = zeros(size(logBlobs,1),1);
    %drugAdd = 60;
    
    % Loop over each cell.. add 1 to every frame in doubling_length column for a cell row once
    % a doubling happens
    for f = 2:size(logBlobs,3)
        for c = 1:size(logBlobs,1)
    
            % If a doubling is recorded for cell
            if logBlobs(c,5,f) > doublingCount(c)
    
                % Add 1 to row c, column 4 from doubling --> end
                logBlobs(c,4,f:end) = logBlobs(c,4,f:end)+1;
            
                % Update cell doubling-counting vector
                doublingCount(c) = logBlobs(c,5,f);
            end
        end
    end
    
    % Smooth the doubling_length
    
    for c = 1:size(logBlobs,1)
        dat1 = logBlobs(c,4,:);
        %dat2 = logBlobs(c,4,drugAdd+1:end);
        logBlobs(c,4,:) = smoothdata(dat1,'movmean',[3,3],'omitnan');
        %logBlobs(c,4,drugAdd+1:end) = smoothdata(dat2,'movmean',[1,1],'omitnan');
    end
    
    % Re-calculate the instantaneous growth from the smoothing doubling_length
    
    for f = 2:size(logBlobs,3)
        for c = 1:size(logBlobs,1)
            logBlobs(c,9,f) = (logBlobs(c,4,f) - logBlobs(c,4,f-1)).*6;
        end
    end
    
    for c = 1:size(logBlobs,1)
        dat1 = logBlobs(c,9,:);
        %dat2 = logBlobs(c,9,drugAdd+1:end);
        logBlobs(c,9,:) = smoothdata(dat1,'movmean',[5,3],'omitnan');
        %logBlobs(c,9,drugAdd+1:end) = smoothdata(dat2,'movmean',[5,5],'omitnan');
    end
    
    % Don't use any timepoints from smoothed curve where cell isn't detected
    for f = 2:size(logBlobs,3)
        for c = 1:size(logBlobs,1)
            if isnan(logBlobs(c,3,f))
                logBlobs(c,9,f) = NaN;
            end
        end
    end
end