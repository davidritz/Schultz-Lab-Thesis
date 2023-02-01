function [blobsGlobalNew,doubling] = doubleCheckMother(blobsGlobal,idx,frame)

    doubFac = 0.7;
    
    currentLen = blobsGlobal(idx,3,frame);
    lastLen = blobsGlobal(idx,3,frame-1);
    
    blobsGlobalNew = blobsGlobal;
    
    % If instantaneous length column decreased enough from frame-1 to frame, count it as doubling
    if lastLen*doubFac > currentLen
    
        % Reflect instant length + doubling in doubling_length column
        blobsGlobalNew(idx,4,frame) = blobsGlobalNew(idx,4,frame-1) + currentLen - (lastLen/2);
    
        doubling = 1;
    else
        doubling = 0;
    end
end





