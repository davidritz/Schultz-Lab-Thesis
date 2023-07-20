% Something a good noodle should never do

% Adapted from https://github.com/SystemsBiologyUniandes/MotherMachineAnalysis/blob/master/correction.ipynb

% Inputted vid should only be 3D i.e. vid(:,:,:,2) --> 2 is GFP flu channel for PA
function [correctSeg,correctFlu] = correctMother(AVG_TotalGREEN,rawVid)

    % Get size data of each frame and number of frames
    numY = size(rawVid,1);
    numX = size(rawVid,2);
    numFrames = size(rawVid,3);
    
    %AVG_TotalGREEN = uint64(AVG_TotalGREEN);
    correctFlu = uint16(zeros(numY,numX,numFrames));
    correctSeg = uint16(zeros(numY,numX,numFrames));

    maxi = mean(mean(AVG_TotalGREEN));


    for frame = 1:numFrames

        img =rawVid(:,:,frame);

        % Applying FFC to GFP images
        for y = 1:numY
            for x = 1:numX
                if img(y,x) > AVG_TotalGREEN(y,x)
                    correctFlu(y,x,frame) = maxi * (img(y,x) - AVG_TotalGREEN(y,x)) / AVG_TotalGREEN(y,x);
                end
            end
        end
    
        % Pixels too bright are assigned to uint16 max, pixels too dim are asigned to 0, and a corresponding 
        % scaled value in between 
        
        % ADD 6/22/23
        % Dynamic threshold for each column of cells
        % Don't look at media channel
        A = correctFlu(1:800,:,frame);
        
        % Find x-location of channels
        winWidth = 20;
        xArr=sum(double(A),1);
        [bla,xCh]=findpeaks(xArr,'MinPeakDistance',40,'MinPeakHeight',10000,'MinPeakProminence',10000);
        for i = 1:size(xCh,2)
            c = max(max(correctFlu(:,xCh(i)-winWidth:xCh(i)+winWidth,frame)));
            % 4.8, 0.99 6/23
            %background = 4*mean(mean(correctFlu(:,xCh(i)-winWidth:xCh(i)+winWidth,frame)));
            background = 110;
            maxThreshold = 0.99 * c;
            jpegMax = 65535;
            for y = 1:numY
                for x = xCh(i)-winWidth:xCh(i)+winWidth
                    if 0.5 * c > maxThreshold 
                        if correctFlu(y,x,frame) < background 
                            correctSeg(y,x,frame) = 0;
                        elseif correctFlu(y,x,frame) > maxThreshold
                            correctSeg(y,x,frame) = jpegMax;
                        else
                            %correctSeg(y,x,frame) = jpegMax * ((correctFlu(y,x,frame) - background) / (maxThreshold-background));
                            correctSeg(y,x,frame) = jpegMax * ((double(correctFlu(y,x,frame) - background)) / double(maxThreshold-background));
                        end
                            
                    else
                        
                        if correctFlu(y,x,frame) < background 
                            correctSeg(y,x,frame) = 0;
                        elseif correctFlu(y,x,frame) > maxThreshold
                            correctSeg(y,x,frame)=jpegMax;
                        else
                            correctSeg(y,x,frame) = jpegMax * ((double(correctFlu(y,x,frame) - background)) / double(maxThreshold-background));
                        end
                    end
        
                end
            end
        end
        if mod(frame,100) == 0
            disp(['Just applied FFC on frame ',num2str(frame)])
        end
    end
end


