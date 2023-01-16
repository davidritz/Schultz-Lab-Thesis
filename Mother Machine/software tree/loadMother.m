function [vid, outName] = loadMother(parFile,fovNum)
    
    % e.g. user inputs FOV 2, will grab 3rd file and 4th file in folder
    fovNumCorr = fovNum*2-1;

    % cd into directory and order files by alphabetical order
    cd(parFile)
    parDir = dir;
    [~,ind]=sort({parDir.name});
    parDir = parDir(ind);
    
    % Remove any blank files from dir variable
    blank = 0;
    for fov=1:size(parDir,1)
        fov = fov-blank;
        if parDir(fov).bytes == 0
            parDir(fov) = [];
            blank=blank+1;
        end
    end
    
    % vid(y,x,frame,flu channel)
    
    fovTot = size(parDir,1)/2;
    %%
    for fov=fovNumCorr:fovNumCorr+1
        fname = [parFile,'\',parDir(fov).name];
        info = imfinfo(fname);
        
        % Get tif data for single frame
        if fov == fovNum
            p={info.Width};
            p=p{1};
            q={info.Height};
            q=q{1};
            frames = size(info,1);
            vid = uint16(zeros(q,p,frames,2));
        end

        for f = 1:frames
            currentImg = imread(fname,f,'Info', info);
            % If CFP
            if fov == fovNum
                vid(:,:,f,1) = currentImg;
            % Else GFP
            else 
                vid(:,:,f,2) = currentImg;
            end
        end

        % Output progress to user
        if fov == fovNum
            disp(['Just loaded GFP FOV ',num2str(fovNum)])
        else
            disp(['Just loaded CFP FOV ',num2str(fovNum)])
        end
    end
    outName = parDir(fov).name;
end







