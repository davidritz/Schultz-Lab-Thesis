function [AVG_Dark,AVG_TotalGREEN,AVG_TotalRED] = loadCorrectMother(correctFile)
    
    % Average (Green+red)
    darkPath = [correctFile,'\AVG_Dark.tif'];
    info = imfinfo(darkPath);
    AVG_Dark = imread(darkPath,'Info', info);
    
    % Green channel (mexXY)
    greenPath = [correctFile,'\AVG_TotalGREEN.tif'];
    info = imfinfo(greenPath);
    AVG_TotalGREEN = imread(greenPath,'Info', info);
    
    % Red channel (constitutive)
    redPath = [correctFile,'\AVG_TotalRED.tif'];
    info = imfinfo(redPath);
    AVG_TotalRED = imread(redPath,'Info', info);

end







