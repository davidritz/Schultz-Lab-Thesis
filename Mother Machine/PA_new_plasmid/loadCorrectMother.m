function [AVG_TotalGREEN,AVG_TotalRED] = loadCorrectMother(correctFile)
    
    % Green channel (constit.)
    greenPath = [correctFile,'\GFP_back.tif'];
    info = imfinfo(greenPath);
    AVG_TotalGREEN = imread(greenPath,'Info', info);
    
    % Red channel (mexXY)
    redPath = [correctFile,'\RFP_back.tif'];
    info = imfinfo(redPath);
    AVG_TotalRED = imread(redPath,'Info', info);

end







