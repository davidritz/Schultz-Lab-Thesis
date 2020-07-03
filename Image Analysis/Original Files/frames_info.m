function [b,g,r,s,p]=frames_info(mask,CH2,CH3,CH4) %2,3,4 -> CFP,GFP,mCherry

D=regionprops(mask,'Centroid','Area','MajorAxisLength');
s=[D.MajorAxisLength]/27;
v=[D.Centroid];
p=v(1:2:length(v)-1);

Wb=regionprops(mask,CH2,'MeanIntensity');
b=[Wb.MeanIntensity];
Wg=regionprops(mask,CH3,'MeanIntensity');
g=[Wg.MeanIntensity];
Wr=regionprops(mask,CH4,'MeanIntensity');
r=[Wr.MeanIntensity];

