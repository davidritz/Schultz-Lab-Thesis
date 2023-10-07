% blobsGlobal (y,x,instant_length,doubling_length,doubling_counter,RFP,GFP,normal GFP,instant_growth(um/min)) 
% for each frame

% Load datasets
load('blobsGlobal_wt_008.mat')
wt=blobsGlobal;
load('blobs_KOmexZ_fixFOV10_12hr_zero_24_30.mat')
koz = blobsGlobal;
load('blobsGlobal_KOmexY8_008.mat')
koy = blobsGlobal;

% Take only frames I segmented
wt_plot=wt(:,:,1:204);
koz_plot=koz(:,:,1:213);
koy_plot=koy(:,:,1:192);

% Drug addition times (frames)
drug_wt = 60;
drug_koy = 120;
drug_koz = 67;
%% Split wt into growing and not growing
clear wt_ng
clear wt_g
% threshold for final growth.. splits wt into growing and not growing
thresh = 0.1;
ng = 1;
g = 1;

logBlobs = conv2Log(wt_plot);

for c = 1:size(logBlobs,1)
    if isnan(logBlobs(c,9,size(logBlobs,3))) || logBlobs(c,9,size(logBlobs,3)) < thresh
        wt_ng(ng,:,:) = wt_plot(c,:,:);
        ng=ng+1;
    else
        wt_g(g,:,:) = wt_plot(c,:,:);
        g=g+1;
    end
end
%% Plot growth all together
[vec1,t1] = graphGrowExp(wt_g,drug_wt,[0 0 0.61, 0.2],1);
[vec2,t2] = graphGrowExp(wt_ng,drug_wt,[0 0.75 1, 0.2],1);
[vec3,t3] = graphGrowExp(koz_plot,drug_koz,[0.8 0.33 0 0.2],2);
[vec4,t4] = graphGrowExp(koy_plot,drug_koy,[1 0.75 0 0.2],3);
avg1 = plot(t1,vec1,'Color',[0 0 0.61],'LineWidth',6);
avg2 = plot(t2,vec2,'Color',[0 0.75 1],'LineWidth',6);
avg3 = plot(t3,vec3,'Color',[0.8 0.33 0],'LineWidth',6);
avg4 = plot(t4,vec4,'Color',[1 0.75 0],'LineWidth',6);
%ylim = y;
%ymax = max(y);
%ymin = min(y);
%yline(0,'r--','LineWidth',4)
plot([t4(end) 12], [avg4(end) avg4(end)],'--','Color',[1 0.75 0],'LineWidth',6)
legend([avg1,avg2,avg3,avg4],'WT growing','WT arrested','Δ{\itmexZ}','Δ{\itmexY}')
set(findall(gcf,'-property','FontSize'),'FontSize',26)
hold off






