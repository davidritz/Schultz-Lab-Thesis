% [graph] = graphMother(blobsGlobal)

% blobsGlobal (y,x,instant_length,doubling_length,doubling_counter,CFP,GFP,normal GFP,instant_growth(um/min)) for each frame
%% Scatterplot of growth after drug exposure plotted against normalized mexXY right before drug exposure

% All normalized GFP signals for all cells right before drug addition (drug
% added frame 13)
%save = blobsGlobal;
blobsGlobal = blobsGlobal(:,:,1:85);

x=blobsGlobal(:,8,12);
x=x(:);

% Total growth of cells after drug addition. Get final value (some cells
% die/arrest/aren't detectable on final frame)

y = zeros(size(blobsGlobal,1),1);
numDiv = zeros(size(blobsGlobal,1),1);
instGrow = zeros(size(blobsGlobal,1),1);
finalMexXY = zeros(size(blobsGlobal,1),1);
totalLen = zeros(size(blobsGlobal,1),size(blobsGlobal,3));

for c = 1:size(blobsGlobal,1)

    % Get all total growth values of cell, c
    arrayOfCells = blobsGlobal(c,4,:);
    arrayOfCells = arrayOfCells(:);

    arrayOfDiv = blobsGlobal(c,5,:);
    arrayOfDiv = arrayOfDiv(:);

    arrayOfGrow = blobsGlobal(c,9,:);
    arrayOfGrow = arrayOfGrow(:);

    arrayOfMexXY = blobsGlobal(c,8,:);
    arrayOfMexXY = arrayOfMexXY(:);
    
    
    totalLen(c,:) = arrayOfCells;

    % Remove all NaN values
    indRmv = ~isnan(arrayOfCells);
    arrayOfCells = arrayOfCells(indRmv);

    indRmv = ~isnan(arrayOfDiv);
    arrayOfDiv = arrayOfDiv(indRmv);

    indRmv = ~isnan(arrayOfMexXY);
    arrayOfMexXY = arrayOfMexXY(indRmv);
    
    % Record number at end minus right before drug addition
    y(c) = arrayOfCells(end)-arrayOfCells(12);
    numDiv(c) = arrayOfDiv(end)-arrayOfDiv(12);
    % log2 of dL/dt growth
    instGrow(c) = (arrayOfGrow(end)-arrayOfGrow(end-2))/30;
    finalMexXY(c) = arrayOfMexXY(end);
end

%%
figure(1)
cmapSize = max(numDiv)+1;
cmap = hot(cmapSize);
numValues = length(numDiv);
markerColors = zeros(numValues, 3);
% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = numDiv(k)+1;
    markerColors(k, :) = cmap(row, :);
end
% Create the scatter plot.
scatter(x, y, 60, markerColors,'filled','MarkerEdgeColor','black');
c = colorbar;
colormap(hot(cmapSize));
clim([0 max(numDiv)])
grid on;

hold on
h = plot(stats);

ylabel('Total growth after drug exposure (μm)')
xlabel('MexXY (norm.) just before drug exposure')
ylabel(c,'Number of divisions after drug exposure')
title('MexXY expression upon drug exposure predicting cell growth of PA14')

% change linear regression colors

% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
% The confidence bounds have 2 handles but only one of 
% the handles contains the legend string.  The first
% line below finds that object and then searches for 
% other objects in the plot that have the same linestyle
% and color. 
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);

% Grab 95% CI data
xCI = cbHandles(1).XData;
yL = cbHandles(1).YData;
yU = cbHandles(2).YData;

% Hide data and automatically generated 95% CI
set([dataHandle;dataHandle.Children], 'Visible','off')
set([cbHandles;cbHandles.Children], 'Visible','off')

% Plot 95% CI
set(fitHandle, 'Color', [0 0.5 0], 'LineWidth', 1)
patch([xCI fliplr(xCI)], [yU fliplr(yL)], 'g', 'facealpha',0.3,'LineStyle','none')

r = round(corr(x,y),3);
p = round(coefTest(stats),4);
lgd = legend('Cell data points','','Trend line','','','95% CI');
xLoc = 2;
yLoc = 14;
text(xLoc,yLoc,['Pearson''s r = ',num2str(r)])
text(xLoc+0.066,yLoc-2,['p-value = ',num2str(p)])

set(findall(gcf,'-property','FontSize'),'FontSize',18)
hold off

%% Bin data according to total growth after drug exposure, plot normalized fluorescence over time

%% Plot final growth rate in doublings per hour by final mexXY (norm)
ind = isnan(instGrow);
instGrow(ind) = 0;

plot(finalMexXY,instGrow,'o')

ylabel('d/dx of growth at end of experiment (last 3 frames) (μm)')
xlabel('Last recorded MexXY (norm.) of cell')
title('Final growth and MexXY (norm.) of cells')




%%
% figure(1)
% g = gramm('x',x,'y',y,'color',numDiv);
% g.geom_point();
% %g.geom_line();
% %g.stat_glm();
% g.set_names('x','Normalized MexXY right before drug exposure','y','Total growth of cell (μm)','color','# of divisions')
% g.draw()
% %set([g.results.geom_point_handle],'MarkerSize',2);
% hold on 
% plot(1,1)
% hold off

%%

Lminus1 = zeros((size(blobsGlobal,3)-1)*size(blobsGlobal,1),1);
L = zeros((size(blobsGlobal,3)-1)*size(blobsGlobal,1),1);

count = 1;
for c = 1:size(blobsGlobal,1)
    for f = 2:size(blobsGlobal,3)
        L(count) = blobsGlobal(c,3,f);
        Lminus1(count) = blobsGlobal(c,3,f-1);
        count = count+1;
    end
end

figure(3)
histogram(log2(L./Lminus1))
%%
% Plot added growth over time for each cell
t = (0:10/60:size(totalLen,2)*(10/60)-(10/60));
figure(5)
for p = 1:size(totalLen,1)
    hold on
    plot(t,(totalLen(p,:)-totalLen(p,1)))
end
hold off

xlabel('Time (hr)')
ylabel('Total growth of cell (μm)')
title('Growth curves of cells')












