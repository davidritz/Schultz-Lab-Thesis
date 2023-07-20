% [graph] = graphMother(blobsGlobal)

% blobsGlobal (y,x,instant_length,doubling_length,doubling_counter,CFP,GFP,normal GFP,instant_growth(um/min)) for each frame
%% Scatterplot of growth after drug exposure plotted against normalized mexXY right before drug exposure

% All normalized GFP signals for all cells right before drug addition (drug
% added frame 13)
%save = blobsGlobal;
blobsGlobal = blobsGlobal(:,:,1:86);

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

figure(4)
histogram(log2(L./Lminus1))
%histogram(L)
%xlim([0 12])
%% Instant length distribution


%% OG method of counting length
% Plot added growth over time for each cell
run = blobsGlobal(:,:,1:86);

for c = 1:size(run,1)
    dat1 = log2(run(c,4,:));
    dat2Smooth = dat1(:);
    run(c,4,:) = dat1;%smoothdata(dat2Smooth,'movmean',[1,1]);
end

% Re-calculate the instantaneous growth from the smoothing doubling_length

for f = 2:size(run,3)
    for c = 1:size(run,1)
        run(c,9,f) = run(c,4,f) - run(c,4,f-1);
    end
end

run(:,9,1) = run(:,9,2);


t = (0:10/60:size(run,3)*(10/60)-(10/60));
figure(5)
for c = 1:size(run,1)
    hold on
    grow = run(c,4,:);
    growthVec = grow(:);
    plot(t,growthVec)
end
hold off
%ylim([-1 3])
xlabel('Time (hr)')
ylabel('Total length of cell (um)')
title('Growth curves of cells')

%%

drugAdd = 243-123;
% blobsGlobal (y,x,instant_length,doubling_length,doubling_counter,CFP,GFP,normal GFP,instant_growth(um/min)) for each frame
logBlobs = blobsGlobal(:,:,:);

% Change instantaneous growth to log2 scale
logBlobs(:,3,:) = log2(logBlobs(:,3,:));

% Copy instant length to doubling_length column
logBlobs(:,4,:) = logBlobs(:,3,:);

doublingCount = zeros(size(logBlobs,1),1);
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
    dat1 = logBlobs(c,4,1:drugAdd);
    dat2 = logBlobs(c,4,drugAdd+1:end);
    logBlobs(c,4,1:drugAdd) = smoothdata(dat1,'movmean',[1,1],'omitnan');
    logBlobs(c,4,drugAdd+1:end) = smoothdata(dat2,'movmean',[1,1],'omitnan');
end

% Re-calculate the instantaneous growth from the smoothing doubling_length

for f = 2:size(logBlobs,3)
    for c = 1:size(logBlobs,1)
        logBlobs(c,9,f) = (logBlobs(c,4,f) - logBlobs(c,4,f-1)).*6;
    end
end

for c = 1:size(logBlobs,1)
    dat1 = logBlobs(c,9,1:drugAdd);
    dat2 = logBlobs(c,9,drugAdd+1:end);
    logBlobs(c,9,1:drugAdd) = smoothdata(dat1,'movmean',[3,3],'omitnan');
    logBlobs(c,9,drugAdd+1:end) = smoothdata(dat2,'movmean',[3,3],'omitnan');
end

% Don't use any timepoints from smoothed curve where cell isn't detected
for f = 2:size(logBlobs,3)
    for c = 1:size(logBlobs,1)
        if isnan(logBlobs(c,3,f))
            logBlobs(c,9,f) = NaN;
        end
    end
end

% Plot the instant growth of each cell against time
% Right now, units of growth are doub. per 10 min.. multiply by 6
% x doub./frame * frame/10 min * 60 min/hr

t = (0:5/60:size(logBlobs,3)*(5/60)-(5/60));
drugAdd = 5/60*(243-123);
t = t-drugAdd;
figure(17)
%for c = 1:size(logBlobs,1)
for c = 1:size(logBlobs,1)
    hold on
    grow = logBlobs(c,9,2:end);
    growthVec = grow(:);
    plot(t(2:end),growthVec)
end

ylim([-0.2 1.2])
xlim([-4 12])
y = ylim;
plot([0 0],[y(1) y(2)],'--k','LineWidth',3)
hold off

xlabel('Time (hr)')
ylabel('Growth (doub./hr)')
title('Growth curves of cells')


set(findall(gcf,'-property','FontSize'),'FontSize',26)
%%
yVal = [];
xVal = [];
for c = 1:size(logBlobs,1)
    dat = squeeze(logBlobs(c,9,2:end));
    dat(isnan(dat)) = 0;
    datT = squeeze(t(2:end));
    yVal = [yVal;dat];
    xVal = [xVal;datT'];
end

figure
heatScat = heatscatter(xVal, yVal, [], 'test.png','10','80');
ylim([-1 3])
%% plot normalized flu vs time for mut and wt

dat1 = wDat(:,:,1:86);
dat2 = mDat(:,:,1:86);

figure(21)
t = (0:10/60:size(logBlobs,3)*(10/60)-(10/60));

% drug added at frame 13
drugAdd = 10/60*(13);

t = t-drugAdd;

hold on

% dat1Mean = mean(dat1(:,8,:),'omitnan');
% dat1Mean = dat1Mean(:);
% plot(t,dat1Mean,'LineWidth',8,'Color',[0.9290 0.6940 0.1250])
% 
dat2Mean = mean(dat2(:,8,:),'omitnan');
dat2Mean = dat2Mean(:);
plot(t,dat2Mean,'LineWidth',8,'Color',[0.4940 0.1840 0.5560])

ylim([1 3])
xlim([-2.1 12])
y=ylim;
plot([0 0],[y(1) y(2)],'--k','LineWidth',3)

%lgd = legend('WT PA14','Δ{\itmexZ} PA14','700 μg/mL spec. (maintained)','Location','northwest');

% for c = 1:size(dat2,1)
%     normFlu = dat1(c,8,:);
%     fluVec = smoothdata(normFlu(:),'movmean',[5 5]);
% %     %fluVec = normFlu(:);
% %     if c > 1
% %         hline = plot(t,fluVec,'LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
% %     else
% %         hline = plot(t,fluVec,'LineWidth',2,'HandleVisibility','off','Color',[0.9290 0.6940 0.1250]);
% %     end
%     hline = plot(t,fluVec,'LineWidth',2,'Color',[0.9290 0.6940 0.1250],'HandleVisibility','off');
%     hline.Color = [hline.Color 0.35];  % alpha=0.1
% end

for c = 1:size(dat2,1)
    normFlu = dat2(c,8,:);
    fluVec = smoothdata(normFlu(:),'movmean',[5 5]);
    %fluVec = normFlu(:);
    hline = plot(t,fluVec,'LineWidth',2,'Color',[0.4940 0.1840 0.5560],'HandleVisibility','off');
    hline.Color = [hline.Color 0.35];  % alpha=0.1
end

%lgd = legend('700 μg/mL spec. (maintained)','WT PA14 cell','Location','northwest');
ylim([1 3])
xlim([-2.1 12])

grid on;
ylabel('Norm. MexXY')
xlabel('Time (hrs)')
set(findall(gcf,'-property','FontSize'),'FontSize',26)
hold off

%%







