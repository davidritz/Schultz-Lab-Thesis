function data=filter_spikes(data0)

% % Define the window size for the moving median
% windowSize = 3;
% 
% % Calculate the moving median
% filtered = medfilt1(data0, windowSize);
% 
% % Detect spikes: Find where the absolute difference is too large
% threshold = 3; % Threshold for spike detection, adjust based on data characteristics
% spikeIdx = abs(data0 - filtered) > threshold * std(data0);
% 
% % Replace spikes with the median filtered value
% data0(spikeIdx) = filtered(spikeIdx);
% data=data0;

% Define the threshold for change
threshold = 0.05;

% Length of data vector
n = length(data0);

% Iterate over the data vector to find and adjust for large changes
for i = 1:n-1
    if abs(data0(i+1) - data0(i)) > threshold
        % Calculate the shift
        shift = data0(i+1) - data0(i);
        
        % Adjust the rest of the vector to disconsider the shift
        data0(i+1:end) = data0(i+1:end) - shift;
    end
end
data=data0;
end