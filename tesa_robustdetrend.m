% tesa_robustdetrend() - performs robust polynomial detrending on EEG data 
%                       within a specified time window, optionally excluding a sub-window
%                       (e.g., to avoid fitting over TMS artefacts or tep components). 
%                       It removes slow trends from the signal while reducing the influence 
%                       of outliers.
%
% Usage:
%   >>  EEG = tesa_robustdetrend(EEG, timeWindow, thresholdSD, polyOrder, excludeWindow);
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   timeWindow      - vector (required). Two-element vector [startTime, endTime] 
%                   defining the segment to detrend in each trial. Time is in milliseconds.
%                   Example: [-1000, 1000]
%   thresholdSD     - integer (required). Standard deviation threshold for
%                   outlier exclusion.
%                   Example: 3
%   polyOrder       - integer (required). Order of the polynomial used for detrending (e.g., 1 = linear, 2 = quadratic).
%                   Example: 1 (linear detrending)
%   excludeWindow   - vectpr (optional) Two-element vector [start end] defining the
%                   time range to exclude from the polynomial fit. Leave empty if unused.
%                   Time is in milliseconds.
%                   Example: [-10,500]
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = tesa_robustdetrend(EEG, [-1000,1000], 3, 1); % performs a linear detrend between -1000 to 1000 ms not including outlier data > 3 standard deviations from the mean.
%   EEG = tesa_robustdetrend(EEG, [-1000,0], 5, 2); % performs a quadratic detrend between -1000 to 0 ms not including outlier data > 5 standard deviations from the mean.
%   EEG = tesa_robustdetrend(EEG, [-1000,1000], 3, 1, [-10,500]); % performs a linear detrend between -1000 to 1000 ms not including outlier data > 3 standard deviations from the mean and excluding data between -10 and 500 ms.
%
% 
% See also:
%   tesa_detrend

% Copyright (C) 2025  Nigel Rogasch, University of Adelaide,
% nigel.rogasch@adelaide.edu.au
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% Author: Nigel Rogasch and ChatGPT (2025), based on Nigel Rogasch's TESA framework conventions.

function EEG = tesa_robustdetrend(EEG,timeWindow,thresholdSD,polyOrder,excludeWindow)

    if nargin < 4
        error('Not enough input arguments.');
    end
    if nargin < 5
        excludeWindow = [];
    end

    % Check EEG structure
    if ~isstruct(EEG) || ~isfield(EEG, 'data') || ~isfield(EEG, 'times')
        error('EEG must be a structure containing fields "data" and "times".');
    end

    % Get limits of EEG time vector
    xmin = min(EEG.times);
    xmax = max(EEG.times);

    % Validate timeWindow
    if ~isnumeric(timeWindow) || numel(timeWindow) ~= 2
        error('timeWindow must be a two-element numeric vector, e.g., [-100 500].');
    elseif timeWindow(1) >= timeWindow(2)
        error('timeWindow(1) must be less than timeWindow(2).');
    elseif timeWindow(1) < xmin || timeWindow(2) > xmax
        error('timeWindow values must fall within the EEG time range [%.2f, %.2f].', xmin, xmax);
    end

    % Validate thresholdSD
    if ~isnumeric(thresholdSD) || ~isscalar(thresholdSD) || thresholdSD <= 0
        error('thresholdSD must be a positive numeric scalar.');
    end

    % Validate polyOrder
    if ~isnumeric(polyOrder) || ~isscalar(polyOrder) || polyOrder < 0 || mod(polyOrder,1) ~= 0
        error('polyOrder must be a non-negative integer scalar.');
    end

    % Validate excludeWindow (if provided)
    if ~isempty(excludeWindow)
        if ~isnumeric(excludeWindow) || numel(excludeWindow) ~= 2
            error('excludeWindow must be a two-element numeric vector, e.g., [-10 500].');
        elseif excludeWindow(1) >= excludeWindow(2)
            error('excludeWindow(1) must be less than excludeWindow(2).');
        elseif excludeWindow(1) < timeWindow(1) || excludeWindow(2) > timeWindow(2)
            error('excludeWindow must fall entirely within timeWindow.');
        end
    end

    % -----------------------
    % Detrending process
    % -----------------------

    % Find indices for timeWindow
    [~, t1] = min(abs(timeWindow(1) - EEG.times));
    [~, t2] = min(abs(timeWindow(2) - EEG.times));

    % Extract data for the time window
    eegData = EEG.data(:, t1:t2, :);
    [nChannels, nTimePoints, nTrials] = size(eegData);

    % Time vector for fitting
    timeVec = EEG.times(t1:t2);
    detrendedData = zeros(size(eegData));

    % Determine exclusion indices (if any)
    if ~isempty(excludeWindow)
        excludeIdx = EEG.times(t1:t2) >= excludeWindow(1) & EEG.times(t1:t2) <= excludeWindow(2);
    else
        excludeIdx = false(size(timeVec));
    end

    % Loop over channels and trials
    for ch = 1:nChannels
        for tr = 1:nTrials
            y = squeeze(eegData(ch, :, tr));  % 1 x time

            % Step 1: Exclude time points from excludeWindow
            validIdx = ~excludeIdx;

            % Step 2: Estimate initial mean & std from valid points
            mu = mean(y(validIdx));
            sigma = std(y(validIdx));

            % Step 3: Identify inliers within ±thresholdSD
            inlierIdx = validIdx & abs(y - mu) <= thresholdSD * sigma;

            % Step 4: Fit polynomial to inliers only
            t_inliers = timeVec(inlierIdx);
            y_inliers = y(inlierIdx);

            p = polyfit(t_inliers, y_inliers, polyOrder);

            % Step 5: Evaluate and subtract polynomial trend
            trend = polyval(p, timeVec);
            detrendedData(ch, :, tr) = y - trend;
        end
    end

    % Replace detrended segment
    EEG.data(:, t1:t2, :) = detrendedData;
end