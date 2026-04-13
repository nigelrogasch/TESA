% tesa_robustdemean - Removes the mean from EEG data in a time window using
%                     robust averaging (excluding outliers) and optional exclusion
%                     of a subwindow.
%
% USAGE:
%   EEG = tesa_robustdemean(EEG, timeWindow, thresholdSD)
%   EEG = tesa_robustdemean(EEG, timeWindow, thresholdSD, excludeWindow)
%
% INPUTS:
%   EEG          - EEGLAB EEG structure containing fields 'data' and 'times'
%   timeWindow   - Two-element vector [tmin tmax] defining the time range (ms)
%   thresholdSD  - Outlier threshold in standard deviations (e.g., 3)
%   excludeWindow - (Optional) Two-element vector [tmin tmax] defining a region
%                   to exclude from the mean calculation (ms)
%
% OUTPUT:
%   EEG - EEG structure with demeaned data in the specified time window
%
% Example:
%   EEG = tesa_robustdemean(EEG, [-500 500], 3, [-50 100]);

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

function EEG = tesa_robustdemean(EEG, timeWindow, thresholdSD, excludeWindow)

% -------------------------------------------------------------------------
% Validate inputs
% -------------------------------------------------------------------------
if nargin < 3
    error('Not enough input arguments. Usage: tesa_robustdemean(EEG, timeWindow, thresholdSD, [excludeWindow])');
end

% Check EEG structure
if ~isstruct(EEG) || ~isfield(EEG, 'data') || ~isfield(EEG, 'times')
    error('EEG must be a structure containing fields "data" and "times".');
end

xmin = min(EEG.times);
xmax = max(EEG.times);

% timeWindow
if ~isnumeric(timeWindow) || numel(timeWindow) ~= 2
    error('timeWindow must be a two-element numeric vector, e.g., [-100 500].');
elseif timeWindow(1) >= timeWindow(2)
    error('timeWindow(1) must be less than timeWindow(2).');
elseif timeWindow(1) < xmin || timeWindow(2) > xmax
    error('timeWindow values must fall within the EEG time range [%.2f, %.2f].', xmin, xmax);
end

% thresholdSD
if ~isnumeric(thresholdSD) || ~isscalar(thresholdSD) || thresholdSD <= 0
    error('thresholdSD must be a positive numeric scalar.');
end

% excludeWindow
if nargin < 4 || isempty(excludeWindow)
    excludeWindow = [];
elseif ~isnumeric(excludeWindow) || numel(excludeWindow) ~= 2
    error('excludeWindow must be a two-element numeric vector, e.g., [-10 100].');
elseif excludeWindow(1) >= excludeWindow(2)
    error('excludeWindow(1) must be less than excludeWindow(2).');
elseif excludeWindow(1) < xmin || excludeWindow(2) > xmax
    error('excludeWindow values must fall within the EEG time range [%.2f, %.2f].', xmin, xmax);
end

% -------------------------------------------------------------------------
% Extract time indices
% -------------------------------------------------------------------------
[~, t1] = min(abs(timeWindow(1) - EEG.times));
[~, t2] = min(abs(timeWindow(2) - EEG.times));

includeIdx = true(size(EEG.times));
includeIdx(1:t1-1) = false;
includeIdx(t2+1:end) = false;

if ~isempty(excludeWindow)
    [~, e1] = min(abs(excludeWindow(1) - EEG.times));
    [~, e2] = min(abs(excludeWindow(2) - EEG.times));
    includeIdx(e1:e2) = false;
end

% -------------------------------------------------------------------------
% Robust demean per channel and trial
% -------------------------------------------------------------------------
eegData = EEG.data(:, :, :);
[nChannels, ~, nTrials] = size(eegData);

for ch = 1:nChannels
    for tr = 1:nTrials
        y = squeeze(eegData(ch, :, tr));

        % Use only included time points
        y_included = y(includeIdx);
        mu = mean(y_included);
        sigma = std(y_included);

        % Identify inliers
        inlierIdx = abs(y_included - mu) <= thresholdSD * sigma;

        % Compute robust mean from inliers
        robustMean = mean(y_included(inlierIdx));

        % Subtract robust mean across entire time window
        eegData(ch, t1:t2, tr) = eegData(ch, t1:t2, tr) - robustMean;
    end
end

EEG.data = eegData;

fprintf('tesa_robustdemean: Demeaned data in [%g, %g] ms', timeWindow(1), timeWindow(2));
if ~isempty(excludeWindow)
    fprintf(' excluding [%g, %g] ms.\n', excludeWindow(1), excludeWindow(2));
else
    fprintf('.\n');
end

end
