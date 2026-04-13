% pop_tesa_robustdetrend() - performs robust polynomial detrending on EEG data 
%                       within a specified time window, optionally excluding a sub-window
%                       (e.g., to avoid fitting over TMS artefacts or tep components). 
%                       It removes slow trends from the signal while reducing the influence 
%                       of outliers.
%
% Usage:
%   >>  EEG = pop_tesa_robustdetrend(EEG); % GUI pop-up
%   >>  EEG = pop_tesa_robustdetrend(EEG, timeWindow, thresholdSD, polyOrder, excludeWindow);
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
%   EEG = pop_tesa_robustdetrend(EEG, [-1000,1000], 3, 1); % performs a linear detrend between -1000 to 1000 ms not including outlier data > 3 standard deviations from the mean.
%   EEG = pop_tesa_robustdetrend(EEG, [-1000,0], 5, 2); % performs a quadratic detrend between -1000 to 0 ms not including outlier data > 5 standard deviations from the mean.
%   EEG = pop_tesa_robustdetrend(EEG, [-1000,1000], 3, 1, [-10,500]); % performs a linear detrend between -1000 to 1000 ms not including outlier data > 3 standard deviations from the mean and excluding data between -10 and 500 ms.
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

function [EEG, com] = pop_tesa_robustdetrend(EEG, timeWindow, thresholdSD, polyOrder, excludeWindow)

% Check inputs or launch GUI if missing
if nargin < 2

    % -------------------------------
    % GUI mode
    % -------------------------------
    
    % Find the xmin and xmax
    xmin = min(EEG.times);
    xmax = max(EEG.times);
    trange = sprintf('%g, %g',xmin, xmax);

    geometry = {1 [1 0.3] [1 0.3] [1 0.3] [1 0.3]};

    uilist = {{'style', 'text', 'string', 'Robust detrend','fontweight','bold'} ...
              {'style', 'text', 'string', 'Time to detrend (ms): start, end'} ...
              {'style', 'edit', 'string', trange} ...
              {'style', 'text', 'string', 'Outlier rejection threshold (standard deviations)'} ...
              {'style', 'edit', 'string', '3'}...
              {'style', 'text', 'string', 'Polynomial order (e.g., 1 = linear, 2 = quadratic)'} ...
              {'style', 'edit', 'string', '1'}...
              {'style', 'text', 'string', 'Time to exclude (ms): start, end (Optional)'} ...
              {'style', 'edit', 'string', ''}...
              };
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Modified band-pass filter with autoregressive padding -- pop_tesa_modifiedbandpassfilter()', 'helpcom', 'pophelp(''pop_tesa_modifiedbandpassfilter'')');
    if isempty(result), return; end;

    % Parse GUI inputs
    try
        timeWindow   = str2num(result{1}); %#ok<ST2NM>
        thresholdSD  = str2double(result{2});
        polyOrder    = str2double(result{3});
        excludeWindow = str2num(result{4}); %#ok<ST2NM>
        if isempty(excludeWindow)
            excludeWindow = [];
        end
    catch
        error('Invalid input format. Please check your entries.');
    end
end

if exist('excludeWindow', 'var') ~= 1
    excludeWindow = [];
end

% -------------------------------
% Run detrending function
% -------------------------------
fprintf('Running robust polynomial detrending...\n');
EEG = tesa_robustdetrend(EEG, timeWindow, thresholdSD, polyOrder, excludeWindow);

% -------------------------------
% Update EEGLAB history
% -------------------------------
if isempty(excludeWindow)
    com = sprintf('EEG = pop_tesa_robustdetrend(EEG, [%g %g], %g, %g);', ...
        timeWindow(1), timeWindow(2), thresholdSD, polyOrder);
else
    com = sprintf('EEG = pop_tesa_robustdetrend(EEG, [%g %g], %g, %g, [%g %g]);', ...
        timeWindow(1), timeWindow(2), thresholdSD, polyOrder, excludeWindow(1), excludeWindow(2));
end

fprintf('Detrending complete.\n');

end