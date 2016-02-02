% tesa_medianfilt() - applies a 1-d median filter of nth-order to remove
%                       muscle and decay artefacts.
%
% Usage:
%   >>  EEG = tesa_medianfilt( EEG, timeWin );
%   >>  EEG = tesa_medianfilt( EEG, timeWin, filtOrd );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   timeWin         - vector with time range for applying median filter [t1,t2]
%   filtOrd         - (optional) integer indicating filter order 
%                       default = [29]
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% See also:
%   SAMPLE, EEGLAB 

% Copyright (C) 2015  Nigel Rogasch, Monash University,
% nigel.rogasch@monash.edu
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

function EEG = tesa_medianfilt( EEG, timeWin, filtOrd )

if nargin < 2
	error('Not enough input arguments.');
end

if nargin < 3
	filtOrd = 29;
end

%Check if filtOrd has been specified
if isempty(filtOrd)
    filtOrd = 29;
end

%Check that two time points have been specified
if size(timeWin,2) ~= 2
	error('Please provide two time values for median filter: [start time, end time]');
end

%Check that time values given for TMS artifact are in range of data
if timeWin(1,1) < EEG.times(1,1) || timeWin(1,2) > EEG.times(1,end)
    error('Time values for median filter are out of data range. Note that time values are in ms.');
end

%find the time values to remove TMS artifact
[val1,tp1] = min(abs(EEG.times-timeWin(1,1)));
[val2,tp2] = min(abs(EEG.times-timeWin(1,2)));

%Apply median filter to TMS data
for a = 1:size(EEG.data,3)
    EEG.data(:,tp1:tp2,a) = medfilt1(double(EEG.data(:,tp1:tp2,a)),filtOrd); 
end

%display message
fprintf('Median filter applied between %d ms and %d ms\n',timeWin(1,1),timeWin(1,2));

end
