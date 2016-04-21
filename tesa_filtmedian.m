% tesa_filtmedian() - applies a 1-d median filter of nth-order to remove
%                       artifacts such as spikes and muscle artifacts. 
%                       Note that spike artifacts can be selected using the 
%                       gui option in tesa_findpulsepeak.
%
% Usage:
%   >>  EEG = tesa_filtmedian( EEG, timeWin, filtOrd, label );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   timeWin         - Vector with time range for applying median filter in ms [t1,t2] (example: [-1,20])
%                      Note that t1 must be 0 or negative and t2 positive
%   filtOrd         - Integer indicating filter order (example: 30). The
%                      filter order determines the number of samples considered 
%                      either side of the data point when calculating the median value. 
%   label           - String indicating which event to filter around (example: 'TMS' or 'spike')
% 
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples:
%   EEG = tesa_filtmedian( EEG, [-2,2], 3, 'spike' ); %Apply a 3rd order median filter to remove small spike artifacts
%
% See also:
%   tesa_filtbutter

% Copyright (C) 2016  Nigel Rogasch, Monash University,
% nigel.rogasch@monash.edu
%
% Authors: Nathan Rose, Nigel Rogasch
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

function EEG = tesa_filtmedian( EEG, timeWin, filtOrd, label )

if nargin < 4
	error('Not enough input arguments.');
end

%Check that two time points have been specified
if size(timeWin,2) ~= 2
	error('Please provide two time values for median filter: [start time (t1), end time (t2)]');
end

%Check that time values given for timeWin are correct
if timeWin(1,1) > 0 
    error('Time values for median filter are incorrect. t1 must be 0 or negative.');
elseif timeWin(1,2) < 0
    error('Time values for median filter are incorrect. t2 must be positive.');
end

data = reshape(EEG.data,size(EEG.data,1),[]);

%Check that event types are strings
for a = 1:size(EEG.event,2)
    if ~ischar(EEG.event(a).type)
        EEG.event(a).type = num2str(EEG.event(a).type);
    end
end

%Check that labels exist
events = [{EEG.event.type}];
labelCheck = unique(events);
labelLog = strcmp(label,labelCheck);
if sum(labelLog) == 0
    error('The event %s is not present in the data.', label)
end

%Find events for median filtering
latency = [EEG.event.latency];
trigind = strcmp(label,events);
lat = latency(trigind);

%Convert time window to samples
tp1 = round(timeWin(1,1)/(1/EEG.srate*1000));
tp2 = round(timeWin(1,2)/(1/EEG.srate*1000));

%Apply median filter to data
for i = 1:length(lat)
    win = lat(i)+tp1:lat(i)+tp2;
%     filt = medfilt1(double(data(:,win)),filtOrd);
    filt = tesa_medfilt(data,win,filtOrd);
    data(:,win) = filt;
end

EEG.data = reshape(data,size(EEG.data,1),size(EEG.data,2),[]);

%display message
fprintf('Median filter applied between %d ms and %d ms\n',timeWin(1,1),timeWin(1,2));

end
