% tesa_removedata() - removes data between defined time points and replaces
%                   data with 0s or the average of a defined period (e.g. a baseline period).
%                   Can run on either continous or epoched data.
%                   Removed time points are stored in EEG.tmscut.
%
% Usage:
%   >>  EEG = tesa_removedata( EEG, cutTimesTMS ); %replace with 0s
%   >>  EEG = tesa_removedata( EEG, cutTimesTMS, replaceTimes , cutEvent );% replace with average of defined period in ms around provided event/s
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   cutTimesTMS     - (required) vector with time range for removing TMS artifact in ms. [t1,t2]
%                       Note that the times are relevant to events, not the epoch. For example, if an event is located at -100
%                       ms in an epoch, running the algorithm with [-10,10] would remove between -110 and -90 ms.
%                       Example: [-10,10]
%   replaceTimes    - (optional) vector with time range for calculating average to replace removed data in ms. [t1,t2]
%                       If not included, data will be replaced with 0s.
%                       As above, times are relative to events, not the epoch.
%                       Example: [-500, -100]
%   cutEvent        - (optional) cell with strings indicating event/s to remove data around. {'string'}. 
%                       If left empty, tesa_removedata will remove data for all events.
%                       Example: {'TMS'}
%                       Example: {'single','paired'}
%                       Example: {'1'}
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = tesa_removedata( EEG, [-10,10] ); %replace with 0s
%   EEG = tesa_removedata( EEG, [-10,10], [-500,-100], {'TMS'} ); %replace with average of defined period in ms around event 'TMS'
%   EEG = tesa_removedata( EEG, [-10,10], [], {'TMS'} ); %replace with 0s around event 'TMS'
% 
% See also:
%   tesa_interpdata 

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

% Change log:
% 15.6.2018: Changed function to work on both continuous and epoched data

function EEG = tesa_removedata( EEG, cutTimesTMS,  replaceTimes, cutEvent )

if nargin < 2
	error('Not enough input arguments. Note that you need to supply a time window to epoch around');
end

if nargin < 3
	replaceTimes = [];
end

if nargin < 4
	cutEvent = [];
end

%Check that two time points have been specified
if size(cutTimesTMS,2) ~= 2
	error('Please provide two time values for TMS removal in ms: [start cut, end cut]');
end

%Convert all events to strings
if ~isempty(EEG.event)
    for x = 1:size(EEG.event,2)
        if isnumeric(EEG.event(x).type)
            EEG.event(x).type = num2str(EEG.event(x).type);
        end
    end
end

% Convert numeric identifiers to strings
for conx = 1:length(EEG.event)
    if isnumeric(EEG.event(conx).type)
        EEG.event(conx).type = num2str(EEG.event(conx).type);
    end
end
for conx = 1:length(cutEvent)
    if isnumeric(cutEvent{conx})
        cutEvent{conx} = num2str(cutEvent{conx});
    end
end

% Event names
evNames = {EEG.event.type};

% Check if cutEvent is empty
if isempty(cutEvent)
    cutEvent = unique(evNames);
end

%Check that event exists
for evx = 1:length(cutEvent)
    if ~strcmp(cutEvent{evx},unique({EEG.event.type}))
        error('Event name ''%s'' does not exist. Please provide name of event present in the data.',cutEvent{evx});
    end
end
    
%Check that two time points have been specified for replaceTimes
if ~isempty(replaceTimes)    
    if size(replaceTimes,2) ~= 2
        error('Please provide two time values for the time period to average in ms: [start, end]');
    end
end

if size(EEG.data,3) > 1

    %Check that time values given for TMS artifact are in range of data
    if cutTimesTMS(1,1) < EEG.times(1,1) || cutTimesTMS(1,2) > EEG.times(1,end)
        error('Time values for TMS artifact removal are out of data range. Note that time values are in ms.');
    end

    %Check that time values given for replacement window are in range of data
    if ~isempty(replaceTimes)
        if replaceTimes(1,1) < EEG.times(1,1) || replaceTimes(1,2) > EEG.times(1,end)
            error('Time values for averaging data (replaceTimes) are out of data range. Note that time values are in ms.');
        end
    end

    %Check that time point 1 is smaller than time point 2 for TMS artifact
    if cutTimesTMS(1,1) > cutTimesTMS(1,2)
        error('t1 is larger than t2 for cutTimesTMS. t1 must be smaller than t2.');
    end

    %Check that time point 1 is smaller than time point 2 for replace times
    if ~isempty(replaceTimes)
        if replaceTimes(1,1) > replaceTimes(1,2)
            error('t1 is larger than t2 for replaceTimes. t1 must be smaller than t2.');
        end
    end
    
end

% Replacement
if ~isempty(replaceTimes)
    replacement = 'average';
else
    replacement = '0s';
end

% Check if data is continuous or epoched and replace data
if size(EEG.data,3) == 1 % Continuous data
    EEGtemp = EEG.data;
    dataType = 'continuous';
elseif size(EEG.data,3) > 1 % Epoched data
    EEGtemp = reshape(EEG.data,size(EEG.data,1),[],1);
    dataType = 'epoched';
end

% Find events
for evx = 1:length(cutEvent)
    evLogTemp(:,evx) = double(strcmp(cutEvent{evx},evNames));
end

% Convet to logical
evLog = sum(evLogTemp,2) > 0;

% Extract latencies
evLat = [EEG.event.latency];
evLatOut = round(evLat(evLog));

% Determine time periods for removal
for x = 1:length(evLatOut)
    cutWin(x,1) = evLatOut(1,x) + round(cutTimesTMS(1,1)*(EEG.srate/1000));
    cutWin(x,2) = evLatOut(1,x) + round(cutTimesTMS(1,2)*(EEG.srate/1000));
    
    % Check if cut window is outside data
    if cutWin(x,1)<0 || cutWin(x,2)<0 
        error('Window for event type ''%s'' is outside data limits. Either reduce the cut window or select a different event.',EEG.event(evLatOut(1,x)).type);
    elseif cutWin(x,1)>size(EEGtemp,2) || cutWin(x,2)>size(EEGtemp,2)
        error('Window for event type ''%s'' is outside data limits. Either reduce the cut window or select a different event.',EEG.event(evLatOut(1,x)).type);
    end
    
end

if isempty(replaceTimes)
    
    %replace TMS artifact data with 0s
    for x = 1:length(evLatOut)
        EEGtemp(:,cutWin(x,1):cutWin(x,2),:) = 0;
    end
    
    if strcmp('continuous',dataType)
        EEG.data = EEGtemp;
    elseif strcmp('epoched',dataType)
        EEG.data = reshape(EEGtemp,size(EEG.data,1),size(EEG.data,2),[]);
    end
    
    %display message
    fprintf('TMS artifact data between %d ms and %d ms replaced with 0\n',cutTimesTMS(1,1),cutTimesTMS(1,2));
    
else
    
    % Determine time periods for baseline data
    for x = 1:length(evLatOut)
        avWin(x,1) = evLatOut(1,x) + round(replaceTimes(1,1)*(EEG.srate/1000));
        avWin(x,2) = evLatOut(1,x) + round(replaceTimes(1,2)*(EEG.srate/1000));
    end
    
    % Replace data with baseline average
    for x = 1:length(evLatOut)
        
        % Calcuate baseline average for each event
        avPer = mean(EEGtemp(:,avWin(x,1):avWin(x,2),:),2);
        
        % Create an empty matrix
        newMat = zeros(size(EEGtemp,1),size(EEGtemp(:,cutWin(x,1):cutWin(x,2),:),2),size(EEGtemp,3));
        
        % Fill empty matrix with average data
        for a = 1:size(newMat,2)
            newMat(:,a,:) = avPer;
        end
        
        % Replace data with baseline average
        EEGtemp(:,cutWin(x,1):cutWin(x,2),:) = newMat;
        
    end
    
    if strcmp('continuous',dataType)
        EEG.data = EEGtemp;
    elseif strcmp('epoched',dataType)
        EEG.data = reshape(EEGtemp,size(EEG.data,1),size(EEG.data,2),[]);
    end
    
    fprintf('TMS artifact data between %d ms and %d ms replaced with average data\n',cutTimesTMS(1,1),cutTimesTMS(1,2));
end

%store the removed values
if ~isfield(EEG, 'tmscut')
    EEG.tmscut(1).cutTimesTMS = cutTimesTMS;
    EEG.tmscut(1).cutEvent = cutEvent;
    EEG.tmscut(1).replaceTimes = replaceTimes;
    EEG.tmscut(1).replacement = replacement;
    EEG.tmscut(1).srate = EEG.srate;
    EEG.tmscut(1).interpolated = 'no';
elseif isfield(EEG, 'tmscut')
    num = size(EEG.tmscut,2)+1;
    EEG.tmscut(num).cutTimesTMS = cutTimesTMS;
    EEG.tmscut(num).cutEvent = cutEvent;
    EEG.tmscut(num).replaceTimes = replaceTimes;
    EEG.tmscut(num).replacement = replacement;
    EEG.tmscut(num).srate = EEG.srate;
    EEG.tmscut(num).interpolated = 'no';
end

end
