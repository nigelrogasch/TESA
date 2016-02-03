% tesa_removedata() - removes data between defined time points and replaces
%                   data with 0s or the average of a defined period (e.g. a baseline period). 
%                   Removed time points are stored in EEG.tmscut.
%
% Usage:
%   >>  EEG = tesa_removedata( EEG, cutTimesTMS ); %replace with 0s
%   >>  EEG = tesa_removedata( EEG, cutTimesTMS, replaceTimes );% replace with average of defined period in ms
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   cutTimesTMS     - (required) vector with time range for removing TMS artifact in ms. [t1,t2]
%                       Example: [-10,10]
%   replaceTimes    - (optional) vector with time range for calculating average to replace removed data in ms. [t1,t2]
%                      If not included, data will be replaced with 0s.
%                       Example: [-500, -100]
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = tesa_removedata( EEG, [-10,10] ); %replace with 0s
%   EEG = tesa_removedata( EEG, [-10,10], [-500,-100] ); %replace with average of defined period in ms
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

function EEG = tesa_removedata( EEG, cutTimesTMS, replaceTimes )

if nargin < 2
	error('Not enough input arguments.');
end

if nargin < 3
	replaceTimes = [];
end

%Check that two time points have been specified
if size(cutTimesTMS,2) ~= 2
	error('Please provide two time values for TMS removal in ms: [start cut, end cut]');
end

%Check that two time points have been specified for replaceTimes
if ~isempty(replaceTimes)    
    if size(replaceTimes,2) ~= 2
        error('Please provide two time values for the time period to average in ms: [start, end]');
    end
end

%Check that time values given for TMS artifact are in range of data
if cutTimesTMS(1,1) < EEG.times(1,1) || cutTimesTMS(1,2) > EEG.times(1,end)
    error('Time values for TMS artifact removal are out of data range. Note that time values are in ms.');
end

%Check that time values given for TMS artifact are in range of data
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

%find the time values to remove TMS artifact
[val1,tp1] = min(abs(EEG.times-cutTimesTMS(1,1)));
[val2,tp2] = min(abs(EEG.times-cutTimesTMS(1,2)));

%find the time values to calculate average to replace data
if ~isempty(replaceTimes)
    [val3,tp3] = min(abs(EEG.times-replaceTimes(1,1)));
    [val4,tp4] = min(abs(EEG.times-replaceTimes(1,2)));
else
    tp3 = [];
    tp4 = [];
end

% Replacement
if ~isempty(replaceTimes)
    replacement = 'average';
else
    replacement = '0s';
end

%store the removed values
if ~isfield(EEG, 'tmscut')
    EEG.tmscut(1).cutTimesTMS = cutTimesTMS;
    EEG.tmscut(1).replaceTimes = replaceTimes;
    EEG.tmscut(1).replacement = replacement;
    EEG.tmscut(1).tmscut1 = tp1;
    EEG.tmscut(1).tmscut2 = tp2;
    EEG.tmscut(1).tmscut3 = tp3;
    EEG.tmscut(1).tmscut4 = tp4;
    EEG.tmscut(1).srate = EEG.srate;
    EEG.tmscut(1).interpolated = 'no';
elseif isfield(EEG, 'tmscut')
    num = size(EEG.tmscut,2)+1;
    EEG.tmscut(num).cutTimesTMS = cutTimesTMS;
    EEG.tmscut(num).replaceTimes = replaceTimes;
    EEG.tmscut(num).replacement = replacement;
    EEG.tmscut(num).tmscut1 = tp1;
    EEG.tmscut(num).tmscut2 = tp2;
    EEG.tmscut(num).tmscut3 = tp3;
    EEG.tmscut(num).tmscut4 = tp4;
    EEG.tmscut(num).srate = EEG.srate;
    EEG.tmscut(num).interpolated = 'no';
end

if isempty(replaceTimes)
    %replace TMS artifact data with 0s
    EEG.data(:,tp1:tp2,:) = 0;

    %display message
    fprintf('TMS artifact data between %d ms and %d ms replaced with 0\n',cutTimesTMS(1,1),cutTimesTMS(1,2));
else
    %replace TMS artifact data with average of defined period
    avPer = mean(EEG.data(:,tp3:tp4,:),2);
    
    newMat = zeros(size(EEG.data,1),size(EEG.data(:,tp1:tp2,:),2),size(EEG.data,3));
    
    for a = 1:size(newMat,2);
        newMat(:,a,:) = avPer;
    end
    
    EEG.data(:,tp1:tp2,:) = newMat;
       
    %display message
    fprintf('TMS artifact data between %d ms and %d ms replaced with average data\n',cutTimesTMS(1,1),cutTimesTMS(1,2));
    
end

end
