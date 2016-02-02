% tesa_removedata() - removes data between defined timepoints and replaces
%                   data with NaNs. Stores removed timepoints.
%
% Usage:
%   >>  EEG = tesa_removedata( EEG, cutTimesTMS, cutTimesRec );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   cutTimesTMS     - vector with time range for removing TMS artifact [t1,t2]
%   cutTimesRec     - (optional) vector with time range for removing recharge artifact [t1,t2]
%    
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

function EEG = tesa_removedata( EEG, cutTimesTMS, cutTimesRec )

if nargin < 2
	error('Not enough input arguments.');
end

if nargin < 3
	cutTimesRec = [];
end

%Check that two time points have been specified
if size(cutTimesTMS,2) ~= 2
	error('Please provide two time values for TMS removal: [start cut, end cut]');
end

if ~isempty(cutTimesRec)
    if size(cutTimesRec,2) ~= 2
        error('Please provide two time values for recharge removal: [start cut, end cut]');
    end
end

%Check that time values given for TMS artifact are in range of data
if cutTimesTMS(1,1) < EEG.times(1,1) || cutTimesTMS(1,2) > EEG.times(1,end)
    error('Time values for TMS artifact removal are out of data range. Note that time values are in ms.');
end

%Check that time values given for recharge artifact are in range of data
if ~isempty(cutTimesRec)
    if cutTimesRec(1,1) < EEG.times(1,1) || cutTimesRec(1,2) > EEG.times(1,end)
        error('Time values for TMS artifact removal are out of data range. Note that time values are in ms.');
    end
end

%find the time values to remove TMS artifact
[val1,tp1] = min(abs(EEG.times-cutTimesTMS(1,1)));
[val2,tp2] = min(abs(EEG.times-cutTimesTMS(1,2)));

%find the time values to remove recharge artifact
if ~isempty(cutTimesRec)
    [val3,tp3] = min(abs(EEG.times-cutTimesRec(1,1)));
    [val4,tp4] = min(abs(EEG.times-cutTimesRec(1,2)));
 
elseif isempty(cutTimesRec)
    tp3 = [];
    tp4 = [];
end

%store the removed values
if ~isfield(EEG, 'tmscut')
    EEG.tmscut(1).cutTimesTMS = cutTimesTMS;
    EEG.tmscut(1).tmscut1 = tp1;
    EEG.tmscut(1).tmscut2 = tp2;
    EEG.tmscut(1).cutTimesRec = cutTimesRec;    
    EEG.tmscut(1).reccut1 = tp3;
    EEG.tmscut(1).reccut2 = tp4;
    EEG.tmscut(1).srate = EEG.srate;
    EEG.tmscut(1).interpolated = 'no';
elseif isfield(EEG, 'tmscut')
    num = size(EEG.tmscut,1)+1;
    EEG.tmscut(num).cutTimesTMS = cutTimesTMS;
    EEG.tmscut(num).tmscut1 = tp1;
    EEG.tmscut(num).tmscut2 = tp2;
    EEG.tmscut(num).cutTimesRec = cutTimesRec;     
    EEG.tmscut(num).reccut1 = tp3;
    EEG.tmscut(num).reccut2 = tp4;
    EEG.tmscut(num).srate = EEG.srate;
    EEG.tmscut(num).interpolated = 'no';
end

%replace TMS artifact data with 0s
EEG.data(:,tp1:tp2,:) = 0;

%display message
fprintf('TMS artifact data between %d ms and %d ms replaced with 0\n',cutTimesTMS(1,1),cutTimesTMS(1,2));

%replace recharge artifact data with 0s
if ~isempty(cutTimesRec)
    EEG.data(:,tp3:tp4,:) = 0;
    fprintf('Recharge artifact data between %d ms and %d ms replaced with 0\n',cutTimesRec(1,1),cutTimesRec(1,2));
end

end
