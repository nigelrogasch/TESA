% pop_tesa_fixevent() - finds TMS pulses by detecting the large TMS artifacts
%                   present in already epoched data. This script is
%                   designed for instances when the recorded events do
%                   not correspond with when the TMS pulse was given.
%                   The script works by extracting a
%                   single channel and finding the time points in which the 
%                   first derivatives exceed a certain threshold (defined 
%                   by 'rate'). Paired pulses can also be detected.
% 
%                   IMPORTANT: If you need to use this function, make sure
%                   that the initial epochs you use are larger than the
%                   final epoch size you desire. If the initial epoch size is
%                   too small, the new epoch window will be out of range 
%                   with the new events. E.g. inital epoch is -100 to 100
%                   and the event is shifted 10 ms so the new 0 now sits
%                   at 10 ms. Re-epoching the data to -100 to 100 won't
%                   work as the new range is effectively -90 to 110. In this case, 
%                   run the initial epoch at -120 to 120 and the epoch can now
%                   be taken.
%                   
%
% Usage:
%   >>  EEG = pop_tesa_fixevent( EEG ); % pop up window
%   >>  EEG = pop_tesa_fixevent( EEG, elec, newEpoch, tmsLabel );
%   >>  EEG = pop_tesa_fixevent( EEG, elec, newEpoch, tmsLabel, 'key1', value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   elec            - [required] string with electrode to use for finding artifact
%   newEpoch        - [required] vector with start and end time of new epoch in
%                   seconds (following pop_epoch convention). 
%                   Example: [-1,1] %For -1 s to 1s epoch
%   tmsLabel        - [required] string indicating the event that requires 
%                   correcting (e.g. 'TMS')  
%                   
% Optional input pairs:
%   'refract',int   - int defines the refractory period (the time for the
%                   TMS artefact to recover below the rate of change). int
%                   is in ms. 
%                   default = 3
%   'rate',int      - int defines the rate of change of the TMS artifact in
%                   uV/ms. 
%                   default = 1e4
%  
% Input pairs for detecting paired pulses
%   'paired','str'  - [required] 'str' - type 'yes' to turn on paired detection
%                   default = 'no'
%   'ISI', [int]    - [required] [int] is a vector defining interstimulus intervals
%                   between conditioning and test pulses. Multiple ISIs can 
%                   be defined as [1,2,...]. 
%                   default = []
%     
% Outputs:
%   EEG             - EEGLAB EEG structure
%
%   % Examples
%   EEG = pop_tesa_fixevent( EEG, 'Cz', [-0.8,0.8], 'TMS' ); %default use
%   EEG = pop_tesa_fixevent( EEG, 'Fz', [-0.7,0.7], 'TMS', 'refract', 4, 'rate', 2e5 ); %user defined input
%   EEG = pop_tesa_fixevent( EEG, 'Cz', [-0.8,0.8], 'LICI', 'paired', 'yes', 'ISI', 100 ); %paired pulse use
%
% See also:
%   pop_tesa_findpulse, pop_tesa_findpulsepeak

% Copyright (C) 2016  Nigel Rogasch, Monash University,
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

 function [EEG com] = pop_tesa_fixevent( EEG, elec, newEpoch, tmsLabel, varargin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
    
    for x =1:size(EEG.chanlocs,2)
        chanAll{x,1} = EEG.chanlocs(x).labels;
    end
    
    % Find default channel (CZ)
    if sum(strcmpi('CZ',chanAll)) > 0
%         chan = chanAll{strcmpi('CZ',chanAll)};
        chanVal = find(strcmpi('CZ',chanAll));
    else
%         chan = chanAll{1,1};
        chanVal = 1;
    end
    
    % Create channel list
    chanOr = chanAll;
    for x = 1:size(chanAll,1)-1
        chanOr{x,1} = [chanOr{x,1}, '|'];
    end
    chanIn = [chanOr{:}];
    
    % Create event list
    for x = 1:size(EEG.event,2)
        if ~ischar(EEG.event(x).type)
            EEG.event(x).type = num2str(EEG.event(x).type);
        end
        trigAll{x,1} = EEG.event(x).type;
    end
    trigUnique = unique(trigAll);
    trigOr = trigUnique;
    for x = 1:size(trigOr,1)-1
        trigOr{x,1} = [trigOr{x,1}, '|'];
    end
    trigIn = [trigOr{:}];
    
    % Find default event ('TMS')
    if sum(strcmpi('TMS',trigUnique)) > 0
        trigVal = find(strcmpi('TMS',trigUnique));
    else
        trigVal = 1;
    end
    
    geometry = {1 [1 0.3] [1 0.3] [1 0.3] [1 0.3] [1 0.3] 1 [1 0.3] [1 0.3]};

    uilist = {{'style', 'text', 'string', 'Fix event position','fontweight','bold'} ...
              {'style', 'text', 'string', 'New epoch length (secs) [required]'} ...
              {'style', 'edit', 'string', '-1, 1'} ...
              {'style', 'text', 'string', 'Electrode for finding artifact'} ...
              {'style', 'popupmenu', 'string', chanIn, 'tag', 'interp', 'Value', chanVal } ...
              {'style', 'text', 'string', 'Refractory period - time for TMS artifact to recover (ms)'} ...
              {'style', 'edit', 'string', '3'}...
              {'style', 'text', 'string', 'Rate of change of TMS artifact - used to find artifact (uV/ms)'} ...
              {'style', 'edit', 'string', '1e4'}...
              {'style', 'text', 'string', 'Event to be adjusted.'} ...
              {'style', 'popupmenu', 'string', trigIn, 'tag', 'interp', 'Value', trigVal } ...
              {} ...
              {'style', 'text', 'string', 'Paired pulse TMS','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Interstimulus interval (ms) [required]'} ...
              {'style', 'edit', 'string', ''}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Fix event position -- pop_tesa_fixevent()', 'helpcom', 'pophelp(''pop_tesa_fixevent'')');
    if isempty(result), return; end;
    
    %Extract data for single pulse artifact find
    newEpoch = str2num(result{1,1});
    elec = chanAll{result{1,2},1};
    refract = str2num(result{1,3});
    rate = str2num(result{1,4});
    tmsLabel = trigUnique{result{1,5},1};
    
    %Check if correct information is provided
    if isempty(refract)
        refract = 3;      
    end 
    if isempty(rate)
        rate = 1e4;      
    end
    if isempty(tmsLabel)
        tmsLabel = 'TMS';      
    end
    
    %Check for paired option
    if result{1,6} == 1 %paired on
        paired = 'yes';
        ISI = str2num(result{1,7});
    end
    
end

%Run script from input
if nargin == 4;
    EEG = tesa_fixevent(EEG,elec,newEpoch,tmsLabel);
    com = sprintf('%s = pop_tesa_fixevent( %s, ''%s'', %s, ''%s'' );', inputname(1), inputname(1), elec, mat2str(newEpoch), tmsLabel );
elseif nargin > 4
    EEG = tesa_fixevent(EEG,elec,newEpoch,tmsLabel,varargin{:});
    com = sprintf('%s = pop_tesa_fixevent( %s,''%s'',%s,''%s'',%s );', inputname(1), inputname(1), elec, mat2str(newEpoch), tmsLabel, vararg2str(varargin) );
end
    
%find artifact and return the string command using pop window info
if nargin < 2
    if result{1,6}==0
        EEG = tesa_fixevent( EEG, elec, newEpoch, tmsLabel, 'refract', refract, 'rate', rate);
        com = sprintf('%s = pop_tesa_fixevent( %s, ''%s'', %s, ''%s'', ''refract'', %s, ''rate'', %s );', inputname(1), inputname(1), elec, mat2str(newEpoch), tmsLabel, mat2str(refract), mat2str(rate));
    elseif result{1,6}==1
        EEG = tesa_fixevent( EEG, elec, newEpoch, tmsLabel, 'refract', refract, 'rate', rate, 'paired', paired, 'ISI', ISI );
        com = sprintf('%s = pop_tesa_fixevent( %s, ''%s'', %s, ''%s'', ''refract'', %s, ''rate'', %s, ''paired'', ''%s'', ''ISI'', %s);', inputname(1), inputname(1), elec, mat2str(newEpoch), tmsLabel, mat2str(refract), mat2str(rate), paired, mat2str(ISI));
    end
end

end
