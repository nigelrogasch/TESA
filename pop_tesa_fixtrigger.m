% pop_tesa_fixtrigger() - finds TMS pulses by detecting the large TMS artifacts
%                   present in already epoched data. This script is
%                   designed for instances when the recoreded triggers do
%                   not correspond with when the TMS pulse was given.
%                   The script works by extracting a
%                   single channel and finding the time points in which the 
%                   first derivatives exceed a certain threshold (defined 
%                   by 'rate'). Paired pulses and repetitive TMS trains can
%                   also be deteceted.
% 
%                   IMPORTANT: If you need to use this function, make sure
%                   that the initial epochs you use are larger than the
%                   final epoch size you desire. If the initial epoch size is
%                   too small, the new epoch window will be out of range of
%                   with the new triggers. E.g. inital epoch is -100 to 100
%                   and the trigger is shifted 10 ms so the new 0 now sits
%                   at 10 ms. Re-epoching the data to -100 to 100 won't
%                   work as the new range is effectively -90 to 110. In this case, 
%                   run the initial epoch at -120 to 120 and the epoch can now
%                   be taken.
%                   
%
% Usage:
%   >>  EEG = pop_tesa_fixtrigger( EEG );
%   >>  EEG = pop_tesa_fixtrigger( EEG, elec, newEpoch );
%   >>  EEG = pop_tesa_fixtrigger( EEG, elec, newEpoch, varargin );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   elec            - string with electrode to use for finding artifact
%   newEpoch        - vector with start and end time of new epoch in
%                   seconds (following pop_epoch convention). 
%                   Example: [-1,1] %For -1 s to 1s epoch
% 
% Optional input pairs:
%   'refract',int   - int defines the refractory period (the time for the
%                   TMS artefact to recover below the rate of change). int
%                   is in ms. 
%                   default = 3
%   'rate',int      - int defines the rate of change of the TMS artifact in
%                   uV/ms. 
%                   default = 1e5
%   'tmsLabel','str'- 'str' is a string for the single TMS label.  
%                   default = 'TMS'
%  
% Input pairs for detecting paired pulses
%   'paired','str'  - required. 'str' - type 'yes' to turn on paired detection
%                   default = 'no'
%   'ISI', [int]    - required. [int] is a vector defining interstimulus intervals
%                   between conditioning and test pulses. Multiple ISIs can 
%                   be defined as [1,2,...]. 
%                   default = []
%   'pairLabel',{'str'} - required if more than 1 ISI. {'str'} is a cell array
%                   containing string labels for different ISI conditions.  
%                   Multiple labels can be defined as {'SICI','LICI',...}.
%                   The number of labels defined must equal the number of
%                   ISI conditions defined.
%                   default = {'TMSpair'}
% 
%  Input pairs for detecting repetitive TMS trains
%  'repetitive','str' - required. 'str' - type 'yes' to turn on repetitive detection
%                   default = 'no'
%   'ITI', int      - required. int defines the inter-train interval in ms.
%                   For example, if a 10 Hz rTMS condition is used with 4s
%                   of stimulation (40 pulses) and 26s of rest, ITI = 2600;
%                   default = []
%   'pulseNum', int - required. int defines the number of pulses in a
%                   train. Using the above example, this would be 40. 
%                   deafult = []
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

 function [EEG com] = pop_tesa_findpulse( EEG, elec, varargin )

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
    
    if sum(strcmpi('CZ',chanAll)) > 0
        chan = chanAll{strcmpi('CZ',chanAll)};
    else
        chan = chanAll{1,1};
    end
    
    geometry = {1 [1 0.3] [1 0.3] [1 0.3] [1 0.3] [1 0.3] 1 [1 0.3] [1 0.3] 1 [1 0.3] 1 1 1 [1 0.3] [1 0.3] [1 0.3]};

    uilist = {{'style', 'text', 'string', 'Fix trigger position','fontweight','bold'} ...
              {'style', 'text', 'string', 'New epoch length (secs) [required]'} ...
              {'style', 'edit', 'string', '-1, 1'} ...
              {'style', 'text', 'string', 'Electrode for finding artifact'} ...
              {'style', 'edit', 'string', chan} ...
              {'style', 'text', 'string', 'Refractory period - time for TMS artifact to recover (ms)'} ...
              {'style', 'edit', 'string', '3'}...
              {'style', 'text', 'string', 'Rate of change of TMS artifact - used to find artifact (uV/ms)'} ...
              {'style', 'edit', 'string', '1e5'}...
              {'style', 'text', 'string', 'Label for TMS pulse (single).'} ...
              {'style', 'edit', 'string', 'TMS'}...
              {} ...
              {'style', 'text', 'string', 'Paired pulse TMS','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'Yes?' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Interstimulus interval (ms) [required]'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', '     Multiple ISIs can be entered as follows [2,15,100]','fontangle','italic'} ...
              {'style', 'text', 'string', 'Label for paired pulses (e.g. SICI) [required]'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', '     Multiple labels entered as follows (e.g. SICI, ICF, LICI)','fontangle','italic'} ...
              {'style', 'text', 'string', '     Number of labels must equal number of ISIs','fontangle','italic'}...
              {} ...
              {'style', 'text', 'string', 'Repetitive TMS','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'Yes?' 'value' 0 'tag' 'rep' } ...
              {'style', 'text', 'string', 'Inter-train interval  (ms) [required]'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'Number of pulses in train [required]'} ...
              {'style', 'edit', 'string', ''}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Fix trigger position -- pop_tesa_fixtrigger()', 'helpcom', 'pophelp(''tesa_fixtrigger'')');
    
    %Check that both paired and repetitive are not on
    if result{1,6}==1 && result{1,9} == 1
        error('tesa_findpulse can not search for both paired and repetitive stimuli within the same file. Please choose one.');
        end6
    
    %Extract data for single pulse artifact find
    newEpoch = str2mat(result{1,1});
    elec = result{1,2};
    refract = str2num(result{1,3});
    rate = str2num(result{1,4});
    tmsLabel = result{1,5};
    
    %Check if correct information is provided
    if isempty(elec)
        error('Electrode name not entered - this is required to find artifact. Script terminated')       
    end
    if isempty(refract)
        refract = 3;      
    end 
    if isempty(rate)
        rate = 1e5;      
    end
    if isempty(tmsLabel)
        tmsLabel = 'TMS';      
    end
    
    %Check for paired option
    if result{1,6} == 1 %paired on
        paired = 'yes';
        ISI = str2num(result{1,7});
        pairLabel = strtrim(strsplit(result{1,8},','));
    end
    
    %Check for repetitive option
    if result{1,9} == 1 %repetitive on
        repetitive = 'yes';
        ITI = str2num(result{1,10});
        pulseNum = str2num(result{1,11});
    end

end

%Run script from input
if nargin == 3;
    EEG = tesa_fixtrigger(EEG,elec,newEpoch);
    com = sprintf('%s = pop_tesa_fixtrigger( %s, %s, %s );', inputname(1), inputname(1), elec, newEpoch );
elseif nargin > 3
    EEG = tesa_fixtrigger(EEG,elec,newEpoch,varargin{:});
    com = sprintf('%s = pop_tesa_fixtrigger( %s, %s, %s, %s );', inputname(1), inputname(1), elec, newEpoch, vararg2str(varargin) );
end
    

%find artifact and return the string command using pop window info
if nargin < 2
    if result{1,6}==0 && result{1,9} == 0
        EEG = tesa_ffixtrigger( EEG, elec, newEpoch, 'refract', refract, 'rate', rate, 'tmsLabel', tmsLabel );
        com = sprintf('%s = pop_tesa_fixtrigger( %s, %s, %s, ''refract'', %s, ''rate'', %s, ''tmsLabel'', %s);', inputname(1), inputname(1), elec, newEpoch, mat2str(refract), mat2str(rate), tmsLabel);
    elseif result{1,6}==1
        if strcmp(pairLabel{1,1},'')
            EEG = tesa_fixtrigger( EEG, elec, newEpoch, 'refract', refract, 'rate', rate, 'tmsLabel', tmsLabel, 'paired', paired, 'ISI', ISI );
            com = sprintf('%s = pop_tesa_fixtrigger( %s, %s, %s, ''refract'', %s, ''rate'', %s, ''tmsLabel'', %s, ''paired'', %s, ''ISI'', %s);', inputname(1), inputname(1), elec, newEpoch, mat2str(refract), mat2str(rate), tmsLabel, paired, mat2str(ISI));
        else
            EEG = tesa_fixtrigger( EEG, elec, newEpoch, 'refract', refract, 'rate', rate, 'tmsLabel', tmsLabel, 'paired', paired, 'ISI', ISI, 'pairLabel', pairLabel );
            com = sprintf('%s = pop_tesa_fixtrigger( %s, %s, %s, ''refract'', %s, ''rate'', %s, ''tmsLabel'', %s, ''paired'', %s, ''ISI'', %s, ''pairLabel'', {%s});', inputname(1), inputname(1), elec, newEpoch, mat2str(refract), mat2str(rate), tmsLabel, paired, mat2str(ISI), result{1,8});
        end
    elseif result{1,9}==1
        EEG = tesa_fixtrigger( EEG, elec, newEpoch, 'refract', refract, 'rate', rate, 'tmsLabel', tmsLabel, 'repetitive', repetitive, 'ITI', ITI, 'pulseNum', pulseNum );
        com = sprintf('%s = pop_tesa_fixtrigger( %s, %s, %s, ''refract'', %s, ''rate'', %s, ''tmsLabel'', %s, ''repetitive'', %s, ''ITI'', %s, ''pulseNum'', %s);', inputname(1), inputname(1), elec, newEpoch, mat2str(refract), mat2str(rate), tmsLabel, repetitive, mat2str(ITI), mat2str(pulseNum));
    end
end

end
