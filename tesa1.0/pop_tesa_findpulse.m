% pop_tesa_findpulse() - finds TMS pulses by detecting the large TMS artifacts
%                   present in the data. This script works by extracting a
%                   single channel and finding the time points in which the 
%                   first derivatives exceed a certain threshold (defined 
%                   by 'rate'). Paired pulses and repetitive TMS trains can
%                   also be deteceted. 
%
% Usage:
%   >>  EEG = pop_tesa_findpulse( EEG ); % pop up window
%   >>  EEG = pop_tesa_findpulse( EEG, elec );
%   >>  EEG = pop_tesa_findpulse( EEG, elec, 'key1', value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   elec            - string defining electrode to use for finding artifact
%                   example: 'Cz'
% 
% Optional input pairs:
%   'refract',int   - int defines the refractory period (the time for the
%                   TMS artifact to recover below the rate of change). int
%                   is an integer in ms. 
%                   default = 3
% 
%   'rate',int      - int defines the rate of change of the TMS artifact in
%                   uV/ms. 
%                   default = 1e4
% 
%   'tmsLabel','str'- 'str' is a string for the single TMS label.  
%                   default = 'TMS'
%
%   'plots','str' - 'on'|'off'. Brings up a plot showing the detected
%                   peaks. Red = detected 
%                   default = 'on'
%  
% Input pairs for detecting paired pulses
%   'paired','str'  - required. 'str' - type 'yes' to turn on paired detection
%                   default = 'no'
% 
%   'ISI', [int]    - required. [int] is a vector defining interstimulus intervals
%                   between conditioning and test pulses. Multiple ISIs can 
%                   be defined as [3,100,...]. 
%                   default = []
% 
%   'pairLabel',{'str'} - required if more than 1 ISI. {'str'} is a cell array
%                   containing string labels for different ISI conditions.  
%                   Multiple labels can be defined as {'SICI','LICI',...}.
%                   The number of labels defined must equal the number of
%                   ISI conditions defined.
%                   default = {'TMSpair'}
% 
% Input pairs for detecting repetitive TMS trains
%  'repetitive','str' - required. 'str' - type 'yes' to turn on repetitive detection
%                   default = 'no'
% 
%   'ITI', int      - required. int defines the inter-train interval in ms.
%                   For example, if a 10 Hz rTMS condition is used with 4s
%                   of stimulation (40 pulses) and 26s of rest, ITI = 2600;
%                   default = []
% 
%   'pulseNum', int - required. int defines the number of pulses in a
%                   train. Using the above example, this would be 40. 
%                   deafult = []
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples
%   EEG = pop_tesa_findpulse( EEG, 'Cz' ); %default use
%   EEG = pop_tesa_findpulse( EEG, 'Fz', 'refract', 4, 'rate', 2e5, 'tmsLabel', 'single','plots','off' ); %user defined input
%   EEG = pop_tesa_findpulse( EEG, 'Cz', 'paired', 'yes', 'ISI', [100],'pairLabel', {'LICI'}); %paired pulse use
%   EEG = pop_tesa_findpulse( EEG, 'Cz', 'repetitive', 'yes', 'ITI', 26, 'pulseNum', 40 ); %rTMS use 
%
% See also:
%   pop_tesa_findpulsepeak, pop_tesa_fixevent

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
    
    geometry = {1 [1 0.3] [1 0.3] [1 0.3] [1 0.3] [1 0.3] 1 [1 0.3] [1 0.3] 1 [1 0.3] 1 1 1 [1 0.3] [1 0.3] [1 0.3]};

    uilist = {{'style', 'text', 'string', 'Find TMS pulses','fontweight','bold'} ...
              {'style', 'text', 'string', 'Electrode for finding artifact'} ...
              {'style', 'popupmenu', 'string', chanIn, 'tag', 'interp', 'Value', chanVal } ...
              {'style', 'text', 'string', 'Refractory period - time for TMS artifact to recover (ms)'} ...
              {'style', 'edit', 'string', '3'}...
              {'style', 'text', 'string', 'Rate of change of TMS artifact - used to find artifact (uV/ms)'} ...
              {'style', 'edit', 'string', '1e4'}...
              {'style', 'text', 'string', 'Label for TMS pulse (single).'} ...
              {'style', 'edit', 'string', 'TMS'}...
              {'style', 'text', 'string', 'Plot identified artifacts'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 1 'tag' 'pair' } ...
              {} ...
              {'style', 'text', 'string', 'Paired pulse TMS','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Interstimulus interval (ms) [required]'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', '     Multiple ISIs can be entered as follows: 2,15,100','fontangle','italic'} ...
              {'style', 'text', 'string', 'Label for paired pulses (e.g. SICI) [default = TMSpair]'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', '     Multiple labels entered as follows: SICI, ICF, LICI','fontangle','italic'} ...
              {'style', 'text', 'string', '     Number of labels must equal number of ISIs','fontangle','italic'}...
              {} ...
              {'style', 'text', 'string', 'Repetitive TMS','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'rep' } ...
              {'style', 'text', 'string', 'Inter-train interval  (ms) [required]'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'Number of pulses in train [required]'} ...
              {'style', 'edit', 'string', ''}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Find TMS pulses -- pop_tesa_findpulse()', 'helpcom', 'pophelp(''pop_tesa_findpulse'')');
    if isempty(result), return; end;
    
    %Check that both paired and repetitive are not on
    if result{1,6}==1 && result{1,9} == 1
        error('pop_tesa_findpulse can not search for both paired and repetitive stimuli within the same file. Please choose one.');
    end
    
    %Extract data for single pulse artifact find
    elec = chanAll{result{1,1},1};
    refract = str2num(result{1,2});
    rate = str2num(result{1,3});
    tmsLabel = result{1,4};
    if result{1,5} == 1
        plots = 'on';
    else
        plots = 'off';
    end
    
    %Check if correct information is provided
    if isempty(elec)
        error('Electrode name not entered - this is required to find artifact. Script terminated')       
    end
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
        pairLabel = strtrim(strsplit(result{1,8},','));
        pairString = ['''',result{1,8},''''];
        pairString = strrep(pairString,',',''',''');
    end
    
    %Check for repetitive option
    if result{1,9} == 1 %repetitive on
        repetitive = 'yes';
        ITI = str2num(result{1,10});
        pulseNum = str2num(result{1,11});
    end

end

%Run script from input
if nargin == 2;
    EEG = tesa_findpulse(EEG,elec);
    com = sprintf('%s = pop_tesa_findpulse( %s, ''%s'' );', inputname(1), inputname(1), elec );
elseif nargin > 2
    EEG = tesa_findpulse(EEG,elec,varargin{:});
    com = sprintf('%s = pop_tesa_findpulse( %s, ''%s'', %s );', inputname(1), inputname(1), elec, vararg2str(varargin) );
end
    
if nargin < 2
    %find artifact and return the string command using pop window info
    if result{1,6}==0 && result{1,9} == 0
        EEG = tesa_findpulse( EEG, elec, 'refract', refract, 'rate', rate, 'tmsLabel', tmsLabel, 'plots', plots);
        com = sprintf('%s = pop_tesa_findpulse( %s, ''%s'', ''refract'', %s, ''rate'', %s, ''tmsLabel'', ''%s'', ''plots'', ''%s'');', inputname(1), inputname(1), elec, mat2str(refract), mat2str(rate), tmsLabel, plots);
    elseif result{1,6}==1
        if strcmp(pairLabel{1,1},'')
            EEG = tesa_findpulse( EEG, elec, 'refract', refract, 'rate', rate, 'tmsLabel', tmsLabel, 'plots', plots, 'paired', paired, 'ISI', ISI );
            com = sprintf('%s = pop_tesa_findpulse( %s, ''%s'', ''refract'', %s, ''rate'', %s, ''tmsLabel'', ''%s'', ''plots'', ''%s'', ''paired'', %s, ''ISI'', %s);', inputname(1), inputname(1), elec, mat2str(refract), mat2str(rate), tmsLabel, plots, paired, mat2str(ISI));
        else
            EEG = tesa_findpulse( EEG, elec, 'refract', refract, 'rate', rate, 'tmsLabel', tmsLabel, 'plots', plots, 'paired', paired, 'ISI', ISI, 'pairLabel', pairLabel );
            com = sprintf('%s = pop_tesa_findpulse( %s, ''%s'', ''refract'', %s, ''rate'', %s, ''tmsLabel'', ''%s'', ''plots'', ''%s'', ''paired'', ''%s'', ''ISI'', %s, ''pairLabel'', {%s});', inputname(1), inputname(1), elec, mat2str(refract), mat2str(rate), tmsLabel, plots, paired, mat2str(ISI), pairString);
        end
    elseif result{1,9}==1
        EEG = tesa_findpulse( EEG, elec, 'refract', refract, 'rate', rate, 'tmsLabel', tmsLabel, 'plots', plots, 'repetitive', repetitive, 'ITI', ITI, 'pulseNum', pulseNum );
        com = sprintf('%s = pop_tesa_findpulse( %s, ''%s'', ''refract'', %s, ''rate'', %s, ''tmsLabel'', ''%s'', ''plots'', ''%s'', ''repetitive'', ''%s'', ''ITI'', %s, ''pulseNum'', %s);', inputname(1), inputname(1), elec, mat2str(refract), mat2str(rate), tmsLabel, plots, repetitive, mat2str(ITI), mat2str(pulseNum));
    end
end

end
