% pop_tesa_findpulse() - finds TMS pulses by detecting the large TMS artifacts
%                   present in the data. This script works by extracting a
%                   single channel and finding the time points in which the 
%                   first derivatives exceed a certain threshold (defined 
%                   by 'rate'). Paired pulses and repetitive TMS trains can
%                   also be deteceted. 
%
% Usage:
%   >>  EEG = pop_tesa_findpulse( EEG );
%   >>  EEG = pop_tesa_findpulse( EEG, elec, varargin );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   elec            - string with electrode to use for finding artifact
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
    
    geometry = {1 [1 0.3] [1 0.3] [1 0.3] [1 0.3] 1 [1 0.3] [1 0.3] 1 [1 0.3] 1 1 1 [1 0.3] [1 0.3] [1 0.3]};

    uilist = {{'style', 'text', 'string', 'Find TMS pulses','fontweight','bold'} ...
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
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Find TMS pulses -- pop_tesa_findpulse()', 'helpcom', 'pophelp(''tesa_findpulse'')');
    
    %Check that both paired and repetitive are not on
    if result{1,5}==1 && result{1,8} == 1
        error('tesa_findpulse can not search for both paired and repetitive stimuli within the same file. Please choose one.');
    end
    
    %Extract data for single pulse artifact find
    elec = result{1,1};
    refract = str2num(result{1,2});
    rate = str2num(result{1,3});
    tmsLabel = result{1,4};
    
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
    if result{1,5} == 1 %paired on
        paired = 'yes';
        ISI = str2num(result{1,6});
        pairLabel = strtrim(strsplit(result{1,7},','));
    end
    
    %Check for repetitive option
    if result{1,8} == 1 %repetitive on
        repetitive = 'yes';
        ITI = str2num(result{1,9});
        pulseNum = str2num(result{1,10});
    end

end

%Run script from input
if nargin == 2;
    EEG = tesa_findpulse(EEG,elec);
    com = sprintf('%s = pop_tesa_findpulse( %s, %s );', inputname(1), inputname(1), elec );
elseif nargin > 2
    EEG = tesa_findpulse(EEG,elec,varargin{:});
    com = sprintf('%s = pop_tesa_findpulse( %s, %s, %s );', inputname(1), inputname(1), elec, vararg2str(varargin) );
end
    
if nargin < 2
    %find artifact and return the string command using pop window info
    if result{1,5}==0 && result{1,8} == 0
        EEG = tesa_findpulse( EEG, elec, 'refract', refract, 'rate', rate, 'tmsLabel', tmsLabel );
        com = sprintf('%s = pop_tesa_findpulse( %s, %s, ''refract'', %s, ''rate'', %s, ''tmsLabel'', %s);', inputname(1), inputname(1), elec, mat2str(refract), mat2str(rate), tmsLabel);
    elseif result{1,5}==1
        if strcmp(pairLabel{1,1},'')
            EEG = tesa_findpulse( EEG, elec, 'refract', refract, 'rate', rate, 'tmsLabel', tmsLabel, 'paired', paired, 'ISI', ISI );
            com = sprintf('%s = pop_tesa_findpulse( %s, %s, ''refract'', %s, ''rate'', %s, ''tmsLabel'', %s, ''paired'', %s, ''ISI'', %s);', inputname(1), inputname(1), elec, mat2str(refract), mat2str(rate), tmsLabel, paired, mat2str(ISI));
        else
            EEG = tesa_findpulse( EEG, elec, 'refract', refract, 'rate', rate, 'tmsLabel', tmsLabel, 'paired', paired, 'ISI', ISI, 'pairLabel', pairLabel );
            com = sprintf('%s = pop_tesa_findpulse( %s, %s, ''refract'', %s, ''rate'', %s, ''tmsLabel'', %s, ''paired'', %s, ''ISI'', %s, ''pairLabel'', {%s});', inputname(1), inputname(1), elec, mat2str(refract), mat2str(rate), tmsLabel, paired, mat2str(ISI), result{1,7});
        end
    elseif result{1,8}==1
        EEG = tesa_findpulse( EEG, elec, 'refract', refract, 'rate', rate, 'tmsLabel', tmsLabel, 'repetitive', repetitive, 'ITI', ITI, 'pulseNum', pulseNum );
        com = sprintf('%s = pop_tesa_findpulse( %s, %s, ''refract'', %s, ''rate'', %s, ''tmsLabel'', %s, ''repetitive'', %s, ''ITI'', %s, ''pulseNum'', %s);', inputname(1), inputname(1), elec, mat2str(refract), mat2str(rate), tmsLabel, repetitive, mat2str(ITI), mat2str(pulseNum));
    end
end

end
