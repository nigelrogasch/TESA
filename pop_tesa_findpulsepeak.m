% pop_tesa_findpulsepeak() - finds TMS pulses by detecting the large TMS artifacts
%                   peaks in the data. This script works by extracting a
%                   single channel and finding the time points in which a peak
%                   above a certain threshold is detected. Different
%                   methods for setting the threshold can be determined and
%                   either the positive or negative peak used (or an
%                   interactive gui which allows the user to select which
%                   peaks are included).
%                   Paired pulses and repetitive TMS trains can
%                   also be deteceted.
%                   This script is an alternative to tesa_findpulse
% 
%                   Note that this script requires a toolbox developed by
%                   Daniel Wagenaar (WagenaarMBL) which is available for
%                   download at: http://www.its.caltech.edu/~daw/teach.html
%                   Once the toolbox is unzipped, please make sure this is
%                   added to your matlab path.
%
% Usage:
%   >>  EEG = pop_tesa_findpulsepeak( EEG );
%   >>  EEG = pop_tesa_findpulsepeak( EEG, elec );
%   >>  EEG = pop_tesa_findpulsepeak( EEG, elec, 'key1', value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   elec            - string with electrode to use for finding artifact
% 
% Optional input pairs:
%   'dtrnd','str'  'poly'|'gradient'|'median'|'linear'|'none'. Defines the type of detrend used
%                   to centre the data, to generate an analytic signal.
%                   default = poly
%   'thrshtype','str'/int - 'dynamic'|'median'|value. Defines the type of
%                   threshold used to determine peaks. Dynamic sets threshold 
%                   to the range points above/below 99.9 percent of data
%                   trace. Median sets threshold as median of points above/below 
%                   99.9 percent of data trace. Value is a user defined integer 
%                   for setting the threshold (in uV). (e.g. 1000)
%                   default = 'dymanic'
%   'wpeaks','str' - 'pos'|'neg'|'gui'. Defines whether to use the
%                   positive or negative peak to define the artifact, or to
%                   use an interactive GUI. For the GUI, a box is generated
%                   by clicking and dragging the mouse cursor. Peaks within
%                   the box will be included as events.
%                   default = 'pos'
%   'plots','str' - 'on'|'off'. Brings up a plot showing the detected
%                   peaks. Black = detected, pink = selected for
%                   definition.
%                   default = 'on'
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
% Examples
%   EEG = pop_tesa_findpulsepeak( EEG, 'Cz' ); %default use
%   EEG = pop_tesa_findpulsepeak( EEG, 'Fz', 'dtrnd', 'linear', 'thrshtype', 'median', 'wpeaks', 'gui', 'plots', 'off', 'tmsLabel', 'single' ); %user defined
%   EEG = pop_tesa_findpulsepeak( EEG, 'Cz', 'paired', 'yes', 'ISI', [100],'pairLabel', {'LICI'}); %paired pulse use
%   EEG = pop_tesa_findpulsepeak( EEG, 'Cz', 'repetitive', 'yes', 'ITI', 26, 'pulseNum', 40 ); %rTMS use 
%
% See also:
%   pop_tesa_findpulse, pop_tesa_fixevent 

% Copyright (C) 2016  Nigel Rogasch, Monash University,
% nigel.rogasch@monash.edu
%
% Authors:
% Caley Sullivan, Monash University, calley.sullivan@monash.edu
% 
% Based on functions developed by Daniel Wagenaar 
%                                  http://www.its.caltech.edu/~daw/teach.html
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

 function [EEG, com] = pop_tesa_findpulsepeak( EEG, elec, varargin )

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
    
    geometry = {1 [1 0.3] [1 0.3] [1 0.3] [1 0.3] [1 0.3] [1 0.3] [1 0.3] 1 [1 0.3] [1 0.3] 1 [1 0.3] 1 1 1 [1 0.3] [1 0.3] [1 0.3]};

    uilist = {{'style', 'text', 'string', 'Find TMS pulses','fontweight','bold'} ...
              { 'Style', 'text', 'string', 'Electrode for finding artifact' }, ...
            { 'Style', 'edit', 'string', chan, 'tag', 'chans' }, ...
            {},  ...
            { 'style' 'pushbutton' 'string'  'select', 'enable' fastif(isempty(EEG.chanlocs), 'off', 'on') ...
            'callback' 'tmpchanlocs = EEG(1).chanlocs; [tmp tmpval] = tesa_chansel({tmpchanlocs.labels}, ''withindex'', ''on'', ''selectionmode'',''single''); set(findobj(gcbf, ''tag'', ''chans''), ''string'',tmpval); clear tmp tmpchanlocs tmpval' }, ...
              {'style', 'text', 'string', 'Type of detrend'} ...
              {'style', 'popupmenu', 'string', 'poly|gradient|median|linear|none' 'tag' 'detrend'}...
              {'style', 'text', 'string', 'Type of thresholding'} ...
              {'style', 'popupmenu', 'string', 'dynamic|median|manual' 'tag' 'thresh'}...
              {'style', 'text', 'string', 'Peak to define artifact'} ...
              {'style', 'popupmenu', 'string', 'positive|negative|gui' 'tag' 'peak'}...
              {'style', 'text', 'string', 'Sanity check plot'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 1 'tag' 'plot' } ...
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
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Find TMS pulses -- pop_tesa_findpulsepeak()', 'helpcom', 'pophelp(''pop_tesa_findpulsepeak'')');
    if isempty(result), return; end;
    
    %Check that both paired and repetitive are not on
    if result{1,7} == 1 && result{1,10} == 1
        error('tesa_findpulsepeak can not search for both paired and repetitive stimuli within the same file. Please choose one.');
    end
    
    %Extract data for single pulse artifact find
    elec = result{1,1};
    
    if result{1,2} == 1;
        dtrnd = 'poly';
    elseif result{1,2} == 2;
        dtrnd = 'gradient';
    elseif result{1,2} == 3;
        dtrnd = 'median';
    elseif result{1,2} == 4;
        dtrnd = 'linear';
    elseif result{1,2} == 5;
        dtrnd = 'none';
    end
    
    if result{1,3} == 1;
        thrshtype = 'dynamic';
    elseif result{1,3} == 2;
        thrshtype = 'median';
    elseif result{1,3} == 3;
        geometry = {[1 0.3]};
        uilist = {{'style', 'text', 'string', 'Threshold for artifact detection (uV)'} ...
                    {'style', 'edit', 'string', '1000'}};
        result2 = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Manual threshold for artifact');       
        thrshtype = str2num(result2{1,1});
    end
    
    if result{1,4} == 1
        wpeaks = 'pos';
    elseif result{1,4} == 2
        wpeaks = 'neg';
    elseif result {1,4} == 3
        wpeaks = 'gui';
    end
    
    if result{1,5} == 1
        plots = 'on';
    elseif result{1,5} == 0
        plots = 'off';
    end
        
    tmsLabel = result{1,6};
    
    %Check if correct information is provided
    if isempty(elec)
        error('Electrode name not entered - this is required to find artifact. Script terminated')       
    end

    if isempty(tmsLabel)
        tmsLabel = 'TMS';      
    end
    
    %Check for paired option
    if result{1,7} == 1 %paired on
        paired = 'yes';
        ISI = str2num(result{1,8});
        pairLabel = strtrim(strsplit(result{1,9},','));
        pairString = ['''',result{1,9},''''];
        pairString = strrep(pairString,',',''',''');
    end
    
    %Check for repetitive option
    if result{1,10} == 1 %repetitive on
        repetitive = 'yes';
        ITI = str2num(result{1,11});
        pulseNum = str2num(result{1,12});
    end

end
    chan=mat2str(elec);
    thrshtype1 = mat2str(thrshtype);
%Run script from input
if nargin == 2;
    EEG = tesa_findpulsepeak(EEG,elec);
    com = sprintf('%s = tesa_findpulsepeak( %s, %s );', inputname(1), inputname(1),chan);
elseif nargin > 2
    EEG = tesa_findpulsepeak(EEG,elec,varargin{:});
    com = sprintf('%s = tesa_findpulsepeak( %s, %s, %s );', inputname(1), inputname(1),chan, vararg2str(varargin) );
end
    

%find artifact and return the string command using pop window info
if nargin < 2
    if result{1,7} == 0 && result{1,10} == 0
        EEG = tesa_findpulsepeak( EEG, elec, 'dtrnd', dtrnd, 'thrshtype', thrshtype, 'wpeaks', wpeaks, 'plots', plots, 'tmsLabel', tmsLabel );
%         if ~ischar(thrshtype)
%             thrshtype1 = mat2str(thrshtype);
%         else
%             thrshtype1 = thrshtype;
%         end
        com = sprintf('%s = tesa_findpulsepeak( %s, %s, ''dtrnd'', ''%s'', ''thrshtype'', %s, ''wpeaks'', ''%s'', ''plots'', ''%s'', ''tmsLabel'', ''%s'');', inputname(1), inputname(1), chan, dtrnd, thrshtype1, wpeaks, plots, tmsLabel);
    elseif result{1,7} == 1
        if strcmp(pairLabel{1,1},'')
            EEG = tesa_findpulsepeak( EEG, elec, 'dtrnd', dtrnd, 'thrshtype', thrshtype, 'wpeaks', wpeaks, 'plots', plots, 'tmsLabel', tmsLabel, 'paired', paired, 'ISI', ISI );
%             if ~ischar(thrshtype)
%                 thrshtype1 = mat2str(thrshtype);
%             else
%                 thrshtype1 = thrshtype;
%             end
            com = sprintf('%s = tesa_findpulsepeak( %s, %s, ''dtrnd'', ''%s'', ''thrshtype'', %s, ''wpeaks'', ''%s'', ''plots'', ''%s'', ''tmsLabel'', ''%s'', ''paired'', ''%s'', ''ISI'', %s);', inputname(1), inputname(1), chan, dtrnd, thrshtype1, wpeaks, plots, tmsLabel, paired, mat2str(ISI));

        else
            EEG = tesa_findpulsepeak( EEG, elec, 'dtrnd', dtrnd, 'thrshtype', thrshtype, 'wpeaks', wpeaks, 'plots', plots, 'tmsLabel', tmsLabel, 'paired', paired, 'ISI', ISI, 'pairLabel', pairLabel);
%             if ~ischar(thrshtype)
%                 thrshtype1 = mat2str(thrshtype);
%             else
%                 thrshtype1 = thrshtype;
%             end
            com = sprintf('%s = tesa_findpulsepeak( %s, %s, ''dtrnd'', ''%s'', ''thrshtype'', %s, ''wpeaks'', ''%s'', ''plots'', ''%s'', ''tmsLabel'', ''%s'', ''paired'', ''%s'', ''ISI'', %s, ''pairLabel'', {%s});', inputname(1), inputname(1), chan, dtrnd, thrshtype1, wpeaks, plots, tmsLabel, paired, mat2str(ISI), pairString);
        end
    elseif result{1,10} == 1
        EEG = tesa_findpulsepeak( EEG, elec, 'dtrnd', dtrnd, 'thrshtype', thrshtype, 'wpeaks', wpeaks, 'plots', plots, 'tmsLabel', tmsLabel, 'repetitive',repetitive, 'ITI', ITI, 'pulseNum', pulseNum);
%         if ~ischar(thrshtype)
%             thrshtype1 = mat2str(thrshtype);
%         else
%             thrshtype1 = thrshtype;
%         end
        com = sprintf('%s = tesa_findpulsepeak( %s, %s, ''dtrnd'', ''%s'', ''thrshtype'', %s, ''wpeaks'', ''%s'', ''plots'', ''%s'', ''tmsLabel'', ''%s'', ''repetitive'', ''%s'', ''ITI'', %s, ''pulseNum'', %s);', inputname(1), inputname(1), chan, dtrnd, thrshtype1, wpeaks, plots, tmsLabel, repetitive, mat2str(ITI), mat2str(pulseNum));
    end
end

end