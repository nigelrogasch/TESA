% pop_tesa_tepextract() - averages over trials to generate a TMS-evoked potential (TEP). 
%                   Either a region-of-interest analysis, which averages
%                   over selected electrodes, or a global field amplitude
%                   analysis (standard deviation across electrodes at each
%                   time point) can be performed. Outputs are
%                   saved in the EEG structure EEG.ROI or EEG.GMFA respectively. 
%                   For TEPs following paired pulses, an additional file can be
%                   specified to subtract from the conditioning pulse and
%                   thereby minimising the impact of on-going activity in  the test
%                   pulse time period.
% 
%                   Finding peaks and returning amplitudes and latencies is
%                   performed with tesa_peakanalysis and tesa_peakoutput.
%
%                   Further reading on importance of correcting paired pulses with TMS-EEG:
%                   'Rogasch N.C. et al (2015) Cortical inhibition of distinct mechanisms in 
%                   the dorsolateral prefrontal cortex is related to working memory performance:    
%                   A TMS–EEG study. Cortex, 64:68-77.'
%                   In particular, see supplementary materials
%
% Usage:
%   >>  EEG = pop_tesa_tepextract( EEG ); %pop up window
%   >>  EEG = pop_tesa_tepextract( EEG, type );
%   >>  EEG = pop_tesa_tepextract( EEG, type, 'key1',value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   type            - string for analysis type. Either 'ROI' or 'GMFA'
%   'elecs','{}'    - [required for type 'ROI'] string or 
%                   cell array defining electrodes to be
%                   averaged for ROI analysis. Type 'all' to average over
%                   all electrodes.
%                   Examples: 'elecs','C3' (single); 'elecs',{'C3','C4','CP1'} (multiple); 
%                   'elecs','all' (all electrodes).
% 
% Optional input pairs (varargin):
%   'tepName','str' - 'str' is a name to identify analysis. This is
%                   useful if multiple different analyses (i.e. ROIs) are to be analysed.
%                   The output will be stored under this name in EEG.ROI or
%                   EEG.GMFA. Note this needs to be a continuous word -
%                   use an underscore ( _ ) instead of a space. 
%                   Example: 'motor', 'parietal', 'right_frontal'
%                   Defaults are: R1,R2,R....
% 
% Optional input pairs for correcting paired pulses (varargin):
%   'pairCorrect','on' - turns on pair correction.
%   'ISI',int       - [required if pairCorrect on] int is an integer which defines 
%                   the inter-stimulus interval between the conditioning and test pulse. 
%                   The single TEP that is subtracted from the paired TEP will
%                   be shifted by this many ms to align with the
%                   conditioning pulse. int is in ms.
%                   Example: 100
%   'fileName','str' - [required if pairCorrect on]'str' is the path and name of the .set file 
%                   to be subtracted from the paired TEP. This should be a single
%                   TEP evoked by a stimulus intensity equivalent to the
%                   conditioning pulse.
%                   Example: 'C:\tmseeg\myfile.set'
%
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples
%   EEG = pop_tesa_tepextract( EEG, 'ROI', 'elecs', {'FC1','FC3','C1','C3'} ); % standard ROI analysis
%   EEG = pop_tesa_tepextract( EEG, 'ROI', 'elecs', {'C1','C3'}, 'tepName','motor' ); % ROI analysis with specific name
%   EEG = pop_tesa_tepextract( EEG, 'ROI', 'elecs', 'C3', 'pairCorrect', 'on', 'ISI', 100, 'fileName', 'C:\tmseeg\LICI.set' ); % paired pulse analysis
%   EEG = pop_tesa_tepextract( EEG, 'GMFA'); % standard GMFA analysis
%   EEG = pop_tesa_tepextract( EEG, 'GMFA', 'pairCorrect', 'on', 'ISI', 100, 'fileName', 'C:\tmseeg\LICI.set' ); % paired pulse analysis for GMFA
% 
% See also:
%   pop_tesa_peakanalysis, pop_tesa_peakoutput 

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

 function [EEG, com] = pop_tesa_tepextract( EEG, type, varargin )

com = '';          

% %check that data is present
% if isempty(EEG.data)
%     error('Data is empty');
% end

% pop up window
% -------------
if nargin < 2
    
    geometry = {1 [0.5 0.4] 1 1 [0.5 0.4] [0.5 0.4] 1 [0.5 0.4] 1 [0.5 0.4] [0.5 0.4] 1};

    uilist = {{'style', 'text', 'string', 'TMS-evoked potential analysis','fontweight','bold'} ...
              {'style', 'text', 'string', 'Select type of analysis'} ...
              {'style', 'popupmenu', 'string', 'Region of interest|Global mean field amplitude' 'tag' 'type' }...
              {}...
              {'style', 'text', 'string', 'For ROI analysis','fontweight','bold'} ...
              {'style', 'text', 'string', 'Electrodes for ROI analysis (e.g. C1 C3)'} ...
              {'style', 'edit', 'string', '', 'tag', 'chans' } ...
              {'style', 'text', 'string', 'Required for ROI analysis','fontAngle','italic'}...
              {'style', 'pushbutton', 'string',  'Select electrodes', 'enable', 'on' ...
               'callback', 'tmpchanlocs = EEG(1).chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on''); set(findobj(gcbf, ''tag'', ''chans''), ''string'',tmpval); clear tmp tmpchanlocs tmpval' }, ...
              {}...
              {'style', 'text', 'string', 'Identifier for ROI analysis (e.g. motor)'} ...
              {'style', 'edit', 'string', ''}...
              {} ...
              {'style', 'text', 'string', 'Paired pulse correction (if required)','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Interstimulus interval (ms) (e.g. 100)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'If on, press OK to choose file for subtraction','fontAngle','italic'}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Extract TEPs -- pop_tesa_tepextract()', 'helpcom', 'pophelp(''pop_tesa_tepextract'')');
    if isempty(result), return; end;
    
    %Extract data for single pulse artifact find
    if result{1,1} == 1
        type = 'ROI';
    elseif result{1,1} == 2
        type = 'GMFA';
    end
    
    elecs = strtrim(strsplit(result{1,2},' '));
    elecsString = ['''',result{1,2},''''];
    elecsString = strrep(elecsString,' ',''',''');
    elecsString = strrep(elecsString,' ','');
    tepName = result{1,3};
    ISI = str2num(result{1,5});
    
    %If paired correction is on
    if result{1,4} == 1;
        pairCorrect = 'on';
        [name1, inputpath] = uigetfile2('*.SET*;*.set', 'Choose file for subtraction');
        fileName = [inputpath,name1];
    end
    
    %Check if correct information is provided
    if strcmp(type,'ROI') && strcmp(elecs{1,1},'')
        error('Electrode name not entered - this is required for ROI analysis.')       
    end
    
end

%Run script from input
if nargin == 2 
    EEG = tesa_tepextract(EEG,type);
    com = sprintf('%s = pop_tesa_tepextract( %s, ''%s'' );', inputname(1), inputname(1), type );
elseif nargin > 2
    EEG = tesa_tepextract(EEG,type,varargin{:});
    com = sprintf('%s = pop_tesa_tepextract( %s, ''%s'', %s );', inputname(1), inputname(1), type, vararg2str(varargin) );
end

%perform roi analysis and return the string command
if nargin < 2
    if strcmp(type,'GMFA') && result{1,4} == 0
        EEG = tesa_tepextract(EEG,type);
        com = sprintf('%s = pop_tesa_tepextract( %s, ''%s'' );', inputname(1), inputname(1), type );
    elseif strcmp(type,'GMFA') && result{1,4} == 1
        EEG = tesa_tepextract(EEG,type,'pairCorrect',pairCorrect,'ISI',ISI,'fileName',fileName);
        com = sprintf('%s = pop_tesa_tepextract( %s, ''%s'', ''pairCorrect'', ''%s'', ''ISI'', %s, ''fileName'', ''%s'');', inputname(1), inputname(1), type, pairCorrect, mat2str(ISI), fileName );
    elseif result{1,4} == 0 && strcmp(tepName,'')
        EEG = tesa_tepextract( EEG, type, 'elecs', elecs);
        com = sprintf('%s = pop_tesa_tepextract( %s, ''%s'', ''elecs'', {%s} );', inputname(1), inputname(1), type, elecsString);
    elseif result{1,4} == 0
        EEG = tesa_tepextract( EEG, type, 'elecs', elecs, 'tepName', tepName );
        com = sprintf('%s = pop_tesa_tepextract( %s, ''%s'', ''elecs'', {%s}, ''tepName'', ''%s'' );', inputname(1), inputname(1), type, elecsString, tepName);
    elseif result{1,4} == 1 && strcmp(tepName,'')
        EEG = tesa_tepextract( EEG, type, 'elecs', elecs, 'pairCorrect', pairCorrect, 'ISI', ISI, 'fileName', fileName );
        com = sprintf('%s = pop_tesa_tepextract( %s, ''%s'', ''elecs'',{%s}, ''pairCorrect'', ''%s'', ''ISI'', %s, ''fileName'', ''%s'');', inputname(1), inputname(1), type, elecsString, pairCorrect, mat2str(ISI), fileName);
    elseif result{1,4} == 1
        EEG = tesa_tepextract( EEG, type, 'elecs', elecs, 'tepName', tepName, 'pairCorrect', pairCorrect, 'ISI', ISI, 'fileName', fileName );
        com = sprintf('%s = pop_tesa_tepextract( %s, ''%s'', ''elecs'',{%s}, ''tepName'', ''%s'', ''pairCorrect'', ''%s'', ''ISI'', %s, ''fileName'', ''%s'' );', inputname(1), inputname(1), type, elecsString, tepName, pairCorrect, mat2str(ISI), fileName);
    end
end

end
