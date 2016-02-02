% pop_ tesa_tepanalysis() - averages over trials and electrodes in a region of
%                   interest to extract TMS-evoked potentials. Outputs are
%                   saved in the EEG structure (EEG.ROI). For TEPs
%                   following paired pulses, an additional file can be
%                   specified to subtract from the conditioning pulse and
%                   thereby minimising the impact of on-going activity in  the test
%                   pulse time period.
%
% Usage:
%   >>  EEG = pop_tesa_roianalysis( EEG ); %for pop up window
%   >>  EEG = pop_tesa_roianalysis( EEG, roi );
%   >>  EEG = pop_tesa_roianalysis( EEG, roi, varargin );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   roi             - string or cell array defining electrodes to be
%                   averaged for ROI analysis. Type 'all' to average over
%                   all electrodes.
%                   Examples: 'C3' (single); {'C3','C4','CP1'} (multiple); 
%                       'all' (all electrodes).
% 
% Optional input pairs (varargin):
%   'roiName','str' - 'str' is a name to identify ROI analysis. This is
%                   useful if multiple different ROIs are to be analysed.
%                   The output will be stored under this name in EEG.ROI.
%                   Example: 'motor'
%                   Defaults are: R1,R2,R....
% 
% Optional input pairs for correcting paired pulses (varargin):
%   'pairCorrect','on' - turns on pair correction.
%   'ISI',int       - int is a number which defines the inter-stimulus
%                   interval between the conditioning and test pulse. The
%                   single TEP that is subtracted from the paired TEP will
%                   be shifted by this many ms to align with the
%                   conditioning pulse. int is in ms.
%                   Example: 100
%   'fileName','str' - 'str' is the path and name of the .set file to be
%                   subtracted from the paired TEP. This should be a single
%                   TEP evoked by a stimulus intensity equivalent to the
%                   conditioning pulse.
%                   Example: 'C:\tmseeg\myfile.set'
% 
%   Further reading on importance of correcting paired pulses with TMS-EEG:
%   Rogasch N.C. et al (2015) Cortical inhibition of distinct mechanisms in 
%   the dorsolateral prefrontal cortex is related to working memory performance:    
%   A TMS–EEG study. Cortex, 64:68-77.
%   In particular, see supplementary materials
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

 function [EEG, com] = pop_tesa_tepextract( EEG, type, varargin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
    
    geometry = {1 [1 0.5] 1 1 [1 0.5] 1 [1 0.5] 1 1 [1 0.5] [1 0.5] 1 1};

    uilist = {{'style', 'text', 'string', 'TMS-evoked potential analysis','fontweight','bold'} ...
              {'style', 'text', 'string', 'Select type of analysis'} ...
              {'style', 'popupmenu', 'string', 'ROI|GMFA' 'tag' 'type' }...
              {}...
              {'style', 'text', 'string', 'For ROI analysis','fontweight','bold'} ...
              {'style', 'text', 'string', 'Electrodes for ROI analysis [required for ROI]'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', '     Examples: C3 (one); C3,FC1,CP1 (multiple); all (all)','fontangle','italic'} ...
              {'style', 'text', 'string', 'Identifier for ROI analysis (optional)'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', '     Example: motor', 'fontangle','italic'} ...
              {} ...
              {'style', 'text', 'string', 'Paired pulse correction (if required)','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'Yes?' 'value' 0 'tag' 'pair' } ...
              {'style', 'text', 'string', 'Interstimulus interval (ms) [required]'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', '     Example: 100','fontangle','italic'}...
              {'style', 'text', 'string', 'Press OK, choose file for subtraction'}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'ROI analysis -- pop_tesa_roianalysis()', 'helpcom', 'pophelp(''tesa_roianalysis'')');
      
    %Extract data for single pulse artifact find
    if result{1,1} == 1
        type = 'ROI';
    elseif result{1,1} == 2
        type = 'GMFA';
    end    
    roi = strtrim(strsplit(result{1,2},','));
    roiName = result{1,3};
    ISI = str2num(result{1,5});
    
    %Check if roi provided
    if strcmp(roi,'')
        error('No electrodes provided for ROI analysis. Please provide at least one electrode.');
    end
    
    %If paired correction is on
    if result{1,4} == 1;
        pairCorrect = 'on';
        [name1, inputpath] = uigetfile2('*.SET*;*.set', 'Choose file for subtraction');
        fileName = [inputpath,name1];
    end
    
    %Check if correct information is provided
    if isempty(roi)
        error('Electrode name not entered - this is required for ROI analysis. Script terminated')       
    end
    
end

%Run script from input
if nargin == 2 && strcmp(type,'GMFA');
    EEG = tesa_tepextract(EEG,type);
    com = sprintf('%s = pop_tesa_tepextract( %s, %s );', inputname(1), inputname(1), type );
elseif nargin > 2
    EEG = tesa_tepextract(EEG,type,varargin{:});
    com = sprintf('%s = pop_tesa_tepextract( %s, %s, %s );', inputname(1), inputname(1), type, vararg2str(varargin) );
end

%perform roi analysis and return the string command
if nargin < 2
    if result{1,4} == 0 
        EEG = tesa_tepextract( EEG, type, 'roi', roi, 'roiName', roiName );
        com = sprintf('%s = pop_tesa_tepextract( %s, %s, ''roi'', %s, ''roiName'', %s );', inputname(1), inputname(1), result{1,2}, result{1,3}, roiName);
    elseif result{1,4} == 1
        EEG = tesa_tepextract( EEG, type, 'roi', roi, 'roiName', roiName, 'pairCorrect', pairCorrect, 'ISI', ISI, 'fileName', fileName );
        com = sprintf('%s = pop_tesa_tepextract( %s, %s, ''roi'',{%s}, ''roiName'', %s, ''pairCorrect'', %s, ''ISI'', %s, ''fileName'', %s );', inputname(1), inputname(1), result{1,2}, result{1,3}, roiName, pairCorrect, mat2str(ISI), fileName);
    end
end

end
