% pop_tesa_peakanalysis() - finds peaks within a time window defined by user 
%                       for either ROI or GMFA analyses.
%                       Either positive or negative peaks are detected.
%                       Peaks are defined as data point which is
%                       larger/smaller than +/- x data points (default = 5,
%                       however this can also be defined by user). Results
%                       are saved in EEG structure under (either EEG.ROI or
%                       EEG.GMFA). If no peak is found in the defined window, 
%                       the amplitude at the latency defined in peak is returned,
%                       and a NaN is returned in latency. The analysis is run on 
%                       all existing outputs from tesa_tepextract (e.g. ROIs or GMFA), unless
%                       the user opts to run the analysis on one specific
%                       ROI or GMFA.
% 
%                       Note that tesa_tepextract must be run prior to this script.
%
% Usage:
%   >>  EEG = pop_tesa_peakanalysis( EEG ); % pop up window
%   >>  EEG = pop_tesa_peakanalysis( EEG, input, direction, peak, peakWin );
%   >>  EEG = pop_tesa_peakanalysis( EEG, input, direction, peak, peakWin, 'key1', value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   input           - string defining whether 'ROI' or 'GMFA' input is used
%   direction       - string defining whether peaks are 'positive' or
%                       'negative'. Use 'positive' for GMFA analysis.
%   peak            - vector defining the peak of interest. 
%                       Example: 25 (one peak), [25,60,..] (multiple peaks)
%   peakWin         - matrix defining the time windows to search for above
%                       peaks. Minimum and maximum values for time window are 
%                       defined as 15,35 and peak definitions separated by ;
%                       Note that the number of peak windows defined should
%                       equal the number of peaks. 
%                       For example if peak = [25,60,180]; peakWin = [15,35;40,80;160,200]
% 
% Optional input pairs
%   'method','str'  - either 'largest' or 'centre' (default = largest). If multiple
%                       peaks are detected in a window, largest will search for
%                       the largest peak within the time window. Centre will 
%                       search for the peak closest to the latency defined
%                       in peak.
%   'samples',int   - int is an integer defining the number of samples
%                       either side of a peak that defines the peak.
%                       Peaks are defined as data point which is
%                       larger/smaller than +/- int data points (default = 5)
%   'tepName','str' - 'str' is a name of a specific ROI to perform the
%                       analysis on. If this is left blank, all ROI/GMFAs
%                       defined by tesa_tepextract are analysed
%                       Example: 'motor'
%                       Defaults are: R1,R2,R....
%     
% Outputs:
%   EEG             - EEGLAB EEG structure
%
%   Examples
%   EEG = pop_tesa_peakanalysis( EEG, 'ROI', 'negative', 100, [80,120] ); %find a negative peak in all ROI analyses at 100 ms searching between 80 and 120 ms.
%   EEG = pop_tesa_peakanalysis( EEG, 'GMFA', 'positive', [30,60,180],[20,40;50,70;170,190] ); %find 3 positive peaks in the GMFA analysis at 30 ms (between 20-40ms), 60 ms (between 50-70 ms), and 180 ms (between 170-190 ms)
%   EEG = pop_tesa_peakanalysis( EEG, 'ROI', 'positive', [25,70], [15,35;60,80], 'method', 'centre', 'samples', 5, 'tepName', 'motor'); %find 2 positive peaks at 25 ms (15-35 ms), and 70 ms (60-80 ms) using the peak closest to the central peak (i.e. 25 ms or 70 ms), defining a peak as a data point that is larger than all data points +/- 5 samples and only for the ROI analysis named 'motor'.
% 
% See also:
%   pop_tesa_tepextract, pop_tesa_peakoutput 

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

 function [EEG com] = pop_tesa_peakanalysis( EEG, input, direction, peak, peakWin, varargin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
    
    if isfield(EEG,'ROI')
        roiFields = fieldnames(EEG.ROI);
        roi = 'yes';
        for a = 1:size(roiFields,1)
            roiFields{a,1} = ['ROI ',roiFields{a,1}];
        end
    else 
        roi = 'no';
    end
    if isfield(EEG,'GMFA')
        gmfaFields = fieldnames(EEG.GMFA);
        gmfa = 'yes';
        for a = 1:size(gmfaFields,1)
            gmfaFields{a,1} = ['GMFA ',gmfaFields{a,1}];
        end
    else 
        gmfa = 'no';
    end
    
    if strcmpi(roi,'no') && strcmpi(gmfa,'no')
        error('There are no ROI or GMFA analyses in this data. Please run tesa_tepextract first.');
    elseif strcmpi(roi,'yes') && strcmpi(gmfa,'yes')
        tepInput = ['ROI all'; roiFields; 'GMFA all'; gmfaFields];
    elseif strcmpi(roi,'yes') && strcmpi(gmfa,'no')
        tepInput = ['ROI all'; roiFields];
    elseif strcmpi(roi,'no') && strcmpi(gmfa,'yes')
        tepInput = ['GMFA all'; gmfaFields];
    end
    
    tepAll = tepInput;
    
    for a = 1:size(tepInput,1)-1
        tepInput{a,1} = [tepInput{a,1},'|'];
    end
    
    tepIn = [tepInput{:}];
      
    geometry = {1 [0.5 0.5 0.5] 1 [0.5 0.5 0.5] 1 [1 0.5] 1 [1 0.5] 1 1 1 1 [1 0.5] [1 0.5]};

    uilist = {{'style', 'text', 'string', 'Find peaks','fontweight','bold'} ...
              {'style', 'text', 'string', 'Input type for analysis'} ...
              {'style', 'popupmenu', 'string', tepIn, 'tag', 'input' } ...
              {}...
              {}...
              {'style', 'text', 'string', 'Direction of peaks'} ...
              {'style', 'popupmenu', 'string', 'positive|negative' 'tag' 'direction' } ...
              {}...
              {} ...
              {'style', 'text', 'string', 'Peaks to find (ms) [required]','fontweight','bold'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', '   Example: 25 (one peak), 25,60,... (multiple peaks)','fontangle','italic'} ...
              {'style', 'text', 'string', 'Time window for peak search (ms) [required]','fontweight','bold'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', '   Example: 15,45 - for one peak; 15-45 ms (P25).','fontangle','italic'} ...
              {'style', 'text', 'string', '   Example: 15,45; 46,75 - for multiple peaks; 15-45 ms (P25), 46-75 ms (P60).','fontangle','italic'} ...
              {} ...
              {'style', 'text', 'string', 'Optional inputs','fontweight','bold'} ...
              {'style', 'text', 'string', 'Method for determining which peak to use'} ...
              {'style', 'popupmenu', 'string', 'largest|centre' 'tag' 'input' } ...
              {'style', 'text', 'string', 'Number of samples for defining peak'} ...
              {'style', 'edit', 'string', '5'}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Find peaks -- pop_tesa_peakanalysis()', 'helpcom', 'pophelp(''pop_tesa_peakanalysis'')');
    if isempty(result), return; end;
    
    %Check that peaks have been entered
    if strcmp(result{1,3},'')
        error('No peaks were defined. Please enter peaks.');
    end
    
    %Check that peak windows have been entered
    if strcmp(result{1,4},'')
        error('No peak windows were defined. Please enter peak windows.');
    end
    
    %Extract data for peak analysis
    tepFor = tepAll{result{1,1},1};
    tepFor = strtrim(strsplit(tepFor,' ')); 
    
    input = tepFor{1,1};
    
    if result{1,2} == 1
        direction = 'positive';
    elseif result{1,2} == 2
        direction = 'negative';
    end
    
    peak = str2num(result{1,3});
    
    peakWin = str2num(result{1,4});
    
    if result{1,5} == 1
        method = 'largest';
    elseif result{1,5} == 2
        method = 'centre';
    end
    
    samples = str2num(result{1,6});
    
    if strcmpi(tepFor{1,2},'all')
        tepName = [];
    else 
        tepName = tepFor{1,2};
    end
    
end

%Run script from input
if nargin == 5
    EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin );
    com = sprintf('%s = pop_tesa_peakanalysis( %s, ''%s'', ''%s'', %s, %s );', inputname(1), inputname(1), input, direction, mat2str(peak), mat2str(peakWin) );
elseif nargin > 5
    EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin, varargin{:} );
    com = sprintf('%s = pop_tesa_peakanalysis( %s, ''%s'', ''%s'' ,%s, %s, %s );', inputname(1), inputname(1), input, direction, mat2str(peak), mat2str(peakWin), vararg2str(varargin) );
else
    if isempty(tepName)
        EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin, 'method', method, 'samples', samples );
        com = sprintf('%s = pop_tesa_peakanalysis( %s, ''%s'', ''%s'', %s, %s, ''method'', ''%s'', ''samples'', %s );', inputname(1), inputname(1), input, direction, mat2str(peak), mat2str(peakWin), method, mat2str(samples));
    else 
        EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin, 'method', method, 'samples', samples, 'tepName', tepName );
        com = sprintf('%s = pop_tesa_peakanalysis( %s, ''%s'', ''%s'', %s, %s, ''method'' ,''%s'', ''samples'', %s, ''tepName'', ''%s'' );', inputname(1), inputname(1), input, direction, mat2str(peak), mat2str(peakWin),  method,  mat2str(samples), tepName );    
    end
end

end
