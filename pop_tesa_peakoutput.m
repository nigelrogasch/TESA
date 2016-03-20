% pop_tesa_peakoutput() - returns the results for the peak analysis in a table
%                       in the workspace. This can be calculated on either the
%                       peak latencies determined using 'tesa_peakanalysis',
%                       or on fixed latencies provided by the user. Users can 
%                       also opt to have the average amplitude incorporating data 
%                       points either side of the peak instead of the absolute
%                       peak amplitude. Either the average of the TEP curve or 
%                       the area under the curve can be calculated (GMFA only). 
%                      
% Usage:
%   >>  output = pop_tesa_peakoutput( EEG ); %pop up window
%   >>  output = pop_tesa_peakoutput( EEG, 'key1', value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs
%   'winType', 'str'  - 'individual' | 'fixed'. Calculates values using either
%                       the latencies determined for each individual
%                       participant (using tesa_peakanalysis) or set latencies 
%                       provided by the user in 'fixedPeak' (see below).
%                       Default = 'individual'  
%   'calcType', 'str' - 'amplitude' | 'area'. Indicates whether to return the
%                       average amplitude of the TEP time series or the area under
%                       the curve (area under curve only for GMFA analysis).
%                       For area under curve, an analysis time window must 
%                       also be entered using averageWin.
%                       Default = 'amplitude'
%   'tepName', 'str'  - string indicating which specific ROI or GMFA TEP to
%                       give output for. If not indicated, output for all
%                       TEP fields will be given.
%                       Examples: 'R1', 'motor' etc.
%   'averageWin', int - Integer describing a time window +/- the peak (in ms)
%                       in which an average amplitude/area will be taken. 
%                       If left empty, the absolute amplitude at the peak latency 
%                       will be returned. A value is required for
%                       calculating area under the curve for GMFA.
%                       Example: 5
%   'fixedPeak', [int]- required for 'winType' fixed. Integer or vector describing 
%                       fixed latencies for calculating average or area under 
%                       the curve.
%                       Examples: [30], [30, 60, 100] etc.
%   'tablePlot','str' - 'on' | 'off'. Plots a table with results from peak
%                       analysis.
%     
% Outputs:
%   output             - table with results from peak analysis
% 
% Examples:
%   output = pop_tesa_peakoutput( EEG ); %returns amplitude at individual latencies for all peaks defined with tesa_peakanalysis
%   output = pop_tesa_peakoutput( EEG, 'averageWin', 5 ); %returns amplitude averaged +/- 5 ms from individual latencies for all peaks defined with tesa_peakanalysis
%   output = pop_tesa_peakoutput( EEG, 'averageWin', 10, 'calcType', 'area', 'tepName','parietal' ); %returns area under the curve +/- 10 ms from individual latencies for all peaks defined in parietal region of interest using tesa_peak analysis
%   output = pop_tesa_peakoutput( EEG, 'winType', 'fixed', 'fixedPeak', [30,60,100,180],'averageWin', 5 ); %returns amplitude averaged +/- 5 ms from peaks given in 'fixedPeak'. It is not necessary to run tesa_peakanalysis for this option
% 
% See also:
%   pop_tesa_tepextract, pop_tesa_peakanalysis 

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


 function [output, com] = pop_tesa_peakoutput( EEG, varargin )

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
    
    tepInput = [];
    if strcmpi(roi,'yes') && strcmpi(gmfa,'yes')
        tepInput = ['All'; roiFields; gmfaFields];
    elseif strcmpi(roi,'yes') && strcmpi(gmfa,'no')
        tepInput = ['All'; roiFields];
    elseif strcmpi(roi,'no') && strcmpi(gmfa,'yes')
        tepInput = ['All'; gmfaFields];
    end
    
    if isempty(tepInput)
        error('There are no ROI or GMFA analyses in this data. Please run tesa_tep_extract first.');
    end
    
    tepAll = tepInput;
    
    for a = 1:size(tepInput,1)-1
        tepInput{a,1} = [tepInput{a,1},'|'];
    end
    
    tepIn = [tepInput{:}];
    
    geometry = {[1 0.5] [1 0.5] 1 [1 0.5] 1 [1 0.5] 1 [1 0.5] 1 1 1 [1 0.5] 1 1};

    uilist = {{'style', 'text', 'string', 'Output peak analysis results','fontweight','bold'} ...
              {'Style', 'checkbox', 'string' 'Table plot: on/off' 'value' 1 'tag' 'pair' } ...
              {'style', 'text', 'string', 'TEP for output'} ...
              {'style', 'popupmenu', 'string', tepIn, 'tag', 'tep' } ...
              {}...
              {'style', 'text', 'string', 'Calculation type (area for GMFA only)'} ...
              {'style', 'popupmenu', 'string', 'amplitude|area', 'tag', 'calcType' } ...
              {}...
              {'style', 'text', 'string', 'Peaks for calculation'} ...
              {'style', 'popupmenu', 'string', 'individual|fixed', 'tag', 'winType' } ...
              {}...                  
              {'style', 'text', 'string', 'Window for averaging over peak (+/- ms) [optional]','fontweight','bold'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', '   Instead of returning the absolute amplitude, entering a value here will return','fontangle','italic'}...
              {'style', 'text', 'string', '   the average value of the window around the peak. e.g. 5 ms either side of peak','fontangle','italic'}...
              {}...
              {'style', 'text', 'string', 'Fixed latencies for calculating peaks (ms)','fontweight','bold'} ...
              {'style', 'edit', 'string', ''}... 
              {'style', 'text', 'string', '   Latencies of peaks for fixed analysis [required for fixed peak calculation]. ','fontangle','italic'}...
              {'style', 'text', 'string', '   Example: 30 - for single peak. 30, 60, 100, 180 - for multiple peaks.','fontangle','italic'}...              
              };
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Output peak -- pop_tesa_peakoutput()', 'helpcom', 'pophelp(''pop_tesa_peakoutput'')');
    output = [];
    if isempty(result), return; end;  
    
    %Run script
    if result{1,1} == 1
        tablePlot = 'on';
    else
        tablePlot = 'off';
    end
    
    tepFor = tepAll{result{1,2},1};
    if ~strcmp(tepFor,'All')
        tepFor = strtrim(strsplit(tepFor,' '));
        tepName = tepFor{1,2};
    else
        tepName = [];
    end
    
    if result{1,3} == 1
        calcType = 'amplitude';
    elseif result{1,3} == 2
        calcType = 'area';
    end
    
    if result{1,4} == 1
        winType = 'individual';
    elseif result{1,4} == 2
        winType = 'fixed';
    end
    
    if strcmp(result{1,5},'')
        averageWin = [];
    else
        averageWin = str2num(result{1,5});
    end
    
    if ~strcmp(result{1,6},'')
        fixedPeak = str2num(result{1,6});
    else
        fixedPeak = [];
    end
    
    if strcmp(calcType,'area') && isempty(averageWin)
        error('A window for averaging around the peak must be entered for calculation type ''area''.');
    end
    
    if strcmp(winType,'fixed') &&  isempty(fixedPeak)
        error('Fixed peak latencies must be included for calculation on fixed peaks.');
    end
    
    if strcmp(winType,'individual') && ~isempty(fixedPeak)
        warning('Calculation type ''individual'' was selected, however values were also entered in fixed peak latencies. Returned values are calculated on individual peak latencies.');
    end
       
end

%Run script from input
if nargin > 1;
    output = tesa_peakoutput(EEG,varargin{:});
    com = sprintf('%s = pop_tesa_peakoutput( %s, %s );', inputname(1), inputname(1), vararg2str(varargin) );
end

%perform roi analysis and return the string command
if nargin == 1;
    if isempty(tepName)
        output = tesa_peakoutput(EEG,'tepName',tepName,'calcType',calcType,'winType',winType,'averageWin',averageWin,'fixedPeak',fixedPeak,'tablePlot',tablePlot);
        com = sprintf('output = pop_tesa_peakoutput( %s, ''tepName'', %s, ''calcType'', ''%s'', ''winType'', ''%s'', ''averageWin'', %s, ''fixedPeak'', %s, ''tablePlot'', ''%s'' );', inputname(1), mat2str(tepName), calcType, winType, mat2str(averageWin), mat2str(fixedPeak), tablePlot);
    else
        output = tesa_peakoutput(EEG,'tepName',tepName,'calcType',calcType,'winType',winType,'averageWin',averageWin,'fixedPeak',fixedPeak,'tablePlot',tablePlot);
        com = sprintf('output = pop_tesa_peakoutput( %s, ''tepName'', ''%s'', ''calcType'', ''%s'', ''winType'', ''%s'', ''averageWin'', %s, ''fixedPeak'', %s, ''tablePlot'', ''%s'' );', inputname(1), tepName, calcType, winType, mat2str(averageWin), mat2str(fixedPeak), tablePlot);
    end
end

end
