% pop_tesa_peakanalysis() - finds peaks within a time window defined by user.
%                       For either ROI or GMFA analyses.
%                       Either positive or negative peaks are detected.
%                       Peaks are defined as data point which is
%                       larger/smaller than +/- x data points (default = 3,
%                       however this can also be defined by user). Results
%                       are saved in EEG structure under (either EEG.ROI or
%                       EEG.GMFA). Note that either pop_tesa_roianalysis or
%                       pop_tesa_gmfaanalysis must be run prior to this
%                       script.
%
% Usage:
%   >>  EEG = pop_tesa_peakanalysis( EEG ) - for pop up window
%   >>  EEG = pop_tesa_peakanalysis( EEG, input, direction, peak, peakWin )
%   >>  EEG = pop_tesa_peakanalysis( EEG, input, direction, peak, peakWin, varargin )
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   input           - string defining whether 'ROI' or 'GMFA' input is used
%   direction       - string defining whether peaks are 'positive' or
%                       'negative'. Use 'positive' for GMFA analysis.
%   peak            - vector defining the peak of interest. 
%                       Example: 25 (one peak), [25,60,..] (multiple peaks)
%   peakWin         - matrix defining the time windows to search for above
%                       peaks. Note that the number of peak windows defined should
%                       equal the number of peaks. 
%                       For example if peak = [25,60,180]; peakWin = [25,45;80,120;180,220]
% 
% Optional input pairs
%   'method','str'  - either 'largest'|'centre (default = largest). If multiple
%                       peaks are detected in a window, largest will search for
%                       the largest peak within the time window. Centre will 
%                       search for the peak closest to the latency defined
%                       in peak.
%   'samples',int   - int is an integer defining the number of samples
%                       either side of a peak that defines the peak.
%                       Peaks are defined as data point which is
%                       larger/smaller than +/- int data points (default = 3)
%   'roiName','str' - 'str' is a name to identify ROI analysis. This is
%                       useful if multiple different ROIs are to be analysed.
%                       The output will be stored under this name in EEG.ROI.
%                       Example: 'motor'
%                       Defaults are: R1,R2,R....
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

 function [EEG com] = pop_tesa_peakanalysis( EEG, input, direction, peak, peakWin, varargin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
      
    geometry = {1 [0.5 0.5 0.5] 1 [0.5 0.5 0.5] 1 [1 0.5] 1 [1 0.5] 1 1 1 1 [1 0.5] [1 0.5]};
    
    

    uilist = {{'style', 'text', 'string', 'Find peaks','fontweight','bold'} ...
              {'style', 'text', 'string', 'Input type for analysis'} ...
              {'style', 'popupmenu', 'string', 'ROI|GMFA' 'tag' 'input' } ...
              {}...
              {}...
              {'style', 'text', 'string', 'Direction of peaks'} ...
              {'style', 'popupmenu', 'string', 'positive|negative' 'tag' 'direction' } ...
              {}...
              {} ...
              {'style', 'text', 'string', 'Peaks to find (ms) [required]','fontweight','bold'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', '   Example: 25 (one peak), [25,60,...] (multiple peaks)','fontangle','italic'} ...
              {'style', 'text', 'string', 'Time window for peak search (ms) [required]','fontweight','bold'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', '   Example: [15,45] - for single peak between 15-45 ms.','fontangle','italic'} ...
              {'style', 'text', 'string', '   Example: [15,45; 46,75] - for multiple peaks between 15-45 ms (P25) and 46-75 ms (P60).','fontangle','italic'} ...
              {} ...
              {'style', 'text', 'string', 'Optional inputs','fontweight','bold'} ...
              {'style', 'text', 'string', 'Method for determining which peak to use (if multiple detected)'} ...
              {'style', 'popupmenu', 'string', 'largest|centre' 'tag' 'input' } ...
              {'style', 'text', 'string', 'Number of samples for defining peak'} ...
              {'style', 'edit', 'string', '3'}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Find peaks -- pop_tesa_peakanalysis()', 'helpcom', 'pophelp(''tesa_peakanalysis'')');
    
    %Check that peaks have been entered
    if strcmp(result{1,3},'')
        error('No peaks were defined. Please enter peaks.');
    end
    
    %Check that peak windows have been entered
    if strcmp(result{1,4},'')
        error('No peak windows were defined. Please enter peak windows.');
    end
    
    %Extract data for single pulse artifact find
    if result{1,1} == 1
        input = 'ROI';
    elseif result{1,1} == 2
        input = 'GMFA';
    end
    
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
    
    if strcmp(input,'ROI')
        if size(fieldnames(EEG.ROI),1) == 1
            temp = fieldnames(EEG.ROI);
            roiName = temp{1,1};
        elseif size(fieldnames(EEG.ROI),1) > 1
            temp1 = fieldnames(EEG.ROI);
            temp2 = strrep(temp1,'''','');
            temp3 = strrep(temp2,',','|');
            
            geometry = {[1 0.5]};
            
            uilist = {{'style', 'text', 'string', 'Select ROI for analysis'} ...
              {'style', 'popupmenu', 'string', temp3 'tag' 'input' }};

           result2 = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Select ROI');
           
           roiName = temp1{result2{1,1},1};
        end
    end
    
    if strcmp(input,'GMFA')
        temp = fieldname(EEG.GMFA);
        roiName = temp{1,1};
    end
    
end

%Run script from input
if nargin == 5
    EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin );
    com = sprintf('%s = pop_tesa_peakanalysis( %s, %s ,%s, %s, %s );', inputname(1), inputname(1), input, direction, num2str(peak), num2str(peakWin) );
elseif nargin > 5
    EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin, varargin{:} );
    com = sprintf('%s = pop_tesa_findpulse( %s, %s, %s ,%s, %s, %s );', inputname(1), inputname(1), input, direction, num2str(peak), num2str(peakWin), vararg2str(varargin) );
else
    EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin, 'method', method, 'samples', samples, 'roiName', roiName );
    com = sprintf('%s = pop_tesa_findpulse( %s, %s, %s, %s, %s, %s ,%s, %s, %s, %s, %s );', inputname(1), inputname(1), input, direction, result{1,3}, result{1,4}, 'method', method, 'sample', result{1,6}, 'roiName', roiName );
end

end
