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

 function com = pop_tesa_plot( EEG, varargin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
    
    geometry = {1 [1 0.3] [1 0.3] 1 [1 0.3] 1 1 [1 0.3] 1 1 [1 0.3] [1 0.3] 1 [1 0.3] 1};

    uilist = {{'style', 'text', 'string', 'Plot TMS-evoked potentials','fontweight','bold'} ...
              {'style', 'text', 'string', 'Time (x-axis) limits in ms'} ...
              {'style', 'edit', 'string', '-100,300'} ...
              {'style', 'text', 'string', 'Amplitude (y-axis) limits in uV'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', '   Example: -30,30 [Leave empty for auto-scale]', 'fontangle','italic'} ...
              {'style', 'text', 'string', 'Type of data to plot'} ...
              {'style', 'popupmenu', 'string', 'data|ROI|GMFA' 'tag' 'data' } ...
              {} ...
              {'style', 'text', 'string', 'For EEG data','fontweight','bold'} ...
              {'style', 'text', 'string', 'Electrode to plot (e.g. Cz) [empty for all]'} ...
              {'style', 'edit', 'string', ''}...
              {} ...
              {'style', 'text', 'string', 'For analysed data (ROI, GMFA)','fontweight','bold'} ...
              {'style', 'text', 'string', 'Region of interest to plot'} ...
              {'style', 'edit', 'string', 'R1'}...
              {'style', 'text', 'string', 'Plot peak analysis results',} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'peaks' } ...
              {} ...
              {'style', 'text', 'string', 'Plot confidence interval',} ...
              {'Style', 'checkbox', 'string' 'on/off' 'value' 0 'tag' 'CI' }...
              {'style', 'text', 'string', '   Single electrode and ROI only','fontangle','italic'}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Plot TEPs -- pop_tesa_plot()', 'helpcom', 'pophelp(''pop_tesa_plot'')');
    
   
    %Extract data from result
    xLim = str2num(result{1,1});
    
    if strcmp(result{1,2},'')
        yLim = [];
    else
        yLim = str2num(result{1,2});
    end
    
    if result{1,3} == 1
        input = 'data';
    elseif result{1,3} == 2
        input = 'ROI';
    elseif result{1,3} == 3
        input = 'GMFA';
    end
    
    if strcmp(result{1,4},'')
        elec = [];
    else
        elec = result{1,4};
    end
    
    roiName = result{1,5};
    
    if result{1,6} == 0
        plotPeak = 'off';
    elseif result{1,7} == 1;
        plotPeak = 'on';
    end
    
    if result{1,7} == 0
        CI = 'off';
    elseif result{1,7} == 1
        CI = 'on';
    end   
   
end

%Run script from input
if nargin > 2
    tesa_plot(EEG,varargin{:});
    com = sprintf('%s = pop_tesa_plot( %s, %s, %s );', inputname(1), inputname(1), vararg2str(varargin) );
elseif nargin < 2
    tesa_plot(EEG,'xLim',xLim,'yLim',yLim,'input',input,'elec',elec,'roiName',roiName,'plotPeak',plotPeak,'CI',CI);
    com = sprintf('%s = pop_tesa_plot( %s, ''xLim'', %s,''yLim'', %s, ''input'', %s, ''elec'', %s, ''roiName'', %s, ''plotPeak'', %s, ''CI'', %s );', inputname(1), inputname(1), mat2str(xLim), mat2str(yLim), input, elec, roiName, plotPeak, CI );
end
    

end
