% pop_tesa_modifiedbandpassfilter() - do a variant of band-pass filtering, 
%                                 specialized to avoid propagating large-amplitude 
%                                 low-latency stimulation artifact to surrounding times
%                                 Auto-regressive prediction is used to pad
%                                 the time periods around where the TMS
%                                 artifact has been removed. Based on
%                                 function from the AARATEP pipeline by
%                                 Chris Cline. 
%
% Usage:
%   >>  EEG = pop_tesa_modifiedbandpassfilter( EEG );
%   >>  EEG = pop_tesa_modifiedbandpassfilter( EEG, 'key1', value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs:
%   'lowCutoff',num   - num defines the low band cut-off frequency for a high-pass filter. 
%                   num is a number in Hz. 
%                   default = 1
% 
%   'highCutoff',num   - num defines the high band cut-off frequency for a low-pass filter. 
%                   num is a number in Hz. 
%                   default = []
% 
%   'filterMethod','str'- 'str' is a filtering method. 'butterworth' | 'eegfiltnew' 
%                   default = 'butterworth'
%
%   'artifactTimespan', [vec] - vec is a vector defining the start and end
%                   time of TMS-pulse artifact removal in seconds. If tesa_removedata
%                   has been run, the most recent cut values will be used
%                   and multiplied by 3. 
%                   default = []
%
%  
%   'pieceWiseTimeToExtend',num - num defines the length of the time window for the 
%                   autoregressive prediction in seconds.
%                   default = 0.5
%
%   'filtOrder', num - defines the filter order (for Butterworth filters)
%                   default = 4
%
%   'filtType', 'str' - defines the filter type. 'auto' | 'bandpass'| 'bandstop' | 'highpass' | 'lowpass'
%                   auto = selects filter type based on lowCutoff and
%                   highCutoff inputs. If both are defined then a band-pass
%                   is used. If only one is defined then either high-pass
%                   or low-pass filter is used accordingly.
%                   default = 'auto'
%
%   Note that additional inputs are available using the
%   tesa_modifiedbandpassfilter function.
%
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples
%   EEG = pop_tesa_modifiedbandpassfilter( EEG ); % calls the pop up window
%   EEG = pop_tesa_modifiedbandpassfilter( EEG, 'lowCutoff', 2, 'pieceWiseTimeToExtend', 0.9 ); %user defined input
%
% See also:
%   pop_tesa_filtbutter

% This script was adapted by Nigel Rogasch for the TESA toolbox. Original
% code is available from:
% https://github.com/chriscline/AARATEPPipeline

% MIT License
% 
% Copyright (c) 2021 Chris Cline
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function [EEG, com] = pop_tesa_modifiedbandpassfilter( EEG, varargin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
        
    geometry = {1 1 [1 0.3] [1 0.3] [1 0.3] [1 0.3] 1 [1 0.3] [1 0.3] [1 0.3]};

    uilist = {{'style', 'text', 'string', 'Modifed band-pass filter with autoregressive padding','fontweight','bold'} ...
              {'style', 'text', 'string', 'Adapted from AARATEP toolbox'} ...
              {'style', 'text', 'string', 'High-pass (Hz)'} ...
              {'style', 'edit', 'string', '1'}...
              {'style', 'text', 'string', 'Low-pass (Hz)'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', 'Filter method'} ...
              {'style', 'popupmenu', 'string', 'butterworth|eegfiltnew', 'tag', 'interp' }...
              {'style', 'text', 'string', 'Vector defining the start and end time of TMS-pulse artifact removal in seconds*.'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', '* If tesa_removedata has been run, the most recent cut values will be used and multiplied by 3.'} ...
              {'style', 'text', 'string', 'Length of the time window for the autoregressive prediction in seconds'} ...
              {'style', 'edit', 'string', '0.5'} ...
              {'style', 'text', 'string', 'Filter order (Butterworth)'} ...
              {'style', 'edit', 'string', '4'}...
              {'style', 'text', 'string', 'Filter type (auto = selects based on high-pass/low-pass inputs)'} ...
              {'style', 'popupmenu', 'string', 'auto|bandpass|bandstop|highpass|lowpass', 'tag', 'interp' }...
              };
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Modified band-pass filter with autoregressive padding -- pop_tesa_modifiedbandpassfilter()', 'helpcom', 'pophelp(''pop_tesa_modifiedbandpassfilter'')');
    if isempty(result), return; end;
    
   
    %Extract data for single pulse artifact find
    if isempty(result{1})
        lowCutoff = [];
        lowCutoffString = '[]';
    else
        lowCutoff = str2num(result{1});
        lowCutoffString = result{1};
    end
    if isempty(result{2})
        highCutoff = [];
        highCutoffString = '[]';
    else
        highCutoff = str2num(result{2});
        highCutoffString = result{2};
    end
    if result{3} == 1
        filterMethod = 'butterworth';
    elseif result{3} == 2
        filterMethod = 'eegfiltnew';
    end
    artifactTimespan = str2num(result{4});
    pieceWiseTimeToExtend = str2num(result{5});
    filtOrder = str2num(result{6});
    if result{7} == 1
        filtType = 'auto';
    elseif result{7} == 2
        filtType = 'bandpass';
    elseif result{7} == 3
        filtType = 'bandstop';
    elseif result{7} == 4
        filtType = 'highpass';
    elseif result{7} == 5
        filtType = 'lowpass';
    end
    
end

%Run script from input
if nargin == 1
    if isempty(artifactTimespan)
        EEG = tesa_modifiedbandpassfilter(EEG,'lowCutoff',lowCutoff,'highCutoff',highCutoff,'filterMethod',filterMethod,'pieceWiseTimeToExtend',pieceWiseTimeToExtend,'filtOrder',filtOrder,'filtType',filtType);
        com = sprintf('%s = pop_tesa_modifiedbandpassfilter( %s, ''lowCutoff'', %s, ''highCutoff'', %s, ''filterMethod'', ''%s'', ''pieceWiseTimeToExtend'', %g, ''filtOrder'' , %g, ''filtType'', %g);', inputname(1), inputname(1), lowCutoffString, highCutoffString, filterMethod, pieceWiseTimeToExtend, filtOrder, filtType );
    else
        EEG = tesa_modifiedbandpassfilter(EEG,'lowCutoff',lowCutoff,'highCutoff',highCutoff,'filterMethod',filterMethod,'artifactTimespan',artifactTimespan,'pieceWiseTimeToExtend',pieceWiseTimeToExtend,'filtOrder',filtOrder,'filtType',filtType);
        com = sprintf('%s = pop_tesa_modifiedbandpassfilter( %s, ''lowCutoff'', %s, ''highCutoff'', %s, ''filterMethod'', ''%s'', ''artifactTimespan'', %g, ''pieceWiseTimeToExtend'', %g, ''filtOrder'' , %g, ''filtType'', %g );', inputname(1), inputname(1), lowCutoffString, highCutoffString, filterMethod, artifactTimespan, pieceWiseTimeToExtend, filtOrder, filtType );
    end
elseif nargin > 2
    EEG = tesa_modifiedbandpassfilter(EEG,varargin{:});
    com = sprintf('%s = pop_tesa_modifiedbandpassfilter( %s, %s );', inputname(1), inputname(1), vararg2str(varargin) );
end
    

end
