% pop_tesa_detrend() - detrends the data by fiting and subtracting a function from
%                   each channel. Either a linear (fitted and subtratcted from 
%                   each trial), exponential or double exponential function 
%                   (fitted to average and subtracted from each trial) can be fitted. 
%                   Note that the Curve Fitting Toolbox is required
%                   to run either the exponential fit or the double
%                   exponential fit.
% 
% Usage:
%   >>  EEG = pop_tesa_detrend( EEG,); %pop up window
%   >>  EEG = pop_tesa_detrend( EEG, detrend, timeWin );
%
% Inputs (required):
%   EEG             - EEGLAB EEG structure
%   detrend         - string with type of detrend to perform; 'linear' |
%                   'exponential' | 'double'
%   timeWin         - vector with time range for detrending [t1,t2]
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples:
%   EEG = pop_tesa_detrend( EEG, 'linear', [11,500]);
%   EEG = pop_tesa_detrend( EEG, 'exponential', [11,500]);
%   EEG = pop_tesa_detrend( EEG, 'double', [11,500]);
% 
% See also:
%   pop_tesa_fastica, pop_tesa_edm, pop_tesa_pcasupress 

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

 function [EEG com] = pop_tesa_detrend( EEG, detrend, timeWin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
   
    geometry = {1 [1 1] 1 1 [1 1]};
    
    uilist = {{'style', 'text', 'string', 'Detrend data','fontweight','bold'} ...
              {'style', 'text', 'string', 'Type of detrend'} ...
              {'style', 'popupmenu', 'string', 'linear|exponential|double' 'tag' 'interp' } ...
              {}...
              {'style', 'text', 'string', 'Time range for detrend','fontweight','bold'} ...
              {'style', 'text', 'string', 'Start, end (ms)'} ...
              {'style', 'edit', 'string', '11, 500'}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Detrend data -- pop_tesa_detrend()', 'helpcom', 'pophelp(''pop_tesa_detrend'')');
    if isempty(result), return; end;
    
    %extract detrend type
    if result{1,1} == 1
        detrend = 'linear';
    elseif result{1,1} == 2 
        detrend = 'exponential'; 
    elseif result{1,1} == 3 
        detrend = 'double'; 
    end
    
    timeWin = str2num(result{1,2});
    
end

%remove data
EEG = tesa_detrend( EEG, detrend, timeWin );

% return the string command
com = sprintf('%s = pop_tesa_detrend( %s, ''%s'', %s );', inputname(1), inputname(1), detrend, mat2str(timeWin));

end
