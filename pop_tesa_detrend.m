% tesa_detrend() - detrends the data by fiting and subtracting either a linear, exponential
%                   or double exponential function to each channel and
%                   trial. Note that the Curve Fitting Toolbox is required
%                   to run either the exponential fit or the double
%                   exponential fit.
%
% Usage:
%   >>  EEG = tesa_detrend( EEG );
%   >>  EEG = tesa_detrend( EEG, detrend, timeWin );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   detrend         - string with type of detrend to perform; either
%                   'linear', 'exponential' or 'double'
%   tineWin         - vector with time range for detrending [t1,t2]
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

 function [EEG com] = pop_tesa_detrend( EEG, detrend, timeWin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
   
    geometry = {[1 1] 1 1 [1 0.5] [1 0.5]};
    
    uilist = {{'style', 'text', 'string', ['Type of detrend']} ...
              {'style', 'popupmenu', 'string', 'linear|exponential|double' 'tag' 'interp' } ...
              {}...
              {'style', 'text', 'string', 'Time range for detrend','fontweight','bold'} ...
              {'style', 'text', 'string', 'Start time (in ms)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'End time (in ms)'} ...
              {'style', 'edit', 'string', ''}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Detrend data -- pop_tesa_detrend()', 'helpcom', 'pophelp(''tesa_detrend'')');
    
    %extract detrend type
    if result{1,1} == 1
        detrend = 'linear';
    elseif result{1,1} == 2 
        detrend = 'exponential'; 
    elseif result{1,1} == 3 
        detrend = 'double'; 
    end
    
    if strcmp(result(1,2),'') || strcmp(result(1,3),'')
        error('Please enter the time range for detrend');  
    else
        timeWin = cellfun(@str2num,(result(:,2:3)));
    end 
    
end

%remove data
EEG = tesa_detrend( EEG, detrend, timeWin );

% return the string command
com = sprintf('%s = pop_tesa_detrend( %s, ''%s'', %s );', inputname(1), inputname(1), detrend, mat2str(timeWin));

end
