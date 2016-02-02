% pop_tesa_removedata() - removes data between defined timepoints and replaces
%                       data with NaNs. Stores removed timepoints.If time range 
%                       is not given, a window pops up
%                       to ask for the value of the additional 
%                       parameters.   
%
% Usage:
%   >>  EEG = pop_tesa_removedata( EEG ); %pop-up window mode
%   >>  EEG = pop_tesa_removedata( EEG, cutTimesTMS, cutTimesRec );

%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   cutTimes        - vector with time range for removing data [t1,t2]
%   cutTimesRec     - (optional) vector with time range for removing recharge artifact [t1,t2]
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

 function [EEG com] = pop_tesa_medianfilt( EEG, timeWin, filtOrd )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

if nargin == 2
	filtOrd = [];
end

% pop up window
% -------------
if nargin < 2
   
    geometry = {1 [1 0.5] [1 0.5] 1 1 [1 0.5]};

    uilist = {{'style', 'text', 'string', 'Time window for median filter','fontweight','bold'} ...
              {'style', 'text', 'string', 'Minimum time for filter (in ms)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'Maximum time for filter (in ms)'} ...
              {'style', 'edit', 'string', ''}...
              {} ...
              {'style', 'text', 'string', 'Filter order (optional)','fontweight','bold'} ...
              {'style', 'text', 'string', 'Order of filter (default = 29)'} ...
              {'style', 'edit', 'string', ''}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Median filter -- pop_tesa_medianfilt()', 'helpcom', 'pophelp(''tesa_medianfilt'')');
    
    %extract TMS artifact removal times. Error if not given.
    if isempty(result)
        error('No times for median filter were entered. Script terminated')     
    else 
        timeWin = cellfun(@str2num,(result(:,1:2)));     
    end
    
    %extract recharge artifact removal times (optional).
    if strcmp(result(1,3),'')
        filtOrd = [];
    elseif ~strcmp(result(1,3),'')
        filtOrd = cellfun(@str2num,(result(:,3)));
    end
    
end

%remove data
EEG = tesa_medianfilt( EEG, timeWin, filtOrd );

% return the string command
if isempty(filtOrd)
    com = sprintf('%s = pop_tesa_medianfilt( %s, %s );', inputname(1), inputname(1), mat2str(timeWin));
elseif ~isempty(filtOrd)
    com = sprintf('%s = pop_tesa_medianfilt( %s, %s, %s );', inputname(1), inputname(1), mat2str(timeWin), mat2str(filtOrd));
end

end
