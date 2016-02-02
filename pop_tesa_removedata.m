% pop_tesa_removedata() - removes data between defined timepoints and replaces
%                       data with 0s. Stores removed timepoints.If time range 
%                       is not given, a window pops up
%                       to ask for the value of the additional 
%                       parameters.   
%
% Usage:
%   >>  EEG = pop_tesa_removedata( EEG ); %pop-up window mode
%   >>  EEG = pop_tesa_removedata( EEG, cutTimesTMS, cutTimesRec );
%
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

 function [EEG com] = pop_tesa_removedata( EEG, cutTimesTMS, cutTimesRec )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

if nargin == 2;
	cutTimesRec = [];
end

% pop up window
% -------------
if nargin < 2
   
    geometry = {1 [1 0.5] [1 0.5] 1 1 [1 0.5] [1 0.5]};

    uilist = {{'style', 'text', 'string', 'Remove TMS artifact','fontweight','bold'} ...
              {'style', 'text', 'string', 'Minimum time to remove TMS artifact (in ms)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'Maximum time to remove TMS artifact (in ms)'} ...
              {'style', 'edit', 'string', ''}...
              {} ...
              {'style', 'text', 'string', 'Remove recharge artifact (optional)','fontweight','bold'} ...
              {'style', 'text', 'string', 'Minimum time to remove recharge artifact (in ms)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'Maximum time to remove recharge artifact (in ms)'} ...
              {'style', 'edit', 'string', ''}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Remove artifact data -- pop_tesa_removedata()', 'helpcom', 'pophelp(''tesa_removedata'')');
    
    %extract TMS artifact removal times. Error if not given.
    if isempty(result)
        error('No times for data removal were entered. Script terminated')     
    else 
        cutTimesTMS = cellfun(@str2num,(result(:,1:2)));     
    end
    
    %extract recharge artifact removal times (optional).
    if strcmp(result(1,3),'')
        cutTimesRec = [];
    elseif ~strcmp(result(1,3),'')
        cutTimesRec = cellfun(@str2num,(result(:,3:4)));
    end
    
end

%remove data
EEG = tesa_removedata( EEG, cutTimesTMS, cutTimesRec );

% return the string command
if isempty(cutTimesRec)
    com = sprintf('%s = pop_tesa_removedata( %s, %s );', inputname(1), inputname(1), mat2str(cutTimesTMS));
elseif ~isempty(cutTimesRec)
    com = sprintf('%s = pop_tesa_removedata( %s, %s, %s );', inputname(1), inputname(1), mat2str(cutTimesTMS), mat2str(cutTimesRec));
end

end
