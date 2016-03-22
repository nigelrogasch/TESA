% pop_tesa_removedata() - removes data between defined time points and replaces
%                   data with 0s or the average of a defined period (e.g. a baseline period). 
%                   Removed time points are stored in EEG.tmscut.
%
% Usage:
%   >>  EEG = pop_tesa_removedata( EEG ); % pop up window
%   >>  EEG = pop_tesa_removedata( EEG, cutTimesTMS ); %replace with 0s
%   >>  EEG = pop_tesa_removedata( EEG, cutTimesTMS, replaceTimes );% replace with average of defined period in ms
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   cutTimesTMS     - (required) vector with time range for removing TMS artifact in ms. [t1,t2]
%                       Example: [-10,10]
%   replaceTimes    - (optional) vector with time range for calculating average to replace removed data in ms. [t1,t2]
%                      If not included, data will be replaced with 0s.
%                       Example: [-500, -100]
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = pop_tesa_removedata( EEG, [-10,10] ); %replace with 0s
%   EEG = pop_tesa_removedata( EEG, [-10,10], [-500,-100] ); %replace with average of defined period in ms
% 
% See also:
%   pop_tesa_interpdata 

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

 function [EEG com] = pop_tesa_removedata( EEG, cutTimesTMS, replaceTimes )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

if nargin == 2;
	replaceTimes = [];
end

% pop up window
% -------------
if nargin < 2
   
    geometry = {1 [1 0.5] 1 1 [1 0.5] 1};

    uilist = {{'style', 'text', 'string', 'Remove TMS artifact','fontweight','bold'} ...
              {'style', 'text', 'string', 'Time to remove TMS artifact (ms): start, end'} ...
              {'style', 'edit', 'string', '-2, 10'} ...
              {} ...
              {'style', 'text', 'string', 'Replace data with baseline average','fontweight','bold'} ...
              {'style', 'text', 'string', 'Baseline period for averaging (ms): start, end'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'Optional (leave blank to replace with 0s)','fontAngle','italic'}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Remove artifact data -- pop_tesa_removedata()', 'helpcom', 'pophelp(''pop_tesa_removedata'')');
    if isempty(result), return; end;
    
    %extract TMS artifact removal times.
    cutTimesTMS = str2num(result{1,1});
    replaceTimes = str2num(result{1,2});
    
end

%remove data
EEG = tesa_removedata( EEG, cutTimesTMS, replaceTimes );

% return the string command
if isempty(replaceTimes)
    com = sprintf('%s = pop_tesa_removedata( %s, %s );', inputname(1), inputname(1), mat2str(cutTimesTMS));
elseif ~isempty(replaceTimes)
    com = sprintf('%s = pop_tesa_removedata( %s, %s, %s );', inputname(1), inputname(1), mat2str(cutTimesTMS), mat2str(replaceTimes));
end

end
