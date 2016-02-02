% pop_tesa_interpdata()  - replaces removed data using interpolated data. If
%                           interpolation type is not given, a window pops up asking which type of
%                           interpolation to use. Note that either tesa_removedata or
%                           pop_tesa_removedata must be ran first.
%
% Usage:
%   >>  EEG = tesa_interpdata( EEG, interpolation, interpWin );
%
% Inputs:
%   EEG                 - EEGLAB EEG structure
%   interpolatation     - string describing type of interpolation, either
%                           'linear' or 'cubic'
%   interpWin           - (optional) vector with times before and after
%                           artefact window for fitting cubic function. 
%                           default = [20,20];
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

 function [EEG com] = pop_tesa_interpdata( EEG, interpolation, interpWin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

if nargin == 2
	interpWin = [];
end

% pop up window
% -------------
if nargin < 2
   
    geometry = {[1 1] 1 1 [1 0.5] [1 0.5]};
    
    uilist = {{'style', 'text', 'string', 'Type of interpolation'} ...
              {'style', 'popupmenu', 'string', 'Linear|Cubic' 'tag' 'interp' } ...
              {}...
              {'style', 'text', 'string', 'Time window for fitting cubic (optional)','fontweight','bold'} ...
              {'style', 'text', 'string', 'Time before removed data (in ms) [default=20]'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'Time after removed data (in ms) [default=20]'} ...
              {'style', 'edit', 'string', ''}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Interpolate removed data -- pop_tesa_interpolatedata()', 'helpcom', 'pophelp(''tesa_interpolatedata'')');
    
    %extract interpolation type
    if result{1,1} == 1
        interpolation = 'linear';
        interpWin = [];
    elseif result{1,1} == 2 
        interpolation = 'cubic';    
    end
    
    %extract interpolation window for cubic (optional).
    if strcmp(interpolation,'cubic')
        if strcmp(result(1,2),'') && strcmp(result(1,3),'')
            interpWin = [];
        elseif ~strcmp(result(1,2),'') && ~strcmp(result(1,3),'')
            interpWin = cellfun(@str2num,(result(:,2:3)));
        else
            error('Please either enter both a before and after window time or leave both fields blank');
        end
    end
end

%remove data
EEG = tesa_interpdata( EEG, interpolation, interpWin );

% return the string command
if isempty(interpWin)
    com = sprintf('%s = pop_tesa_interpdata( %s, ''%s'' );', inputname(1), inputname(1), interpolation);
elseif ~isempty(interpWin)
    com = sprintf('%s = pop_tesa_interpdata( %s, ''%s'', %s );', inputname(1), inputname(1), interpolation, mat2str(interpWin));
end

end
