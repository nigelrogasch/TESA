% pop_tesa_interpdata()     - replaces removed data using interpolated data.
%                           Note that either tesa_removedata or
%                           pop_tesa_removedata must be ran prior to this function.
% Usage:
%   >>  EEG = pop_tesa_interpdata( EEG ); % Pop up window
%   >>  EEG = pop_tesa_interpdata( EEG, interpolation );
%   >>  EEG = pop_tesa_interpdata( EEG, interpolation, interpWin );
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
%   EEG                 - EEGLAB EEG structure
% 
% Examples
%   EEG = pop_tesa_interpdata( EEG, 'linear' ); %replaces missing data with linear interpolation. Linear function is fitted on data point before and after missing data.
%   EEG = pop_tesa_interpdata( EEG, 'cubic', [50,50] ); %replaces mising data with cubic interpolation. Cubic is fitted on data 50 ms before and 50 ms after missing data
%
% See also:
%   pop_tesa_removedata

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
   
    geometry = {1 [1 0.5] 1 1 [1 0.3] [1 0.3]};
    
    uilist = {{'style', 'text', 'string', 'Interpolate missing data','fontweight','bold'} ...
              {'style', 'text', 'string', 'Type of interpolation'} ...
              {'style', 'popupmenu', 'string', 'Linear|Cubic' 'tag' 'interp' } ...
              {}...
              {'style', 'text', 'string', 'Time window for fitting cubic (optional)','fontweight','bold'} ...
              {'style', 'text', 'string', 'Time before removed data (in ms)'} ...
              {'style', 'edit', 'string', '20'} ...
              {'style', 'text', 'string', 'Time after removed data (in ms)'} ...
              {'style', 'edit', 'string', '20'}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Interpolate removed data -- pop_tesa_interpdata()', 'helpcom', 'pophelp(''pop_tesa_interpdata'')');
    if isempty(result), return; end;
    
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
