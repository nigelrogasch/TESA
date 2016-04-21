% pop_tesa_filtmedian() - applies a 1-d median filter of nth-order to remove
%                       artifacts such as spikes and muscle artifacts. 
%                       Note that spike artifacts can be selected using the 
%                       gui option in tesa_findpulsepeak.
%
% Usage:
%   >>  EEG = pop_tesa_filtmedian( EEG ); %Pop up window
%   >>  EEG = pop_tesa_filtmedian( EEG, timeWin, filtOrd, label ); %Default use
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   timeWin         - Vector with time range for applying median filter in ms [t1,t2] (example: [-1,20])
%                       Note that t1 must be 0 or negative and t2 positive
%   filtOrd         - Integer indicating filter order (example: 30). The
%                      filter order determines the number of samples considered 
%                      either side of the data point when calculating the median value. 
%   label           - String indicating which event to filter around (example: 'TMS' or 'spike')
% 
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples:
%   EEG = pop_tesa_filtmedian( EEG, [-2,2], 3, 'spike' ); %Apply a 3rd order median filter to remove small spike artifacts
%
% See also:
%   pop_tesa_filtbutter

% Copyright (C) 2016  Nigel Rogasch, Monash University,
% nigel.rogasch@monash.edu
%
% Authors: Nathan Rose, Nigel Rogasch
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

 function [EEG com] = pop_tesa_filtmedian( EEG, timeWin, filtOrd, label )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

if nargin > 1
    if nargin ~= 4
        error('Not enough input arguments.');
    end
end

% Create event list
for x = 1:size(EEG.event,2)
    if ~ischar(EEG.event(x).type)
        EEG.event(x).type = num2str(EEG.event(x).type);
    end
    trigAll{x,1} = EEG.event(x).type;
end
trigUnique = unique(trigAll);
trigOr = trigUnique;
for x = 1:size(trigOr,1)-1
    trigOr{x,1} = [trigOr{x,1}, '|'];
end
trigIn = [trigOr{:}];

% Find default event ('TMS')
if sum(strcmpi('TMS',trigUnique)) > 0
    trigVal = find(strcmpi('TMS',trigUnique));
else
    trigVal = 1;
end

% pop up window
% -------------
if nargin < 2
   
    geometry = {1 [1 0.5] [1 0.5] [1 0.5]};

    uilist = {{'style', 'text', 'string', 'Median filter','fontweight','bold'} ...
              {'style', 'text', 'string', 'Time window (required): t1, t2'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', 'Filter order (required)'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', 'Event to be filtered.'} ...
              {'style', 'popupmenu', 'string', trigIn, 'tag', 'interp', 'Value', trigVal }};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Median filter -- pop_tesa_filtmedian()', 'helpcom', 'pophelp(''pop_tesa_filtmedian'')');
    if isempty(result), return; end;
    
    timeWin = str2num(result{1,1});
    filtOrd = str2num(result{1,2});
    label = trigUnique{result{1,3},1};
    
end

%remove data
EEG = tesa_filtmedian( EEG, timeWin, filtOrd, label );

% return the string command
com = sprintf('%s = pop_tesa_medianfilt( %s, %s, %s, ''%s'' );', inputname(1), inputname(1), mat2str(timeWin), mat2str(filtOrd), label);

end
