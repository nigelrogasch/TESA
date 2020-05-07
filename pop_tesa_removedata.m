% pop_tesa_removedata() - removes data between defined time points and replaces
%                   data with 0s or the average of a defined period (e.g. a baseline period).
%                   Can run on either continous or epoched data.
%                   Removed time points are stored in EEG.tmscut.
%
% Usage:
%   >>  EEG = pop_tesa_removedata( EEG ); % pop up window
%   >>  EEG = pop_tesa_removedata( EEG, cutTimesTMS ); %replace with 0s
%   >>  EEG = pop_tesa_removedata( EEG, cutTimesTMS, replaceTimes , cutEvent );% replace with average of defined period in ms around provided event/s
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   cutTimesTMS     - (required) vector with time range for removing TMS artifact in ms. [t1,t2]
%                       Note that the times are relevant to events, not the epoch. For example, if an event is located at -100
%                       ms in an epoch, running the algorithm with [-10,10] would remove between -110 and -90 ms.
%                       Example: [-10,10]
%   replaceTimes    - (optional) vector with time range for calculating average to replace removed data in ms. [t1,t2]
%                       If not included, data will be replaced with 0s.
%                       As above, times are relative to events, not the epoch.
%                       Example: [-500, -100]
%   cutEvent        - (optional) cell with strings indicating event/s to remove data around. {'string'}. 
%                       If left empty, tesa_removedata will remove data for all events.
%                       Example: {'TMS'}
%                       Example: {'single','paired'}
%                       Example: {'1'}
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = pop_tesa_removedata( EEG ); % pop up window
%   EEG = pop_tesa_removedata( EEG, [-10,10] ); %replace with 0s
%   EEG = pop_tesa_removedata( EEG, [-10,10], [-500,-100], {'TMS'} ); %replace with average of defined period in ms around event 'TMS'
%   EEG = pop_tesa_removedata( EEG, [-10,10], [], {'TMS'} ); %replace with 0s around event 'TMS'
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

% Change log:
% 20.6.2018: Changed function to work on both continuous and epoched data


 function [EEG com] = pop_tesa_removedata( EEG, cutTimesTMS, replaceTimes, cutEvent )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

if nargin < 1
    error('Not enough inputs. Please include the EEGLAB EEG structure.');
end

if nargin == 1
    cutTimesTMS = [];
	replaceTimes = [];
    cutEvent = [];
end

if nargin == 2
	replaceTimes = [];
    cutEvent = [];
end

if nargin == 3
    cutEvent = [];
end

% Generate cutString input
if ~isempty(cutEvent) && exist('cutString','var') ~= 1
    cutString = '';
    for inx = 1:length(cutEvent)
        cutString = [cutString,' ''',cutEvent{inx},''' '];
    end
end

% pop up window
% -------------
if nargin < 2
   
   cbevent = ['if ~isfield(EEG.event, ''type'')' ...
               '   errordlg2(''No type field'');' ...
               'else' ...
               '   tmpevent = EEG.event;' ...
               '   if isnumeric(EEG.event(1).type),' ...
               '        [tmps,tmpstr] = pop_chansel(unique([ tmpevent.type ]));' ...
               '   else,' ...
               '        [tmps,tmpstr] = pop_chansel(unique({ tmpevent.type }));' ...
               '   end;' ...
               '   if ~isempty(tmps)' ...
               '       set(findobj(''parent'', gcbf, ''tag'', ''events''), ''string'', tmpstr);' ...
               '   end;' ...
               'end;' ...
               'clear tmps tmpevent tmpv tmpstr tmpfieldnames;' ];
    
    geometry = {1 [1 0.5] 1 1 [1 0.5] 1 1 1 [2 2 0.5]};

    uilist = {{'style', 'text', 'string', 'Remove TMS artifact','fontweight','bold'} ...
              {'style', 'text', 'string', 'Time to remove TMS artifact (ms): start, end'} ...
              {'style', 'edit', 'string', '-2, 10'} ...
              {} ...
              {'style', 'text', 'string', 'Replace data with baseline average','fontweight','bold'} ...
              {'style', 'text', 'string', 'Baseline period for averaging (ms): start, end'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'Optional (leave blank to replace with 0s)','fontAngle','italic'}...
              {} ...
              {'style', 'text', 'string', 'Replace data around specific events','fontweight','bold'} ...
              {'style', 'text', 'string', 'Choose events (optional)' } ...
              {'style', 'edit', 'string', '', 'tag', 'events' } ...
              {'style', 'pushbutton', 'string', '...', 'callback', cbevent }};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Remove artifact data -- pop_tesa_removedata()', 'helpcom', 'pophelp(''pop_tesa_removedata'')');
    if isempty(result), return; end;
    
    %extract TMS artifact removal times.
    cutTimesTMS = str2num(result{1,1});
    replaceTimes = str2num(result{1,2});
    if strcmp(result{1,3},'') || isempty(result{1,3})
        cutEvent = [];
    else
        cutEvent = strtrim(strsplit(result{1,3},' '));
        cutString = ['''',result{1,3},''''];
        cutString = strrep(cutString,' ',''',''');
        cutString = strrep(cutString,' ','');
    end
    
end
 
%remove data
EEG = tesa_removedata( EEG, cutTimesTMS, replaceTimes, cutEvent );

% return the string command
if isempty(replaceTimes) && isempty(cutEvent)
    com = sprintf('%s = pop_tesa_removedata( %s, %s );', inputname(1), inputname(1), mat2str(cutTimesTMS));
elseif ~isempty(replaceTimes) && isempty(cutEvent)
    com = sprintf('%s = pop_tesa_removedata( %s, %s, %s );', inputname(1), inputname(1), mat2str(cutTimesTMS), mat2str(replaceTimes));
elseif isempty(replaceTimes) && ~isempty(cutEvent)
    com = sprintf('%s = pop_tesa_removedata( %s, %s, [], {%s} );', inputname(1), inputname(1), mat2str(cutTimesTMS), cutString);
elseif ~isempty(replaceTimes) && ~isempty(cutEvent)
    com = sprintf('%s = pop_tesa_removedata( %s, %s, %s, {%s} );', inputname(1), inputname(1), mat2str(cutTimesTMS), mat2str(replaceTimes), cutString);
end

end
