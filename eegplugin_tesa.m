% eegplugin_tesa() - EEGLAB plugin for analysing TMS-EEG data
%
% Usage:
%   >> eegplugin_tesa(fig, try_strings, catch_strings);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   try_string    - [struct] "try" strings for menu callbacks.
%   catch_string  - [struct] "catch" strings for menu callbacks.
%

% Copyright (C) 2015 Nigel Rogasch, Monash University,
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

function vers = eegplugin_tesa(fig, try_strings, catch_strings)

    vers = 'tesa1.1';
    if nargin < 3
        error('eegplugin_tesa requires 3 arguments');
    end

    % find tools menu
    % ---------------
    menu = findobj(fig, 'tag', 'tools'); 

    % menu callbacks
    % --------------
    comremovedata = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_removedata(EEG);' catch_strings.new_and_hist];
    comsortcomps = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_sortcomps(EEG);' catch_strings.new_and_hist];
    compcacompress = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_pcacompress(EEG);' catch_strings.new_and_hist];
    cominterpdata = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_interpdata(EEG);' catch_strings.new_and_hist];

    % create menus
    % -------------------------
    submenu = uimenu( menu, 'Label', 'TMS-EEG signal analyser (TESA)', 'separator', 'on');
    uimenu( submenu, 'Label', 'Remove artifact data'  , 'CallBack', comremovedata);
    uimenu( submenu, 'Label', 'Sort components'  , 'CallBack', comsortcomps, 'separator', 'on');
    uimenu( submenu, 'Label', 'PCA compress'  , 'CallBack', compcacompress);
    uimenu( submenu, 'Label', 'Interpolate missing data'  , 'CallBack', cominterpdata, 'separator', 'on');