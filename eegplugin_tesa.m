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

% Copyright (C) 2016 Nigel Rogasch, Monash University,
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

    vers = 'tesa1.1.1';
    if nargin < 3
        error('eegplugin_tesa requires 3 arguments');
    end

    % find tools menu
    % ---------------
    menu = findobj(fig, 'tag', 'tools'); 

    % menu callbacks
    % --------------
    comfindpulse = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_findpulse(EEG);' catch_strings.add_to_hist];
    comfindpulsealt = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_findpulsepeak(EEG);' catch_strings.add_to_hist];
    comfixevent = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_fixevent(EEG);' catch_strings.add_to_hist];
    comremovedata = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_removedata(EEG);' catch_strings.new_and_hist];
    cominterpdata = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_interpdata(EEG);' catch_strings.new_and_hist];
    comfastica = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_fastica(EEG);' catch_strings.new_and_hist];
    comcompselect = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_compselect(EEG);' catch_strings.new_and_hist];
    comcompplot = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_compplot(EEG);' catch_strings.new_and_hist];
    comedm = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_edm(EEG);' catch_strings.new_and_hist];
    compcacompress = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_pcacompress(EEG);' catch_strings.new_and_hist];
    compcasuppress = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_pcasuppress(EEG);' catch_strings.new_and_hist];
    comdetrend = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_detrend(EEG);' catch_strings.new_and_hist];
    comSOUND = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_sound(EEG);' catch_strings.new_and_hist];
    comSSPSIR = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_sspsir(EEG);' catch_strings.new_and_hist];
    comfiltbutter = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_filtbutter(EEG);' catch_strings.new_and_hist];
    comfiltmedian = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_filtmedian(EEG);' catch_strings.new_and_hist];
    comtepextract = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_tepextract(EEG);' catch_strings.add_to_hist];
    compeakanalysis = [try_strings.no_check '[EEG LASTCOM] = pop_tesa_peakanalysis(EEG);' catch_strings.add_to_hist];
    compeakoutput = [try_strings.no_check '[output LASTCOM] = pop_tesa_peakoutput(EEG);' catch_strings.add_to_hist];
    complot = [try_strings.no_check '[LASTCOM] = pop_tesa_plot(EEG);' catch_strings.add_to_hist];
    compeakoutputgroup = [try_strings.no_check '[output LASTCOM] = pop_tesa_peakoutputgroup(EEG);' catch_strings.add_to_hist];
    complotgroup = [try_strings.no_check '[LASTCOM] = pop_tesa_plotgroup(EEG);' catch_strings.add_to_hist];

    % create menus
    % -------------------------
    submenu = uimenu( menu, 'Label', 'TMS-EEG signal analyser (TESA)', 'separator', 'on');
    uimenu( submenu, 'Label', 'Find TMS pulse'  , 'CallBack', comfindpulse);
    uimenu( submenu, 'Label', 'Find TMS pulse (alternative)'  , 'CallBack', comfindpulsealt);
    uimenu( submenu, 'Label', 'Fix TMS pulse latency'  , 'CallBack', comfixevent);
    uimenu( submenu, 'Label', 'Remove artifact data'  , 'CallBack', comremovedata, 'separator', 'on');
    uimenu( submenu, 'Label', 'Interpolate removed data'  , 'CallBack', cominterpdata);
    uimenu( submenu, 'Label', 'FastICA'  , 'CallBack', comfastica, 'separator', 'on' );
    uimenu( submenu, 'Label', 'Component classification (TESA)'  , 'CallBack', comcompselect);
    uimenu( submenu, 'Label', 'Plot and remove components'  , 'CallBack', comcompplot);
    uimenu( submenu, 'Label', 'Enhanced deflation method (EDM)'  , 'CallBack', comedm);
    uimenu( submenu, 'Label', 'PCA compression'  , 'CallBack', compcacompress);
    uimenu( submenu, 'Label', 'PCA suppression'  , 'CallBack', compcasuppress);
    uimenu( submenu, 'Label', 'Detrend data'  , 'CallBack', comdetrend);
    uimenu( submenu, 'Label', 'SOUND algorithm'  , 'CallBack', comSOUND, 'separator', 'on');
    uimenu( submenu, 'Label', 'SSP-SIR'  , 'CallBack', comSSPSIR);
    uimenu( submenu, 'Label', 'Butterworth filter'  , 'CallBack', comfiltbutter, 'separator', 'on');
    uimenu( submenu, 'Label', 'Median filter'  , 'CallBack', comfiltmedian);
    uimenu( submenu, 'Label', 'Extract TEPs'  , 'CallBack', comtepextract, 'separator', 'on');
    uimenu( submenu, 'Label', 'Find and analyse TEP peaks'  , 'CallBack', compeakanalysis);
    uimenu( submenu, 'Label', 'Output TEP peak analysis'  , 'CallBack', compeakoutput);
    uimenu( submenu, 'Label', 'Plot data'  , 'CallBack', complot);
    uimenu( submenu, 'Label', 'Output TEP peak analysis (group)'  , 'CallBack', compeakoutputgroup, 'separator', 'on');
    uimenu( submenu, 'Label', 'Plot data (group)'  , 'CallBack', complotgroup);
    