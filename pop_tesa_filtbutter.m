% pop_tesa_filtbutter() filters the data using a zero-phase butterworth
%                   filter. Either a band pass, band stop, high pass, or low pass 
%                   filter can be implemented. The filter order is defined by the user.
%                   This function uses the matlab butter and filtfilt functions.
%
% Usage:
%   >>  EEG = pop_tesa_filtbutter( EEG );
%   >>  EEG = pop_tesa_filtbutter( EEG, high, low, ord, type );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   high            - integer (required). Highpass frequency (allow frequencies
%                   above this value with bandpass or high pass)
%                   Example: 1
%   low             - integer (required). Lowpass frequency (allow
%                   frequencies below this value with bandpass or low pass)
%                   Example: 50
%   ord             - integer (required). Filter order.
%                   Example: 4 (designs a fourth order butterworth filter)
%   type            - 'str' (required). 'bandpass' | 'bandstop' | 'highpass' | 'lowpass'. Designs either
%                   a zero-phase band pass, band stop, high pass, or low pass butterworth filter
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = pop_tesa_filtbutter( EEG, 1, 100, 4, 'bandpass' ); %zero-phase, 4th-order bandpass butterworth filter between 1-100 Hz
%   EEG = pop_tesa_filtbutter( EEG, 48, 52, 4, 'bandstop' ); %zero-phase, 4th-order bandstop butterworth filter between 48-52 Hz
%   EEG = pop_tesa_filtbutter( EEG, 1, [], 4, 'highpass' ); %zero-phase, 4th-order high pass butterworth filter allowing frequencies above 1 Hz
%   EEG = pop_tesa_filtbutter( EEG, [], 45, 4, 'lowpass' ); %zero-phase, 4th-order low pass butterworth filter allowing frequencies below 45 Hz
% 
% See also:
%   eegfiltnew

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

% Change log
% 20.6.2018 - included low and high pass filter options
% 20.6.2018 - changed to filtering across concatenated trials instead of
% individual trials to try and minimise edge artifacts and speed up script

 function [EEG com] = pop_tesa_filtbutter( EEG, high, low, ord, type )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
       
    geometry = {[0.6 0.3 0.3] [0.6 0.3 0.3] [0.6 0.6] [0.6 0.6]};

    uilist = {{'style', 'text', 'string', 'Apply butterworth filter','fontweight','bold'} ...
              {'style', 'text', 'string', 'High-pass'} ...
              {'style', 'text', 'string', 'Low-pass'} ...
              {'style', 'text', 'string', 'Frequencies for filtering (Hz)'} ...
              {'style', 'edit', 'string', '1'} ...
              {'style', 'edit', 'string', '100'} ...
              {'style', 'text', 'string', 'Filter order [required]'} ...
              {'style', 'edit', 'string', '4'}...
              {'style', 'text', 'string', 'Filter type.'} ...
              {'style', 'popupmenu', 'string', 'band-pass|band-stop|high-pass|low-pass', 'tag', 'interp' }};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Filter data -- pop_tesa_filtbutter()', 'helpcom', 'pophelp(''pop_tesa_filtbutter'')');
    if isempty(result), return; end;
    
    %Extract data 
    high = str2num(result{1,1});
    low = str2num(result{1,2});
    ord = str2num(result{1,3});
    if result{1,4} == 1
        type = 'bandpass';
    elseif result{1,4} == 2
       	type = 'bandstop';
    elseif result{1,4} == 3
        type = 'highpass';
    elseif result{1,4} == 4
        type = 'lowpass';
    end
    
end

if (strcmp(type,'bandpasss') || strcmp(type,'bandstop')) && (isempty(high) || isempty(low))
    error('Please enter both a high-pass and low-pass value.');
elseif isempty(ord)
    error('Please enter a filter order value.');
end

%Run script from input
EEG = tesa_filtbutter(EEG,high,low,ord,type);
com = sprintf('%s = pop_tesa_filtbutter( %s, %d, %d, %d, ''%s'' );', inputname(1), inputname(1), high, low, ord, type );

end
