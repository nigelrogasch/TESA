% tesa_filtbutter() - filters the data using a zero-phase butterworth
%                   filter. Either a bandpass or bandstop filter can be
%                   implemented. The filter order is defined by the user.
%                   This function uses the matlab butter and filtfilt
%                   functions.
%
% Usage:
%   >>  EEG = tesa_filtbutter( EEG, high, low, ord, type );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   high            - integer (required). Highpass frequency (allow frequencies
%                   above this value with bandpass)
%                   Example: 1
%   low             - integer (required). Lowpass frequency (allow
%                   frequencies below this value with bandpass)
%                   Example: 50
%   ord             - integer (required). Filter order.
%                   Example: 4 (designs a fourth order butterworth filter)
%   type            - 'str' (required). 'bandpass' | 'bandstop'. Designs either
%                   a zero-phase bandpass or bandstop butterworth filter
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = tesa_filtbutter( EEG, 1, 100, 4, 'bandpass' ); %zero-phase, 4th-order bandpass butterworth filter between 1-100 Hz
%   EEG = tesa_filtbutter( EEG, 48, 52, 4, 'bandstop' ); %zero-phase, 4th-order bandstop butterworth filter between 48-52 Hz
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

function EEG = tesa_filtbutter( EEG, high, low, ord, type )

if nargin < 5
	error('Not enough input arguments.');
end

%Check inputs
if ~isnumeric(high)
    error('Input ''high'' needs to be a number, not a string. e.g. 1, not ''1''.')
elseif ~isnumeric(low)
    error('Input ''low'' needs to be a number, not a string. e.g. 100, not ''100''.')
elseif high > low
    error('Input ''high'' needs to be less than input ''low''.')
elseif ~isnumeric(ord)
    error('Input ''ord'' needs to be a number, not a string. e.g. 4, not ''4''.')
elseif ~(strcmp(type,'bandpass') || strcmp(type,'bandstop'))
    error('Input ''type'' needs to be either ''bandpass'' or ''bandstop''.')
end

Fs = EEG.srate;
ordIn = ord/2;

if strcmpi(type,'bandstop')
    type = 'stop';
end

[z1 p1] = butter(ordIn, [high low]./(Fs/2), type);

data = double(EEG.data);
temp = NaN(size(data,1),EEG.pnts,size(data,3));
for x = 1:size(data,1) 
    for y = 1:size(data,3)
        dataFilt1 = filtfilt(z1,p1,data(x,:,y));
        temp(x,:,y) = dataFilt1;
    end
end 

EEG.data = temp;

%display message
fprintf('Data filtered using a %s zero-phase butterworth filter (order = %d) between %d and %d Hz.\n',type,ord,high,low);
    

end
