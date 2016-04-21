% tesa_medfilt()      - applies a median filter over a given data window.
%                           The median filter takes in to account data either side of 
%                           window according to the filter order and does not zero pad.
% 
% Usage:
%   >>  out = tesa_medfilt( data, win, filtOrd );
%
% Inputs:
%   data                 - data matrix n x m where n is channels and m is
%                       all trial conactenated in to a single time course.
%   win                 - win is a n x m data matrix where n is channels and m
%                       are the time points for filtering
%   filtOrd             - integer describing the filter order. The filter
%                       order determines how many samples are used for the
%                       median filter. The filter order follow the Matlab
%                       median filter convention for odds and evens: 
%                       When filtOrd is odd, y(k) is the median of x(k-(filtOrd-1)/2:k+(filtOrd-1)/2).
%                       When filtOrd is even, y(k) is the median of x(k-filtOrd/2:k+(filtOrd/2)-1).
% 
% Outputs:
%   out                 - data matrix n x m the same size as win with the
%                       filtered data.
%
% See also:
%   tesa_filtmedian 

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

function out = tesa_medfilt(data,win,filtOrd)

out = [];
for a = 1:size(win,2);
    if mod(filtOrd,2) == 0
        filtWin = data(:,win(1,a)-30/2 : win(1,a)+filtOrd/2-1);
    else
        filtWin = data(:,win(1,a)-(filtOrd-1)/2 : win(1,a)+(filtOrd-1)/2);
    end
    out(:,a) = median(filtWin')';
end