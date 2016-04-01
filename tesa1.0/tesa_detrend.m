% tesa_detrend() - detrends the data by fiting and subtracting a function from
%                   each channel. Either a linear (fitted and subtratcted from 
%                   each trial), exponential or double exponential function 
%                   (fitted to average and subtracted from each trial) can be fitted. 
%                   Note that the Curve Fitting Toolbox is required
%                   to run either the exponential fit or the double
%                   exponential fit.
% 
% Usage:
%   >>  EEG = tesa_detrend( EEG, detrend, timeWin );
%
% Inputs (required):
%   EEG             - EEGLAB EEG structure
%   detrend         - string with type of detrend to perform; 'linear' |
%                   'exponetial' | 'double' 
%   timeWin         - vector with time range for detrending [t1,t2]
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples:
%   EEG = tesa_detrend( EEG, 'linear', [11,500]);
%   EEG = tesa_detrend( EEG, 'exponential', [11,500]);
%   EEG = tesa_detrend( EEG, 'double', [11,500]);
% 
% See also:
%   tesa_fastica, tesa_edm, tesa_pcasupress 

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

function EEG = tesa_detrend( EEG, detrend, timeWin )

if nargin < 3
	error('Not enough input arguments.');
end

%Check that two time points have been specified
if size(timeWin,2) ~= 2
	error('Please provide two time values for the detrend window: [start, end]');
end

%Check that time values given for detrend window are in range of data
if timeWin(1,1) < EEG.times(1,1) || timeWin(1,2) > EEG.times(1,end)
    error('Time values for the detrend window are out of data range. Note that time values are in ms.');
end

%Check that the type of detrend is give correctly
if ~(strcmp(detrend,'linear') || strcmp(detrend,'exponential') || strcmp(detrend,'double'))
    error('The type of detrend to apply is incorrect. Please use either ''linear'', ''exponential'' or ''double''.');
end
    
%find the time values to detrend
[val1,tp1] = min(abs(EEG.times-timeWin(1,1)));
[val2,tp2] = min(abs(EEG.times-timeWin(1,2)));

%Check that the detrend window does not over lap with removed data
if sum(mean(EEG.data(:,tp1,:),3)) == 0 || sum(mean(EEG.data(:,tp2,:),3)) == 0
    error('The window for detrending contains 0s. Please ensure that this window does not overlap with the window of removed data.');
end
    
%Extract data and times for fitting
data = double(EEG.data(:,tp1:tp2,:));
time = double(EEG.times(1,tp1:tp2));

%##### Fit linear curve ######

if strcmp(detrend,'linear');
    dataC = zeros(size(data,1),size(data,2),size(data,3));
    fprintf('Start linear detrend.\n');
    for a = 1:size(data,1)
        for b = 1:size(data,3)
            linP = polyfit(time,data(a,:,b),1); %fit the data
            dataC(a,:,b) = data(a,:,b) - polyval(linP,time); %correct the data
        end 
    end
    %display message
    fprintf('Linear detrend complete.\n');
end

%##### Fit exponential curve #####

if strcmp(detrend,'exponential');
    v = ver;
    if ~sum(strcmp('Curve Fitting Toolbox',extractfield(v,'Name')'))
        error('You need the Curve Fitting Toolbox to run this function');
    end
    dataC = zeros(size(data,1),size(data,2),size(data,3));
    fprintf('Start exponential detrend.\n');
    dataMean = mean(data,3);
    for a = 1:size(data,1)
        expP = fit(time',dataMean(a,:)','exp1');% search for solution
        for b = 1:size(data,3)
            dataC(a,:,b) = data(a,:,b) - expP(time)';% correct the data
        end
%         fprintf('%s detrend is done\n',EEG.chanlocs(a).labels);
    end
    fprintf('Exponential detrend complete.\n');
end
     
%##### Fit double exponential #####
if strcmp(detrend,'double');
    v = ver;
    if ~sum(strcmp('Curve Fitting Toolbox',extractfield(v,'Name')'))
        error('You need the Curve Fitting Toolbox to run this function');
    end
    dataC = zeros(size(data,1),size(data,2),size(data,3));
    fprintf('Start double exponential detrend.\n');
    dataMean = mean(data,3);
    for a = 1:size(data,1)
        exp2P = fit(time',dataMean(a,:)','exp2');
        for b = 1:size(data,3)
            dataC(a,:,b) = data(a,:,b) - exp2P(time)';% correct the data
        end
%         fprintf('%s detrend is done\n',EEG.chanlocs(a).labels);
    end
    fprintf('Double exponential detrend complete.\n');
end

%Insert corrected data
EEG.data(:,tp1:tp2,:) = dataC;

end
