% tesa_peakanalysis() - finds peaks within a time window defined by user.
%                       For either ROI or GMFA analyses.
%                       Either positive or negative peaks are detected.
%                       Peaks are defined as data point which is
%                       larger/smaller than +/- x data points (default = 3,
%                       however this can also be defined by user). Results
%                       are saved in EEG structure under (either EEG.ROI or
%                       EEG.GMFA). Note that either pop_tesa_roianalysis or
%                       pop_tesa_gmfaanalysis must be run prior to this
%                       script.
%
% Usage:
%   >>  EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin )
%   >>  EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin, varargin )
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   input           - string defining whether 'ROI' or 'GMFA' input is used
%   direction       - string defining whether peaks are 'positive' or
%                       'negative'. Use 'positive' for GMFA analysis.
%   peak            - vector defining the peak of interest. 
%                       Example: 25 (one peak), [25,60,..] (multiple peaks)
%   peakWin         - matrix defining the time windows to search for above
%                       peaks. Note that the number of peak windows defined should
%                       equal the number of peaks. 
%                       For example if peak = [25,60,180]; peakWin = [25,45;80,120;180,220]
% 
% Optional input pairs
%   'method','str'  - either 'largest' or 'centre' (default = largest). If multiple
%                       peaks are detected in a window, largest will search for
%                       the largest peak within the time window. Centre will 
%                       search for the peak closest to the latency defined
%                       in peak.
%   'samples',int   - int is an integer defining the number of samples
%                       either side of a peak that defines the peak.
%                       Peaks are defined as data point which is
%                       larger/smaller than +/- int data points (default = 3)
%   'roiName','str' - 'str' is a name to identify ROI analysis. This is
%                       useful if multiple different ROIs are to be analysed.
%                       The output will be stored under this name in EEG.ROI.
%                       Example: 'motor'
%                       Defaults are: R1,R2,R....
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

function EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin, varargin )

if nargin < 5
	error('Not enough input arguments.');
end

%define defaults
options = struct('method','largest','samples',3,'roiName','R1');

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('EXAMPLE needs key/value pairs')
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
   inpName = pair{1}; % make case insensitive

   if any(strcmpi(inpName,optionNames))%looks for known options and replaces these in options
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end

%Check that direction is entered correctly
if ~(strcmpi(direction,'positive') || strcmpi(direction,'negative'))
    error('Direction must be either ''positive'' or ''negative''.')
end

%Check that input is entered correctly
if ~(strcmpi(input,'ROI') || strcmpi(input,'GMFA'))
    error('Input must be either ''ROI'' or ''GMFA''.')
end

%Check that input exists
if strcmpi(input,'ROI')
    if ~isfield(EEG,'ROI')
        error('ROI does not exist. Please run tesa_roianalysis first.')
    end
    input = upper(input);
elseif strcmpi(input,'GMFA')
    if ~isfield(EEG,'GMFA')
        error('GMFA does not exist. Please run tesa_gmfaanalysis first.')
    end
    input = upper(input);
end

%Check that peak is within epoch
for a = 1:size(peak,2) 
    if peak(1,a) < EEG.times(1,1) || peak(1,a) > EEG.times(1,end)
        error('%d is outside of the epoch range (%d - %d ms). Script terminated.',peak(1,a),EEG.times(1,1),EEG.times(1,end));
    end
end

%Check that number of peaks match number of peakWins
if size(peak,2) ~= size(peakWin,1)
    error('The number of time windows defined does not match the number of peaks. For example if peak = [25,60,180]; peakWin = [25,45;80,120;180,220]');
end

%Check that peakWin is entered properly
for a = 1:size(peakWin,1)
    if size(peakWin,2) ~= 2 || isempty(peakWin(a,1)) || isempty(peakWin(a,2))
        error('The time window for searching for peaks is not entered correctly. Please follow this formula: [25,45;80,120;180,220]');
    end
end

%Check that peakpeakWin is within epoch
for a = 1:size(peakWin,1) 
    if peakWin(a,1) < EEG.times(1,1) || peakWin(a,1) > EEG.times(1,end)
        error('%d is outside of the epoch range (%d - %d ms). Script terminated.',peakWin(a,1),EEG.times(1,1),EEG.times(1,end));
    elseif peakWin(a,2) < EEG.times(1,1) || peakWin(a,2) > EEG.times(1,end)
        error('%d is outside of the epoch range (%d - %d ms). Script terminated.',peakWin(a,2),EEG.times(1,1),EEG.times(1,end));
    end
end

%Check that peak is with peakWin
for a=1:size(peak,2)
    if peak(1,a) < peakWin(a,1) || peak(1,a) > peakWin(a,2)
        error('The peak at %d is outside of the defined time window to search for the peak (%d - %d ms). Script termninated.',peak(1,a),peakWin(a,1),peakWin(a,2));
    end
end

%Create peaks field
for a = 1:size(peak,2)
    if strcmpi(direction,'positive')
        peakName{a,1} = ['P' num2str(peak(1,a))];
    elseif strcmpi(direction,'negative')
        peakName{a,1} = ['N' num2str(peak(1,a))];
    end
    EEG.(input).(options.roiName).(peakName{a,1}).peak = peak(1,a);
    EEG.(input).(options.roiName).(peakName{a,1}).minWin = peakWin(a,1);
    EEG.(input).(options.roiName).(peakName{a,1}).maxWin = peakWin(a,2);
end
        

%find the sample values for peaks
for a = 1:size(peak,2)
    [val,tp(a,1)] = min(abs(EEG.times-peak(1,a)));
    [val,tpW(a,1)] = min(abs(EEG.times-peakWin(a,1)));
    [val,tpW(a,2)] = min(abs(EEG.times-peakWin(a,2)));
end

%Find peaks (defined as point where defined number of samples [default = 3] either
%side of the time point are are either smaller [positive] or larger [negative] in amplitude)

for a = 1:size(tp,1)
    latHold = [];
    num = 1;
    for b = tpW(a,1):tpW(a,2)
        for c = 1:options.samples
            tPlus(c,1) = EEG.(input).(options.roiName).tseries(1,b) - EEG.(input).(options.roiName).tseries(1,b+c);
            tMinus(c,1) = EEG.(input).(options.roiName).tseries(1,b) - EEG.(input).(options.roiName).tseries(1,b-c);
        end
        
        if strcmpi(direction,'positive')
            tPlusLog = tPlus > 0;
            tMinusLog = tMinus > 0;
        elseif strcmpi(direction,'negative')
            tPlusLog = tPlus < 0;
            tMinusLog = tMinus < 0;
        end
            
        if sum(tPlusLog) + sum(tMinusLog) == options.samples*2;
            latHold(num,1) = b;
            num = num+1;
        end
    end
    
%Determines whether peak was found or not and calculates either 1)largest
%peak in the window; or 2) peak closest to the peak name (defined by method
% - default is 'largest').
%If no peak was found, the latency in the peak name (e.g. 25 ms for P25) is used.

    if size(latHold,1) == 1;
        EEG.(input).(options.roiName).(peakName{a,1}).found = 'yes';
        EEG.(input).(options.roiName).(peakName{a,1}).lat = EEG.times(1,latHold);
        EEG.(input).(options.roiName).(peakName{a,1}).amp = EEG.(input).(options.roiName).tseries(1,latHold);
    elseif isempty(latHold);
        EEG.(input).(options.roiName).(peakName{a,1}).found = 'no';
        EEG.(input).(options.roiName).(peakName{a,1}).lat = NaN;
        EEG.(input).(options.roiName).(peakName{a,1}).amp = EEG.(input).(options.roiName).tseries(1,tp(a,1));
    elseif size(latHold,1) > 1;
        if strcmpi(options.method,'largest')
            temp = EEG.(input).(options.roiName).tseries(1,latHold);
            [val, tempWin] = max(temp);
            EEG.(input).(options.roiName).(peakName{a,1}).found = 'yes';
            EEG.(input).(options.roiName).(peakName{a,1}).lat = EEG.times(1,latHold(tempWin,1));
            EEG.(input).(options.roiName).(peakName{a,1}).amp = EEG.(input).(options.roiName).tseries(1,latHold(tempWin,1));
        elseif strcmpi(options.method,'centre')
            diff = abs(tp(a,1)-latHold);
            sortMat =[diff latHold];
            sorted = sortrows(sortMat);
            EEG.(input).(options.roiName).(peakName{a,1}).found = 'yes';
            EEG.(input).(options.roiName).(peakName{a,1}).lat = EEG.times(1,sorted(1,2));
            EEG.(input).(options.roiName).(peakName{a,1}).amp = EEG.(input).(options.roiName).tseries(1,sorted(1,2));
        end
    end
    
    %Display message
    if strcmp(EEG.(input).(options.roiName).(peakName{a,1}).found, 'yes')
        fprintf('%s %s peak found with latency of %d ms and amplitude of %d uV.\n',input,peakName{a,1},EEG.(input).(options.roiName).(peakName{a,1}).lat,EEG.(input).(options.roiName).(peakName{a,1}).amp);
    elseif strcmp(EEG.(input).(options.roiName).(peakName{a,1}).found, 'no')
        fprintf('%s %s peak not found. Amplitude at %d ms returned.\n',input,peakName{a,1},peak(1,a));
    end
    
end

end
