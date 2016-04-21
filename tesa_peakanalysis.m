% tesa_peakanalysis() - finds peaks within a time window defined by user 
%                       for either ROI or GMFA analyses.
%                       Either positive or negative peaks are detected.
%                       Peaks are defined as data point which is
%                       larger/smaller than +/- x data points (default = 5,
%                       however this can also be defined by user). Results
%                       are saved in EEG structure under (either EEG.ROI or
%                       EEG.GMFA). If no peak is found in the defined window, 
%                       the amplitude at the latency defined in peak is returned,
%                       and a NaN is returned in latency. The analysis is run on 
%                       all existing outputs from tesa_tepextract (e.g. ROIs or GMFA), unless
%                       the user opts to run the analysis on one specific
%                       ROI or GMFA.
% 
%                       Note that tesa_tepextract must be run prior to this script.
%
% Usage:
%   >>  EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin )
%   >>  EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin, 'key1', value1... )
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   input           - string defining whether 'ROI' or 'GMFA' input is used
%   direction       - string defining whether peaks are 'positive' or
%                       'negative'. Use 'positive' for GMFA analysis.
%   peak            - vector defining the peak of interest. 
%                       Example: 25 (one peak), [25,60,..] (multiple peaks)
%   peakWin         - matrix defining the time windows to search for above
%                       peaks. Minimum and maximum values for time window are 
%                       defined as 15,35 and peak definitions separated by ;
%                       Note that the number of peak windows defined should
%                       equal the number of peaks. 
%                       For example if peak = [25,60,180]; peakWin = [15,35;40,80;160,200]
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
%                       larger/smaller than +/- int data points (default = 5)
%   'tepName','str' - 'str' is a name of a specific ROI to perform the
%                       analysis on. If this is left blank, all ROI/GMFAs
%                       defined by tesa_tepextract are analysed
%                       Example: 'motor'
%                       Defaults are: R1,R2,R....
%     
% Outputs:
%   EEG             - EEGLAB EEG structure
%
%   Examples
%   EEG = tesa_peakanalysis( EEG, 'ROI', 'negative', 100, [80,120] ); %find a negative peak in all ROI analyses at 100 ms searching between 80 and 120 ms.
%   EEG = tesa_peakanalysis( EEG, 'GMFA', 'positive', [30,60,180],[20,40;50,70;170,190] ); %find 3 positive peaks in the GMFA analysis at 30 ms (between 20-40ms), 60 ms (between 50-70 ms), and 180 ms (between 170-190 ms)
%   EEG = tesa_peakanalysis( EEG, 'ROI', 'positive', [25,70], [15,35;60,80], 'method', 'centre', 'samples', 5, 'tepName', 'motor'); %find 2 positive peaks at 25 ms (15-35 ms), and 70 ms (60-80 ms) using the peak closest to the central peak (i.e. 25 ms or 70 ms), defining a peak as a data point that is larger than all data points +/- 5 samples and only for the ROI analysis named 'motor'.
% 
% See also:
%   tesa_tepextract, tesa_peakoutput 

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

function EEG = tesa_peakanalysis( EEG, input, direction, peak, peakWin, varargin )

if nargin < 5
	error('Not enough input arguments.');
end

%define defaults
options = struct('method','largest','samples',5,'tepName',[]);

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
        error('ROI does not exist. Please run tesa_peakanalysis first.')
    end
    input = upper(input);
elseif strcmpi(input,'GMFA')
    if ~isfield(EEG,'GMFA')
        error('GMFA does not exist. Please run tesa_peakanalysis first.')
    end
    input = upper(input);
end

%Check that peak is within epoch
for a = 1:size(peak,2) 
    if peak(1,a) < EEG.times(1,1) || peak(1,a) > EEG.times(1,end)
        error('%d is outside of the epoch range (%d to %d ms). Script terminated.',peak(1,a),EEG.times(1,1),EEG.times(1,end));
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

%Check that peakWin is within epoch
for a = 1:size(peakWin,1) 
    if peakWin(a,1) < EEG.times(1,1) || peakWin(a,1) > EEG.times(1,end)
        error('%d is outside of the epoch range (%d to %d ms). Script terminated.',peakWin(a,1),EEG.times(1,1),EEG.times(1,end));
    elseif peakWin(a,2) < EEG.times(1,1) || peakWin(a,2) > EEG.times(1,end)
        error('%d is outside of the epoch range (%d to %d ms). Script terminated.',peakWin(a,2),EEG.times(1,1),EEG.times(1,end));
    end
end

%Check that peak is with peakWin
for a=1:size(peak,2)
    if peak(1,a) < peakWin(a,1) || peak(1,a) > peakWin(a,2)
        error('The peak at %d is outside of the defined time window to search for the peak (%d to %d ms). Script termninated.',peak(1,a),peakWin(a,1),peakWin(a,2));
    end
end

%If a specific TEP is nominated, checks that TEP exists
if ~isempty(options.tepName)
  if ~isfield(EEG.(input),options.tepName)
    error ('The ''tepName'' ''%s'' does not exist. Please re-run tesa_tepextract to extrat this ROI/GMFA.',options.tepName);
  end
end

if isempty(options.tepName)
    teps = fieldnames(EEG.(input));
else
    teps{1,1} = options.tepName;
end

for x = 1:size(teps,1)
    
    %Create peaks field
    for a = 1:size(peak,2)
        if strcmpi(direction,'positive')
            peakName{a,1} = ['P' num2str(peak(1,a))];
        elseif strcmpi(direction,'negative')
            peakName{a,1} = ['N' num2str(peak(1,a))];
        end
        EEG.(input).(teps{x,1}).(peakName{a,1}).peak = peak(1,a);
        EEG.(input).(teps{x,1}).(peakName{a,1}).minWin = peakWin(a,1);
        EEG.(input).(teps{x,1}).(peakName{a,1}).maxWin = peakWin(a,2);
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
                tPlus(c,1) = EEG.(input).(teps{x,1}).tseries(1,b) - EEG.(input).(teps{x,1}).tseries(1,b+c);
                tMinus(c,1) = EEG.(input).(teps{x,1}).tseries(1,b) - EEG.(input).(teps{x,1}).tseries(1,b-c);
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
            EEG.(input).(teps{x,1}).(peakName{a,1}).found = 'yes';
            EEG.(input).(teps{x,1}).(peakName{a,1}).lat = EEG.times(1,latHold);
            EEG.(input).(teps{x,1}).(peakName{a,1}).amp = EEG.(input).(teps{x,1}).tseries(1,latHold);
        elseif isempty(latHold);
            EEG.(input).(teps{x,1}).(peakName{a,1}).found = 'no';
            EEG.(input).(teps{x,1}).(peakName{a,1}).lat = NaN;
            EEG.(input).(teps{x,1}).(peakName{a,1}).amp = EEG.(input).(teps{x,1}).tseries(1,tp(a,1));
        elseif size(latHold,1) > 1;
            if strcmpi(options.method,'largest')
                temp = EEG.(input).(teps{x,1}).tseries(1,latHold);
                if strcmpi(direction,'positive')
                    [val, tempWin] = max(temp);
                elseif strcmpi(direction,'negative')
                    [val, tempWin] = min(temp);
                end
                EEG.(input).(teps{x,1}).(peakName{a,1}).found = 'yes';
                EEG.(input).(teps{x,1}).(peakName{a,1}).lat = EEG.times(1,latHold(tempWin,1));
                EEG.(input).(teps{x,1}).(peakName{a,1}).amp = EEG.(input).(teps{x,1}).tseries(1,latHold(tempWin,1));
            elseif strcmpi(options.method,'centre')
                diff = abs(tp(a,1)-latHold);
                sortMat =[diff latHold];
                sorted = sortrows(sortMat);
                EEG.(input).(teps{x,1}).(peakName{a,1}).found = 'yes';
                EEG.(input).(teps{x,1}).(peakName{a,1}).lat = EEG.times(1,sorted(1,2));
                EEG.(input).(teps{x,1}).(peakName{a,1}).amp = EEG.(input).(teps{x,1}).tseries(1,sorted(1,2));
            end
        end

        %Display message
        if strcmp(EEG.(input).(teps{x,1}).(peakName{a,1}).found, 'yes')
            fprintf('%s %s %s peak found with latency of %d ms and amplitude of %d uV.\n',input,teps{x,1},peakName{a,1},EEG.(input).(teps{x,1}).(peakName{a,1}).lat,EEG.(input).(teps{x,1}).(peakName{a,1}).amp);
        elseif strcmp(EEG.(input).(teps{x,1}).(peakName{a,1}).found, 'no')
            fprintf('%s %s %s peak not found. Amplitude at %d ms returned.\n',input,teps{x,1},peakName{a,1},peak(1,a));
        end

    end
    
end

end
