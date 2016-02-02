% tesa_peakoutput() - returns the results for the peak analysis in a table
%                       in the workspace. Users can also opt to have the
%                       average amplitude incorporating data points either 
%                       side of the peak instead of the absolute
%                       peak amplitude.
%
% Usage:
%   >>  output = tesa_peakoutput( EEG )
%   >>  output = tesa_peakoutput( EEG, varargin )
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs
%   'averageWin',int  - integer describing a time window +/- the peak (in
%                       ms) in which an average amplitude will be taken. If left
%                       empty, the absolute amplitude at the peak will be
%                       returned.
%     
% Outputs:
%   output             - table with results from peak analysis
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

function output = tesa_peakoutput( EEG, varargin )

if nargin < 1
	error('Not enough input arguments.');
end

%define defaults
options = struct('averageWin',0);

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

%Check that fields are present
if ~(isfield(EEG,'ROI') || isfield(EEG,'GMFA'))
    error('There is no ROI or GMFA analyses performed on this data. Please run pop_tesa_tepextract and then pop_tesa_peakanalysis.');
end

output = [];

if isfield(EEG, 'ROI')
    roiNum = fieldnames(EEG.ROI);
    for a = 1:size(roiNum,1)
        fieldNum = fieldnames(EEG.ROI.(roiNum{a,1}));
        peakNum = fieldNum(logical(strncmpi(fieldNum,'P',1) + strncmpi(fieldNum,'N',1)));
        for b = 1:size(peakNum,1)
            if isempty(output)
                num = 1;
            else 
                num = size(output,2)+1;
            end
            output(num).analysis = 'ROI';
            output(num).name = roiNum{a,1};
            output(num).peak = peakNum{b,1};
            output(num).found = EEG.ROI.(roiNum{a,1}).(peakNum{b,1}).found;
            output(num).lat = EEG.ROI.(roiNum{a,1}).(peakNum{b,1}).lat;
            if options.averageWin == 0
                output(num).amp = EEG.ROI.(roiNum{a,1}).(peakNum{b,1}).amp;
            else
                if isnan(output(num).lat)
                    findLat = EEG.ROI.(roiNum{a,1}).(peakNum{b,1}).peak;
                else
                    findLat = EEG.ROI.(roiNum{a,1}).(peakNum{b,1}).lat;
                end
                [val,tp1] = min(abs(EEG.ROI.(roiNum{a,1}).time-(findLat-options.averageWin)));
                [val,tp2] = min(abs(EEG.ROI.(roiNum{a,1}).time-(findLat+options.averageWin)));
                output(num).amp = mean(EEG.ROI.(roiNum{a,1}).tseries(:,tp1:tp2));
            end
        end    
    end
end

if isfield(EEG, 'GMFA')
    roiNum = fieldnames(EEG.GMFA);
    for a = 1:size(roiNum,1)
        fieldNum = fieldnames(EEG.GMFA.(roiNum{a,1}));
        peakNum = fieldNum(logical(strncmpi(fieldNum,'P',1) + strncmpi(fieldNum,'N',1)));
        for b = 1:size(peakNum,1)
            if isempty(output)
                num = 1;
            else 
                num = size(output,2)+1;
            end
            output(num).analysis = 'GMFA';
            output(num).name = roiNum{a,1};
            output(num).peak = peakNum{b,1};
            output(num).found = EEG.GMFA.(roiNum{a,1}).(peakNum{b,1}).found;
            output(num).lat = EEG.GMFA.(roiNum{a,1}).(peakNum{b,1}).lat;
            if options.averageWin == 0
                output(num).amp = EEG.GMFA.(roiNum{a,1}).(peakNum{b,1}).amp;
            else
                if isnan(output(num).lat)
                    findLat = EEG.GMFA.(roiNum{a,1}).(peakNum{b,1}).peak;
                else
                    findLat = EEG.GMFA.(roiNum{a,1}).(peakNum{b,1}).lat;
                end
                [val,tp1] = min(abs(EEG.GMFA.(roiNum{a,1}).time-(findLat-options.averageWin)));
                [val,tp2] = min(abs(EEG.GMFA.(roiNum{a,1}).time-(findLat+options.averageWin)));
                output(num).amp = mean(EEG.GMFA.(roiNum{a,1}).tseries(:,tp1:tp2));
            end

        end    
    end
end

%Check that output contains something
if isempty(output)
    error('Peak analyses were not performed on this data. Please run pop_tesa_peakanalysis first.');
end

%Display message
fprintf('Peak analysis results returned in workspace.\n');

end
