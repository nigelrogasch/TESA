% tesa_peakoutputgroup() - returns the results for the peak analysis across 
%                       a group of particiapants in a table in the workspace 
%                       and in a figure. This can be calculated on either the
%                       peak latencies determined using 'tesa_peakanalysis',
%                       or on fixed latencies provided by the user. Users can 
%                       also opt to have the average amplitude incorporating data 
%                       points either side of the peak instead of the absolute
%                       peak amplitude. Either the average of the TEP curve or 
%                       the area under the curve (GMFA only) can be calculated. 
%                       
%                       To use this function, all of the participant files
%                       must be in the same folder and all files must have
%                       undergone identical tesa_peak_analysis runs.  No 
%                       other files should be in this folder. Open
%                       one file from the folder and run the function on
%                       this file and all files will be analysed.
%                      
% Usage:
%   >>  output = tesa_peakoutputgroup( EEG, tepType, tepName )
%   >>  output = tesa_peakoutput( EEG, tepType, tepName, 'key1', value1... )
% 
% Inputs:
%   EEG               - EEGLAB EEG structure
%   tepType           - 'ROI' | 'GMFA' Indicate whether to perform analysis on
%                       ROI or GMFA
%   tepName           - string indicating which specific ROI or GMFA TEP to
%                       give output for.
%                       Examples: 'R1', 'motor' etc. 
% 
% Optional input pairs
%   'winType', 'str'  - 'individual' | 'fixed'. Calculates values using either
%                       the latencies determined for each individual
%                       participant (using tesa_peakanalysis) or set latencies 
%                       provided by the user in 'fixedPeak' (see below).
%                       Default = 'individual'  
%   'calcType', 'str' - 'amplitude' | 'area'. Indicates whether to return the
%                       average amplitude of the TEP time series or the area under
%                       the curve (area under curve only for GMFA analysis).
%                       For area under curve, an analysis time window must 
%                       also be entered using averageWin.
%                       Default = 'amplitude'
% 
%   'averageWin', int - Integer describing a time window +/- the peak (in ms)
%                       in which an average amplitude/area will be taken. 
%                       If left empty, the absolute amplitude at the peak latency 
%                       will be returned. A value is required for
%                       calculating area under the curve for GMFA.
%                       Example: [] - return value at peak; 5 - return
%                       value averaged 5 ms either side of peak.
%   'fixedPeak', [int]- required for 'winType' fixed. Integer or vector describing 
%                       fixed latencies for calculating average or area under 
%                       the curve.
%                       Examples: [30], [30, 60, 100] etc.
%   'tablePlot','str' - 'on' | 'off'. Plots a table with results from peak
%                       analysis.
% Outputs:
%   output             - table with results from peak analysis
% 
% Examples:
%   output = tesa_peakoutputgroup( EEG, 'ROI', 'frontal' ); %returns amplitude at individual latencies for all peaks defined with tesa_peakanalysis in ROI named frontal.
%   output = tesa_peakoutputgroup( EEG, 'ROI', 'frontal', 'averageWin', 5 ); %returns amplitude averaged +/- 5 ms from individual latencies for all peaks defined with tesa_peakanalysis  in ROI named frontal
%   output = tesa_peakoutputgroup( EEG, 'GMFA', 'R1', 'averageWin', 10, 'calcType', 'area'); %returns area under the curve +/- 10 ms from individual latencies for all peaks defined in GMFA using tesa_peak analysis
%   output = tesa_peakoutputgroup( EEG, 'ROI', 'frontal', 'winType', 'fixed', 'fixedPeak', [30,60,100,180],'averageWin', 5 ); %returns amplitude averaged +/- 5 ms from peaks given in 'fixedPeak'. It is not necessary to run tesa_peakanalysis for this option
% 
% See also:
%   tesa_tepextract, tesa_peakanalysis 

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

function output = tesa_peakoutputgroup( EEG, tepType, tepName, varargin )

if nargin < 3
	error('Not enough input arguments.');
end

%define defaults
options = struct('winType','individual','calcType','amplitude','averageWin',[],'fixedPeak',[],'tablePlot','on');

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

if sum(isfield(EEG,tepType)) == 0
    error('Analysis type ''%s'' is not present in the current data. Please run tesa_tepextract first.',tepType);
end

if sum(isfield(EEG.(tepType),tepName)) == 0
    error('TEP of name ''%s'' is not present in the analysis type ''%s''. Please run tesa_tepextract first.',tepName,tepType);
end

if strcmp(options.calcType,'area') && ~isfield(EEG,'GMFA')
    error ('GMFA analysis has not been performed. Area under the curve can only be calculated for GMFA. Please run tesa_tepextract and perform GMFA analysis.');
end

if strcmp(options.calcType,'area') && ~isempty(tepName)
    if ~isfield(EEG.GMFA,tepName)
        error ('GMFA analysis has not been performed for ''tepName'' ''%s''. Area under the curve can only be calculated for GMFA. Please run tesa_tepextract and perform GMFA analysis.',tepName);
    end
end

if strcmp(options.calcType,'area') && isempty(options.averageWin)
    error('For area under curve, an analysis time window must also be entered using averageWin. For example ''averageWin'', 5.');
end

%Get information from other files
fileInfo = dir([EEG.filepath,filesep,'*.set']);

%Sort the file names by numerical markers
names = {fileInfo.name};
maxlen = max(cellfun(@length, names));
padname = @(s) sprintf(['%0' num2str(maxlen) 's'], s);
namesPadded = cellfun(padname, names, 'UniformOutput', false);
[~, sortOrder] = sort(namesPadded);
fileInfo = fileInfo(sortOrder);

fileNames = [];
for x = 1:size(fileInfo,1)
    fileNames{x,1} = fileInfo(x).name;
end

%Set rownames
rnames = [];
for x = 1:size(fileNames,1)
    rnames{1,x} = ['Participant ', num2str(x)];
end

%Set columnames
cnames = [];
if strcmp(options.winType,'individual')
    teps = fieldnames(EEG.(tepType).(tepName));
    tepLog = strncmp('P',teps,1) | strncmp('N',teps,1);
    clabels = teps(tepLog)';
    numAll = arrayfun(@(x) str2num(clabels{x}(2:end)), 1:size(clabels,2));
    [Y,Z] = sort(numAll);
    clabels = clabels(Z);
    for x = 1:size(clabels,2)
        if strcmp(options.calcType,'amplitude')
            clabelsamp{1,x} = [clabels{1,x},'amp'];
        elseif strcmp(options.calcType,'area')
            clabelsamp{1,x} = [clabels{1,x},'area'];
        end
        clabelslat{1,x} = [clabels{1,x},'lat'];
    end
    cnames = [clabelsamp,clabelslat];
elseif strcmp(options.winType,'fixed')
    teps = fieldnames(EEG.(tepType).(tepName));
    tepLog = strncmp('P',teps,1) | strncmp('N',teps,1);
    clabels = teps(tepLog)';
    options.fixedPeak = sort(options.fixedPeak);
    for x = 1:size(options.fixedPeak,2)
        cnames{1,x} = num2str(options.fixedPeak(1,x));
    end
end

%Open files and run tesa_peakoutput
d = [];
for x = 1:size(fileNames,1)
    EEG = pop_loadset( 'filename', fileNames{x,1}, 'filepath', EEG.filepath);
    
    if sum(isfield(EEG,tepType)) == 0
        error('Analysis type ''%s'' is not present in the following data: %s.',tepType, fileNames{x,1});
    end

    if sum(isfield(EEG.(tepType),tepName)) == 0
        error('TEP of name ''%s'' is not present in the analysis type ''%s'' in the following data: %s.',tepName,tepType, fileNames{x,1});
    end
    
    %Run script
    output = tesa_peakoutput(EEG,'tepName',tepName,'calcType',options.calcType,'winType',options.winType,'averageWin',options.averageWin,'fixedPeak',options.fixedPeak,'tablePlot','off');
    
    output1(x).participant = rnames{1,x};
    
    %For calcType individual
    if strcmp(options.winType,'individual')
        for y = 1:size(output,2)
            analysis{1,y} = output(y).analysis;
            peak{1,y} = output(y).peak;
        end

        for a = 1:size(clabels,2)
            present = [];
            for y = 1:size(output,2)
                if strcmp(output(y).analysis,tepType) && strcmp(output(y).peak,clabels(1,a))
                    if strcmp(options.calcType,'amplitude')
                        d(x,a) = output(y).amp;
                        output1(x).(clabelsamp{1,a}) = output(y).amp;
                    elseif strcmp(options.calcType,'area')
                        d(x,a) = output(y).area;
                        output1(x).(clabelsamp{1,a}) = output(y).area;
                    end
                    d(x,a+size(clabels,2)) = output(y).lat;
                    output1(x).(clabelslat{1,a}) = output(y).lat;
                    present = 1;
                end
            end
            if isempty(present)
                error('The following data file does not have the peak %s for %s %s: %s', clabels{1,a}, tepType, tepName, fileNames{x,1});
            end
        end
        
    %For calcType fixed
    elseif strcmp(options.winType,'fixed')
        for a = 1:size(cnames,2)
            for y = 1:size(output,2)
                if strcmp(output(y).analysis,tepType) && strcmp(num2str(output(y).peak),cnames{1,a})
                    if strcmp(options.calcType,'amplitude')
                        d(x,a) = output(y).amp;
                        output1(x).(clabels{1,a}) = output(y).amp;
                    elseif strcmp(options.calcType,'area')
                        d(x,a) = output(y).area;
                        output1(x).(clabels{1,a}) = output(y).area;
                    end
                end
            end
        end
    end
    
end

if strcmpi(options.tablePlot,'on')
    %Plot table with output
    f = figure;
    t = uitable(f,'Data',d,...
                'ColumnName',cnames,... 
                'RowName',rnames);
    t.Position(3) = t.Extent(3);
    t.Position(4) = t.Extent(4);
    f.Position(3) = t.Extent(3)+40;
    f.Position(4) = t.Extent(4)+40;
    f.Name = 'Peak analysis output for group';
    f.NumberTitle = 'off';
    
    hT = uicontrol('style', 'text',... 
    'string', [tepType,' ',tepName,' peak ',options.calcType,' (',options.winType,')'],... 
    'BackgroundColor',f.Color,...
    'fontWeight','bold');
    hT.Position(2) = t.Position(4)+20;
    hT.Position(3) = f.Position(3);
%     'position', compPos,...
        
end

output = output1;

end