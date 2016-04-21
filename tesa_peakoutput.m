% tesa_peakoutput() - returns the results for the peak analysis in a table
%                       in the workspace and in a figure. This can be calculated on either the
%                       peak latencies determined using 'tesa_peakanalysis',
%                       or on fixed latencies provided by the user. Users can 
%                       also opt to have the average amplitude incorporating data 
%                       points either side of the peak instead of the absolute
%                       peak amplitude. Either the average of the TEP curve or 
%                       the area under the curve (GMFA only) can be calculated. 
%                      
% Usage:
%   >>  output = tesa_peakoutput( EEG )
%   >>  output = tesa_peakoutput( EEG, 'key1', value1... )
%
% Inputs:
%   EEG             - EEGLAB EEG structure
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
%   'tepName', 'str'  - string indicating which specific ROI or GMFA TEP to
%                       give output for. If not indicated, output for all
%                       TEP fields will be given.
%                       Examples: 'R1', 'motor' etc.
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
%     
% Outputs:
%   output             - table with results from peak analysis
% 
% Examples:
%   output = tesa_peakoutput( EEG ); %returns amplitude at individual latencies for all peaks defined with tesa_peakanalysis
%   output = tesa_peakoutput( EEG, 'averageWin', 5 ); %returns amplitude averaged +/- 5 ms from individual latencies for all peaks defined with tesa_peakanalysis
%   output = tesa_peakoutput( EEG, 'averageWin', 10, 'calcType', 'area', 'tepName','parietal' ); %returns area under the curve +/- 10 ms from individual latencies for all peaks defined in parietal region of interest using tesa_peak analysis
%   output = tesa_peakoutput( EEG, 'winType', 'fixed', 'fixedPeak', [30,60,100,180],'averageWin', 5 ); %returns amplitude averaged +/- 5 ms from peaks given in 'fixedPeak'. It is not necessary to run tesa_peakanalysis for this option
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

function output = tesa_peakoutput( EEG, varargin )

if nargin < 1
	error('Not enough input arguments.');
end

%define defaults
options = struct('winType','individual','calcType','amplitude','tepName',[],'averageWin',[],'fixedPeak',[],'tablePlot','on');

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

if strcmp(options.calcType,'area') && ~isfield(EEG,'GMFA')
    error ('GMFA analysis has not been performed. Area under the curve can only be calculated for GMFA. Please run tesa_tepextract and perform GMFA analysis.');
end

if strcmp(options.calcType,'area') && ~isempty(options.tepName)
    if ~isfield(EEG.GMFA,options.tepName)
        error ('GMFA analysis has not been performed for ''tepName'' ''%s''. Area under the curve can only be calculated for GMFA. Please run tesa_tepextract and perform GMFA analysis.',options.tepName);
    end
end

if strcmp(options.calcType,'area') && isempty(options.averageWin)
    error('For area under curve, an analysis time window must also be entered using averageWin. For example ''averageWin'', 5.');
end

% ANALYSIS BASED ON INDIVIDUAL PEAKS
%Check that fields are present
if strcmp(options.winType,'individual')
    
    cnames = {'Peak','Found?','Latency','Amplitude'};
    
    if ~(isfield(EEG,'ROI') || isfield(EEG,'GMFA'))
        error('There is no ROI or GMFA analyses performed on this data. Please run pop_tesa_tepextract and then pop_tesa_peakanalysis.');
    end      
    
    output = [];

    if isfield(EEG, 'ROI')
        if isempty(options.tepName)
            roiNum = fieldnames(EEG.ROI);
        else
            roiNum{1,1} = options.tepName;
        end
        for a = 1:size(roiNum,1)
            if strcmp(options.calcType,'amplitude')
                if isfield(EEG.ROI,roiNum{a,1})
                    fieldNum = fieldnames(EEG.ROI.(roiNum{a,1}));
                    peakNum = fieldNum(logical(strncmpi(fieldNum,'P',1) + strncmpi(fieldNum,'N',1)));
                    numAll = arrayfun(@(x) str2num(peakNum{x}(2:end)), 1:size(peakNum,1));
                    [Y,Z] = sort(numAll);
                    peakNum = peakNum(Z);
                    for b = 1:size(peakNum,1)
                        if isempty(output)
                            num = 1;
                        else 
                            num = size(output,2)+1;
                        end
                        output(num).analysis = 'ROI';
                        output(num).type = 'individual';
                        output(num).name = roiNum{a,1};
                        output(num).peak = peakNum{b,1};
                        output(num).found = EEG.ROI.(roiNum{a,1}).(peakNum{b,1}).found;
                        output(num).lat = EEG.ROI.(roiNum{a,1}).(peakNum{b,1}).lat;
                        if isempty(options.averageWin)
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
                        rnames{1,num} = [output(num).analysis,' ',output(num).name];
                        d{num,1} = output(num).peak;
                        d{num,2} = output(num).found;
                        d{num,3} = output(num).lat;
                        d{num,4} = output(num).amp;
                        
%                     elseif strcmp(options.calcType,'area')
%                         if options.averageWin == 0
%                             output(num).area = trapz(1,EEG.ROI.(roiNum{a,1}).(peakNum{b,1}).amp);
%                         else
%                             if isnan(output(num).lat)
%                                 findLat = EEG.ROI.(roiNum{a,1}).(peakNum{b,1}).peak;
%                             else
%                                 findLat = EEG.ROI.(roiNum{a,1}).(peakNum{b,1}).lat;
%                             end
%                             [val,tp1] = min(abs(EEG.ROI.(roiNum{a,1}).time-(findLat-options.averageWin)));
%                             [val,tp2] = min(abs(EEG.ROI.(roiNum{a,1}).time-(findLat+options.averageWin)));
%                             output(num).area = trapz(EEG.ROI.(roiNum{a,1}).time(1,tp1:tp2),EEG.ROI.(roiNum{a,1}).tseries(:,tp1:tp2));
%                         end          
                    end   
                end
            end
        end
    end

    if isfield(EEG, 'GMFA')
        if isempty(options.tepName)
            roiNum = fieldnames(EEG.GMFA);
        else
            roiNum{1,1} = options.tepName;
        end
        for a = 1:size(roiNum,1)
             if isfield(EEG.GMFA,roiNum{a,1})
                fieldNum = fieldnames(EEG.GMFA.(roiNum{a,1}));
                peakNum = fieldNum(logical(strncmpi(fieldNum,'P',1) + strncmpi(fieldNum,'N',1)));
                numAll = arrayfun(@(x) str2num(peakNum{x}(2:end)), 1:size(peakNum,1));
                [Y,Z] = sort(numAll);
                peakNum = peakNum(Z);
                for b = 1:size(peakNum,1)
                    if isempty(output)
                        num = 1;
                    else 
                        num = size(output,2)+1;
                    end
                    output(num).analysis = 'GMFA';
                    output(num).type = 'individual';
                    output(num).name = roiNum{a,1};
                    output(num).peak = peakNum{b,1};
                    output(num).found = EEG.GMFA.(roiNum{a,1}).(peakNum{b,1}).found;
                    output(num).lat = EEG.GMFA.(roiNum{a,1}).(peakNum{b,1}).lat;
                    if strcmp(options.calcType,'amplitude')
                        if isempty(options.averageWin)
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
                    elseif strcmp(options.calcType,'area')
                        if isempty(options.averageWin)
                            output(num).area = trapz(1,EEG.GMFA.(roiNum{a,1}).(peakNum{b,1}).amp);
                        else
                            if isnan(output(num).lat)
                                findLat = EEG.GMFA.(roiNum{a,1}).(peakNum{b,1}).peak;
                            else
                                findLat = EEG.GMFA.(roiNum{a,1}).(peakNum{b,1}).lat;
                            end
                            [val,tp1] = min(abs(EEG.GMFA.(roiNum{a,1}).time-(findLat-options.averageWin)));
                            [val,tp2] = min(abs(EEG.GMFA.(roiNum{a,1}).time-(findLat+options.averageWin)));
                            output(num).area = trapz(EEG.GMFA.(roiNum{a,1}).time(1,tp1:tp2),EEG.GMFA.(roiNum{a,1}).tseries(:,tp1:tp2));
                        end
                        cnames{1,4}='area';
                    end
                    rnames{1,num} = [output(num).analysis,' ',output(num).name];
                    d{num,1} = output(num).peak;
                    d{num,2} = output(num).found;
                    d{num,3} = output(num).lat;
                    if strcmp(options.calcType,'amplitude')
                        d{num,4} = output(num).amp;
                    elseif strcmp(options.calcType,'area') 
                        d{num,4} = output(num).area;
                    end
                end    
            end
        end
    end
end
    
% ANALYSIS BASED ON FIXED PEAKS
if strcmp(options.winType,'fixed')
    
    cnames = {'Peak','Amplitude'};
    
    if ~(isfield(EEG,'ROI') || isfield(EEG,'GMFA'))
        error('There is no ROI or GMFA analyses performed on this data. Please run pop_tesa_tepextract.');
    end
    
    if isempty(options.fixedPeak)
        error('For option ''winType'', ''fixed'' please also provide ''fixedPeak'',[int]. For example [30, 60, 100] to calculate peak amplitude/area at 30 ms, 60 ms and 100 ms.');
    end
    
    if size(options.fixedPeak,1) > 1
        error('Please use , for entering ''fixedPeak'' - not ;. For example ''fixedPeak'', [30, 60, 100]');
    end
    
    output = [];
    
    options.fixedPeak = sort(options.fixedPeak);

    if isfield(EEG, 'ROI')
        if strcmp(options.calcType,'amplitude')
            if isempty(options.tepName)
                roiNum = fieldnames(EEG.ROI);
            else
                roiNum{1,1} = options.tepName;
            end
            for a = 1:size(roiNum,1)
                 if isfield(EEG.ROI,roiNum{a,1})
    %             fieldNum = fieldnames(EEG.ROI.(roiNum{a,1}));
    %             peakNum = fieldNum(logical(strncmpi(fieldNum,'P',1) + strncmpi(fieldNum,'N',1)));
                    for b = 1:size(options.fixedPeak,2)
                        if isempty(output)
                            num = 1;
                        else 
                            num = size(output,2)+1;
                        end
                        output(num).analysis = 'ROI';
                        output(num).type = 'fixed';
                        output(num).name = roiNum{a,1};
                        output(num).peak = options.fixedPeak(1,b);
        %                 output(num).found = EEG.ROI.(roiNum{a,1}).(peakNum{b,1}).found;
        %                 output(num).lat = EEG.ROI.(roiNum{a,1}).(peakNum{b,1}).lat;
                        [val,tp] = min(abs(EEG.ROI.(roiNum{a,1}).time-options.fixedPeak(1,b)));
                        if isempty(options.averageWin)
                            output(num).amp = EEG.ROI.(roiNum{a,1}).tseries(1,tp);
                        else
                            [val,tp1] = min(abs(EEG.ROI.(roiNum{a,1}).time-(options.fixedPeak(1,b)-options.averageWin)));
                            [val,tp2] = min(abs(EEG.ROI.(roiNum{a,1}).time-(options.fixedPeak(1,b)+options.averageWin)));
                            output(num).amp = mean(EEG.ROI.(roiNum{a,1}).tseries(:,tp1:tp2));
                        end
                        
                        rnames{1,num} = [output(num).analysis,' ',output(num).name];
                        d{num,1} = output(num).peak;
                        d{num,2} = output(num).amp;
%                     elseif strcmp(options.calcType,'area')
%                         if options.averageWin == 0
%                             output(num).area = trapz(1,EEG.ROI.(roiNum{a,1}).tseries(1,tp));
%                         else
%                             [val,tp1] = min(abs(EEG.ROI.(roiNum{a,1}).time-(tp-options.averageWin)));
%                             [val,tp2] = min(abs(EEG.ROI.(roiNum{a,1}).time-(tp+options.averageWin)));
%                             output(num).area = trapz(EEG.ROI.(roiNum{a,1}).time(1,tp1:tp2),EEG.ROI.(roiNum{a,1}).tseries(:,tp1:tp2));
%                         end          
                    end 
                end
            end
        end
    end

    if isfield(EEG, 'GMFA')
        if isempty(options.tepName)
            roiNum = fieldnames(EEG.GMFA);
        else
            roiNum{1,1} = options.tepName;
        end
        for a = 1:size(roiNum,1)
             if isfield(EEG.GMFA,roiNum{a,1})
%             fieldNum = fieldnames(EEG.GMFA.(roiNum{a,1}));
%             peakNum = fieldNum(logical(strncmpi(fieldNum,'P',1) + strncmpi(fieldNum,'N',1)));
                for b = 1:size(options.fixedPeak,2)
                    if isempty(output)
                        num = 1;
                    else 
                        num = size(output,2)+1;
                    end
                    output(num).analysis = 'GMFA';
                    output(num).type = 'fixed';
                    output(num).name = roiNum{a,1};
                    output(num).peak = options.fixedPeak(1,b);
    %                 output(num).found = EEG.GMFA.(roiNum{a,1}).(peakNum{b,1}).found;
    %                 output(num).lat = EEG.GMFA.(roiNum{a,1}).(peakNum{b,1}).lat;
                    [val,tp] = min(abs(EEG.GMFA.(roiNum{a,1}).time-options.fixedPeak(1,b)));
                    if strcmp(options.calcType,'amplitude')
                        if isempty(options.averageWin)
                            output(num).amp = EEG.GMFA.(roiNum{a,1}).tseries(1,tp);
                        else
                            [val,tp1] = min(abs(EEG.GMFA.(roiNum{a,1}).time-(options.fixedPeak(1,b)-options.averageWin)));
                            [val,tp2] = min(abs(EEG.GMFA.(roiNum{a,1}).time-(options.fixedPeak(1,b)+options.averageWin)));
                            output(num).amp = mean(EEG.GMFA.(roiNum{a,1}).tseries(:,tp1:tp2));
                        end
                    elseif strcmp(options.calcType,'area')
                        if isempty(options.averageWin)
                            output(num).area = trapz(1,EEG.GMFA.(roiNum{a,1}).tseries(1,tp));
                        else
                            [val,tp1] = min(abs(EEG.GMFA.(roiNum{a,1}).time-(options.fixedPeak(1,b)-options.averageWin)));
                            [val,tp2] = min(abs(EEG.GMFA.(roiNum{a,1}).time-(options.fixedPeak(1,b)+options.averageWin)));
                            output(num).area = trapz(EEG.GMFA.(roiNum{a,1}).time(1,tp1:tp2),EEG.GMFA.(roiNum{a,1}).tseries(:,tp1:tp2));
                        end 
                        cnames{1,2}='area';
                    end
                    rnames{1,num} = [output(num).analysis,' ',output(num).name];
                    d{num,1} = output(num).peak;
                    if strcmp(options.calcType,'amplitude')
                        d{num,2} = output(num).amp;
                    elseif strcmp(options.calcType,'area') 
                        d{num,2} = output(num).area;
                    end
                end
            end
        end
    end
end
    
%Check that output contains something
if isempty(output)
    error('Peak analyses were not performed on this data. Please run tesa_tep_peakanalysis first, or used fixed latencies.');
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
    f.Name = 'Peak analysis output';
    f.NumberTitle = 'off';
end

%Display message
fprintf('Peak analysis results returned as ''output'' in workspace.\n');

end
