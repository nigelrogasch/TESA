% tesa_plot() - plots the TMS-evoked activity averaged over trials 
%
% Usage:
%   >>  EEG = tesa_plot( EEG, varargin );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs:
%   'refract',int   - int defines the refractory period (the time for the
%                   TMS artefact to recover below the rate of change). int
%                   is in ms. 
%                   default = 3
%   'rate',int      - int defines the rate of change of the TMS artifact in
%                   uV/ms. 
%                   default = 1e5
%   'tmsLabel','str'- 'str' is a string for the single TMS label.  
%                   default = 'TMS'
%  
% Input pairs for detecting paired pulses
%   'paired','str'  - required. 'str' - type 'yes' to turn on paired detection
%                   default = 'no'
%   'ISI', [int]    - required. [int] is a vector defining interstimulus intervals
%                   between conditioning and test pulses. Multiple ISIs can 
%                   be defined as [1,2,...]. 
%                   default = []
%   'pairLabel',{'str'} - required if more than 1 ISI. {'str'} is a cell array
%                   containing string labels for different ISI conditions.  
%                   Multiple labels can be defined as {'SICI','LICI',...}.
%                   The number of labels defined must equal the number of
%                   ISI conditions defined.
%                   default = {'TMSpair'}
% 
%  Input pairs for detecting repetitive TMS trains
%  'repetitive','str' - required. 'str' - type 'yes' to turn on repetitive detection
%                   default = 'no'
%   'ITI', int      - required. int defines the inter-train interval in ms.
%                   For example, if a 10 Hz rTMS condition is used with 4s
%                   of stimulation (40 pulses) and 26s of rest, ITI = 2600;
%                   default = []
%   'pulseNum', int - required. int defines the number of pulses in a
%                   train. Using the above example, this would be 40. 
%                   deafult = []
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

function tesa_plot( EEG, varargin )

if nargin < 1
	error('Not enough input arguments.');
end

%define defaults
options = struct('xLim',[-100,300],'yLim',[],'elec',[],'CI','off','input','data','roiName','R1','plotPeak','off');

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

t = figure;

%Plot figure
if strcmpi(options.input,'data')
    if isempty(options.elec)
        plot(EEG.times,mean(EEG.data,3),'b'); hold on;
    elseif ~isempty(options.elec);
        for z = 1:EEG.nbchan;
            chan{1,z} = EEG.chanlocs(1,z).labels;
        end;
        elecNum = find(strcmpi(options.elec,chan));
        plot(EEG.times,mean(EEG.data(elecNum,:,:),3)); hold on;
    end
elseif strcmpi(options.input,'ROI')
    plot(EEG.ROI.(options.roiName).time,EEG.ROI.(options.roiName).tseries); hold on;
elseif strcmpi(options.input,'GMFA')
    plot(EEG.GMFA.(options.roiName).time,EEG.GMFA.(options.roiName).tseries); hold on;
end

%Figure settings
if isempty(options.yLim)
    set(gca,'XLim',options.xLim,'box','off','tickdir','out')
elseif ~isempty(options.yLim)
    set(gca,'XLim',options.xLim,'YLim',options.yLim,'box','off','tickdir','out')
end

if strcmpi(options.input,'GMFA')
    ylabel('GMFA (\muV)');
else
    ylabel('Amplitude (\muV)');
end
xlabel('Time (ms)');

if strcmp(options.input,'ROI')
    title('Region of interest')
elseif strcmp(options.input,'GMFA')
    title('Global mean field amplitude')
end

%Confidence intervals (if requested)
if strcmp(options.CI,'on')
    if strcmpi(options.input,'data')
        if isempty(options.elec)
            error('Confidence intervals only available for single channel plots')
        end
        M = mean(EEG.data(elecNum,:,:),3); 
        CI = 1.96*(std(EEG.data(elecNum,:,:),0,3)./(sqrt(size(EEG.data,3))));
        f = fill([EEG.times,fliplr(EEG.times)],[M-CI,fliplr(M+CI)],'b');
        set(f,'FaceAlpha',0.3);set(f,'EdgeColor', 'none');
    elseif strcmpi(options.input,'ROI')
        M = EEG.ROI.(options.roiName).tseries; 
        CI = EEG.ROI.(options.roiName).CI;
        f = fill([EEG.times,fliplr(EEG.times)],[M-CI,fliplr(M+CI)],'b');
        set(f,'FaceAlpha',0.3);set(f,'EdgeColor', 'none');
    elseif strcmpi(options.input,'GMFA')
        error('Confidence intervals are not available for GMFA');
    end
end

%Plot peaks found
if strcmpi(options.plotPeak,'on')
    if strcmpi(options.input,'data')
        error('Peaks are unavailable. Please run pop_tesa_tepextract and pop_tesa_peakanalysis first.');
    else
        peaksTemp = fieldnames(EEG.(options.input).(options.roiName));
        peakLog = strncmp('P',peaksTemp,1) | strncmp('N',peaksTemp,1);
        peak = peaksTemp(peakLog)';
    end

    Ysize = get(gca,'YLim');
    Yheight = (Ysize(1,2)-Ysize(1,1))./10;
    for z=1:size(peak,2);
        if strcmp(options.input,'ROI')
            if strcmp(EEG.(options.input).(options.roiName).(peak{1,z}).found,'yes')
                if strncmp('P',peak{1,z},1)
                    tempColour = 'r';
                else
                    tempColour = 'g';
                end
                tempLat = EEG.(options.input).(options.roiName).(peak{1,z}).lat;
            elseif strcmp(EEG.(options.input).(options.roiName).(peak{1,z}).found,'no')
                tempColour = 'k';
                tempLat = EEG.(options.input).(options.roiName).(peak{1,z}).peak;
            end
        elseif strcmp(options.input,'GMFA')
            if strcmp(EEG.(options.input).(options.roiName).(peak{1,z}).found,'yes')
                if mod(z,2)
                    tempColour = 'r';
                else
                   tempColour = 'g';
                end
                tempLat = EEG.(options.input).(options.roiName).(peak{1,z}).lat;
            elseif strcmp(EEG.(options.input).(options.roiName).(peak{1,z}).found,'no')
                tempColour = 'k';
                tempLat = EEG.(options.input).(options.roiName).(peak{1,z}).peak;
            end
        end
                

        plot(tempLat,EEG.(options.input).(options.roiName).(peak{1,z}).amp,'+','MarkerEdgeColor',tempColour,'MarkerSize',20,'LineWidth',2);%hold on;
        rectangle('position',[EEG.(options.input).(options.roiName).(peak{1,z}).minWin, EEG.(options.input).(options.roiName).(peak{1,z}).amp-(Yheight./2),EEG.(options.input).(options.roiName).(peak{1,z}).maxWin-EEG.(options.input).(options.roiName).(peak{1,z}).minWin,Yheight], 'LineStyle','--','EdgeColor',tempColour);
        pos=[EEG.(options.input).(options.roiName).(peak{1,z}).peak,EEG.(options.input).(options.roiName).(peak{1,z}).amp+((Yheight/2)+(Yheight/4))];
        text(EEG.(options.input).(options.roiName).(peak{1,z}).peak,0,peak{1,z},'Position',pos, 'HorizontalAlignment', 'center');
        centreLine = EEG.(options.input).(options.roiName).(peak{1,z}).amp-(Yheight/2):0.01:EEG.(options.input).(options.roiName).(peak{1,z}).amp+(Yheight/2);
        centreLat = ones(1,size(centreLine,2))*EEG.(options.input).(options.roiName).(peak{1,z}).peak;
        plot(centreLat,centreLine,'LineStyle','--','Color',tempColour);
    end
end

%Plot timing of TMS pulse
plot([0 0], get(gca,'YLim'),'k--');

fprintf('TEP plot generated. \n');

end
