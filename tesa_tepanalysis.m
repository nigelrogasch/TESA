% tesa_dataanalysis() - averages over trials and electrodes in a region of
%                   interest to extract TMS-evoked potentials. Outputs are
%                   saved in the EEG structure (EEG.ROI). For TEPs
%                   following paired pulses, an additional file can be
%                   specified to subtract from the conditioning pulse and
%                   thereby minimising the impact of on-going activity in  the test
%                   pulse time period.
%
% Usage:
%   >>  EEG = tesa_dataanalysis( EEG, type );
%   >>  EEG = tesa_dataanalysis( EEG, type, varargin );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   type            - string for analysis type. Either 'ROI' or 'GMFA'
%   'roi','{}'      - [required for type 'ROI'] input put includingstring or 
%                   cell array defining electrodes to be
%                   averaged for ROI analysis. Type 'all' to average over
%                   all electrodes.
%                   Examples: 'roi','C3' (single); 'roi',{'C3','C4','CP1'} (multiple); 
%                       'roi','all' (all electrodes).
% 
% Optional input pairs (varargin):
%   'roiName','str' - 'str' is a name to identify ROI analysis. This is
%                   useful if multiple different ROIs are to be analysed.
%                   The output will be stored under this name in EEG.ROI.
%                   Example: 'motor'
%                   Defaults are: R1,R2,R....
% 
% Optional input pairs for correcting paired pulses (varargin):
%   'pairCorrect','on' - turns on pair correction.
%   'ISI',int       - int is a number which defines the inter-stimulus
%                   interval between the conditioning and test pulse. The
%                   single TEP that is subtracted from the paired TEP will
%                   be shifted by this many ms to align with the
%                   conditioning pulse. int is in ms.
%                   Example: 100
%   'fileName','str' - 'str' is the path and name of the .set file to be
%                   subtracted from the paired TEP. This should be a single
%                   TEP evoked by a stimulus intensity equivalent to the
%                   conditioning pulse.
%                   Example: 'C:\tmseeg\myfile.set'
% 
%   Further reading on importance of correcting paired pulses with TMS-EEG:
%   Rogasch N.C. et al (2015) Cortical inhibition of distinct mechanisms in 
%   the dorsolateral prefrontal cortex is related to working memory performance:    
%   A TMS–EEG study. Cortex, 64:68-77.
%   In particular, see supplementary materials
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

function EEG = tesa_tepanalysis( EEG, type, varargin )

if nargin < 2
	error('Not enough input arguments.');
end

%define defaults
options = struct('roi',[],'roiName',[],'pairCorrect','off','ISI',[],'fileName',[]);

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

%If ROI selected check that roi exists
if strcmp(type,'ROI')
    if isempty(options.roi)
        error('No electrodes given for region of interest analysis. Please provide electrodes e.g. ''roi'',{''C3'',''FC1''}.');
    end
end

%Convert roi to cell
if strcmp('char',class(options.roi))
    options.roi = {options.roi};
end

%Check for other ROI analysis
if isempty(options.roiName)  
    if ~isfield(EEG,'ROI')
        options.roiName = 'R1';
    else
        options.roiName = ['R',num2str(size(EEG.ROI,2)+1)];
    end
end

%Calculate average over trials
avgTrials=nanmean(EEG.data,3); %Calculate average over trials

if strcmp(type,'ROI')

    %Extract electrodes to be averaged
    e = struct2cell(EEG.chanlocs);
    elec = squeeze(e(1,1,:));

    if strcmp(options.roi{1,1}, 'all')
        timeSeries = avgTrials;
    else
        for a = 1:size(options.roi,2)
            if isempty(find(strcmp(options.roi{1,a},elec)))
                error('%s is not present in the current file. Please choose again.',options.roi{1,a});
            else
                eNum (1,a) = find(strcmp(options.roi{1,a},elec));
            end
        end

        timeSeries = avgTrials(eNum,:);
    end

    %Averages over timeseries
    EEG.ROI.(options.roiName).tseries = nanmean(timeSeries,1);
    EEG.ROI.(options.roiName).chans = options.roi;
    if strcmp(options.roi,'all')
        EEG.ROI.(options.roiName).chans = elec';
    end

    %display message
    fprintf('Region of interest extracted and saved in EEG.ROI.%s\n',options.roiName);

    %Perform correction for paired pulse analysis
    if strcmp(options.pairCorrect, 'on')

        %Check for correct inputs
        if isempty(options.ISI)
            error('Inter-stimulus interval not defined. Please provide as e.g. ''ISI'',3 in function.');
        end
        if isempty(options.fileName)
            error('File containing data for subtraction not defiend. Please provide as e.g. ''fileName'', ''C:\tmseeg\myfile.set'' in function.');
        end

        %Load file for correction
        EEG1 = pop_loadset('filename',options.fileName);

        %Check that channels for ROI correction match between files
        e1 = struct2cell(EEG1.chanlocs);
        elec1 = squeeze(e1(1,1,:));

        if strcmp(options.roi{1,1},'all')
            if sum(strcmp(elec,elec1)) ~= size(elec,1)
                error('The channels in the existing file and the file to be subtracted are not the same. Please ensure these match. Script terminated.');
            end
        else
            for a = 1:size(options.roi,2)
                if sum(strcmp(options.roi{1,a},elec1)) == 0
                    error('Electrode %s does not exist in file to be subtracted. Script terminated.',options.roi{1,a});
                end
            end
        end

        %Check that sampling rate matches between files
        if EEG.srate ~= EEG1.srate
            error('Sampling rates are not matched between existing file and file for subtraction. Script terminated.');
        end

        %Perform ROI analysis on file

        %Calculate average over trials
        avgTrials1=nanmean(EEG1.data,3); %Calculate average over trials

        %Extract electrodes to be averaged
        if strcmp(options.roi{1,1}, 'all')
            timeSeries1 = avgTrials1;
        else
            e1 = struct2cell(EEG1.chanlocs);
            elec1 = squeeze(e1(1,1,:));

            eNum1 = [];
            for a = 1:size(options.roi,2)
                eNum1 (1,a) = find(strcmp(options.roi{1,a},elec1));
            end

            timeSeries1 = avgTrials1(eNum1,:);
        end

        %Shift new time series to align with pulse to be subtracted
        ISIS = (EEG1.srate/1000)*options.ISI; %convert ISI to samples
        timeSeries1(:,1:ISIS) = [];
        temp = zeros(1,ISIS);
        tseries1 = nanmean(timeSeries1,1);
        sub = [tseries1, temp];

        %Subtract corrected time series from existing time series
        EEG.ROI.(options.roiName).tseries = EEG.ROI.(options.roiName).tseries - sub;
        EEG.ROI.(options.roiName).corrected = 'yes';

        %display message
        fprintf('Region of interest extracted and saved in EEG.ROI.%s\n',options.roiName);
        fprintf('Paired pulse correction applied for ISI of %d ms\n',options.ISI);

    end
end

if strcmp(type,'GMFA')
       %Calculates GMFA (standard deviation across electrodes at each time point)
    EEG.GMFA.R1.tseries = std(avgTrials,1);

    %display message
    fprintf('GMFA extracted and saved in EEG.GMFA\n');

    %Perform correction for paired pulse analysis
    if strcmp(options.pairCorrect, 'on')

        %Check for correct inputs
        if isempty(options.ISI)
            error('Inter-stimulus interval not defined. Please provide as e.g. ''ISI'',3 in function.');
        end
        if isempty(options.fileName)
            error('File containing data for subtraction not defiend. Please provide as e.g. ''fileName'', ''C:\tmseeg\myfile.set'' in function.');
        end

        %Load file for correction
        EEG1 = pop_loadset('filename',options.fileName);

        %Check that channels for GMFA correction match between files
        e = struct2cell(EEG.chanlocs);
        elec = squeeze(e(1,1,:));

        e1 = struct2cell(EEG1.chanlocs);
        elec1 = squeeze(e1(1,1,:));

        if sum(strcmp(elec,elec1)) ~= size(elec,1)
            error('The channels in the existing file and the file to be subtracted are not the same. Please ensure these match. Script terminated.');
        end

        %Check that sampling rate matches between files
        if EEG.srate ~= EEG1.srate
            error('Sampling rates are not matched between existing file and file for subtraction. Script terminated.');
        end

        %Perform GMFA analysis on file

        %Calculate average over trials
        avgTrials1=nanmean(EEG1.data,3); %Calculate average over trials

        %Shift new time series to align with pulse to be subtracted
        ISIS = (EEG1.srate/1000)*options.ISI; %convert ISI to samples
        avgTrials1(:,1:ISIS) = [];
        temp = zeros(size(avgTrials1,1),ISIS);
        sub = [avgTrials1, temp];
        corrected = avgTrials - sub;

        %Subtract corrected time series from existing time series
        EEG.GMFA.R1.tseries = std(corrected);
        EEG.GMFA.R1.corrected = 'yes';

        %display message
        fprintf('GMFA extracted and saved in EEG.GMFA\n');
        fprintf('Paired pulse correction applied for ISI of %d ms\n',options.ISI);

    end
end

end
