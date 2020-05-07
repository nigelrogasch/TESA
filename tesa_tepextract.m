% tesa_tepextract() - averages over trials to generate a TMS-evoked potential (TEP). 
%                   Either a region-of-interest analysis, which averages
%                   over selected electrodes, or a global field amplitude
%                   analysis (standard deviation across electrodes at each
%                   time point) can be performed. Outputs are
%                   saved in the EEG structure EEG.ROI or EEG.GMFA respectively. 
%                   For TEPs following paired pulses, an additional file can be
%                   specified to subtract from the conditioning pulse and
%                   thereby minimising the impact of on-going activity in  the test
%                   pulse time period.
% 
%                   Finding peaks and returning amplitudes and latencies is
%                   performed with tesa_peakanalysis and tesa_peakoutput.
%
%                   Further reading on importance of correcting paired pulses with TMS-EEG:
%                   'Rogasch N.C. et al (2015) Cortical inhibition of distinct mechanisms in 
%                   the dorsolateral prefrontal cortex is related to working memory performance:    
%                   A TMS–EEG study. Cortex, 64:68-77.'
%                   In particular, see supplementary materials
%
% Usage:
%   >>  EEG = tesa_tepextract( EEG, type );
%   >>  EEG = tesa_tepextract( EEG, type, 'key1',value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   type            - string for analysis type. Either 'ROI' or 'GMFA'
%   'elecs','{}'    - [required for type 'ROI'] string or 
%                   cell array defining electrodes to be
%                   averaged for ROI analysis. Type 'all' to average over
%                   all electrodes.
%                   Examples: 'elecs','C3' (single); 'elecs',{'C3','C4','CP1'} (multiple); 
%                   'elecs','all' (all electrodes).
% 
% Optional input pairs (varargin):
%   'tepName','str' - 'str' is a name to identify analysis. This is
%                   useful if multiple different analyses (i.e. ROIs) are to be analysed.
%                   The output will be stored under this name in EEG.ROI or EEG.GMFA
%                   Example: 'motor' or 'parietal'
%                   Defaults are: R1,R2,R....
% 
% Optional input pairs for correcting paired pulses (varargin):
%   'pairCorrect','on' - turns on pair correction.
%   'ISI',int       - [required if pairCorrect on] int is an integer which defines 
%                   the inter-stimulus interval between the conditioning and test pulse. 
%                   The single TEP that is subtracted from the paired TEP will
%                   be shifted by this many ms to align with the
%                   conditioning pulse. int is in ms.
%                   Example: 100
%   'fileName','str' - [required if pairCorrect on]'str' is the path and name of the .set file 
%                   to be subtracted from the paired TEP. This should be a single
%                   TEP evoked by a stimulus intensity equivalent to the
%                   conditioning pulse.
%                   Example: 'C:\tmseeg\myfile.set'
%
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples
%   EEG = tesa_tepextract( EEG, 'ROI', 'elecs', {'FC1','FC3','C1','C3'} ); % standard ROI analysis
%   EEG = tesa_tepextract( EEG, 'ROI', 'elecs', {'C1','C3'}, 'tepName','motor' ); % ROI analysis with specific name
%   EEG = tesa_tepextract( EEG, 'ROI', 'elecs', 'C3', 'pairCorrect', 'on', 'ISI', 100, 'fileName', 'C:\tmseeg\LICI.set' ); % paired pulse analysis
%   EEG = tesa_tepextract( EEG, 'GMFA'); % standard GMFA analysis
% 
% See also:
%   tesa_peakanalysis, tesa_peakoutput 

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

% Change log:
% 19.12.2018: Fixed error in extracting channel labels

function EEG = tesa_tepextract( EEG, type, varargin )

if nargin < 2
	error('Not enough input arguments.');
end

%define defaults
options = struct('elecs',[],'tepName',[],'pairCorrect','off','ISI',[],'fileName',[]);

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
    if isempty(options.elecs)
        error('No electrodes given for region of interest analysis. Please provide electrodes e.g. ''elecs'',{''C3'',''FC1''}.');
    end
end

%Convert roi to cell
if strcmp('char',class(options.elecs))
    options.elecs = {options.elecs};
end


if strcmp(type,'ROI')
    
    %Check for other ROI analysis
    if isempty(options.tepName)  
        if ~isfield(EEG,'ROI')
            options.tepName = 'R1';
        else
            options.tepName = ['R',num2str(size(fieldnames(EEG.ROI),1)+1)];
        end
    end
    
    %Replace spaces with underscore
    if ~isempty(options.tepName)
        options.tepName = strrep(options.tepName,' ','_');
    end

    %Check if ROI name already exists
    if isfield(EEG,'ROI')
        if isfield(EEG.ROI,options.tepName)
            error('tepName ''%s'' already exists. Please choose another name',options.tepName);
        end
    end

    %Extract electrodes to be averaged
%     e = struct2cell(EEG.chanlocs);
%     elec1 = squeeze(e(1,1,:));
    elec = {EEG.chanlocs.labels}';

    if strcmp(options.elecs{1,1}, 'all')
        timeSeries = EEG.data;
    else
        eNum = [];
        missing = [];
        for a = 1:size(options.elecs,2)
            if isempty(find(strcmp(options.elecs{1,a},elec)))
                warning('%s is not present in the current file. Electrode not included in average.',options.elecs{1,a});
                missing{1,size(missing,2)+1} = options.elecs{1,a};
            else
                eNum (1,size(eNum,2)+1) = find(strcmp(options.elecs{1,a},elec));
            end

        end
        
        if isempty(eNum)
            error('None of the electrodes selected for the ROI are present in the data.');
        end

        timeSeries = EEG.data(eNum,:,:);
    end

    %Averages over timeseries
    timeSeriesIn = nanmean(timeSeries,1);
    
    EEG.ROI.(options.tepName).tseries = nanmean(timeSeriesIn,3);
    EEG.ROI.(options.tepName).chans = options.elecs;
    if ~strcmp(options.elecs{1,1}, 'all')
        EEG.ROI.(options.tepName).missing = missing;
    end
    EEG.ROI.(options.tepName).time = EEG.times;
    if strcmp(options.elecs,'all')
        EEG.ROI.(options.tepName).chans = elec';
    end
    EEG.ROI.(options.tepName).CI = 1.96*(std(timeSeriesIn,0,3)./(sqrt(size(timeSeriesIn,3))));
    
    %display message
    if strcmp(options.pairCorrect, 'off')
        fprintf('Region of interest extracted and saved in EEG.ROI.%s\n',options.tepName);
    end

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

        if strcmp(options.elecs{1,1},'all')
            if size(elec1,1) ~= size(elec,1)
                error('The channels in the existing file and the file to be subtracted are not the same. Please ensure these match. Script terminated.');
            end
        else
            for a = 1:size(options.elecs,2)
                if sum(strcmp(options.elecs{1,a},elec1)) == 0
                    warning('%s is not present in the single file. Electrode not included in average.',options.elecs{1,a});
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
        if strcmp(options.elecs{1,1}, 'all')
            timeSeries1 = avgTrials1;
        else
            e1 = struct2cell(EEG1.chanlocs);
            elec1 = squeeze(e1(1,1,:));

            eNum1 = [];
            for a = 1:size(options.elecs,2)
                eNum1 (1,a) = find(strcmp(options.elecs{1,a},elec1));
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
        EEG.ROI.(options.tepName).tseries = EEG.ROI.(options.tepName).tseries - sub;
        EEG.ROI.(options.tepName).corrected = 'yes';
        EEG.ROI.(options.tepName).chans = options.elecs;
        EEG.ROI.(options.tepName).time = EEG.times;
        if strcmp(options.elecs,'all')
            EEG.ROI.(options.tepName).chans = elec';
        end
        EEG.ROI.(options.tepName).CI = EEG.ROI.(options.tepName).CI - sub;
        
        %display message
        fprintf('Region of interest extracted and saved in EEG.ROI.%s\n',options.tepName);
        fprintf('Paired pulse correction applied for ISI of %d ms\n',options.ISI);

    end
end

if strcmp(type,'GMFA')
    
    %Check for other GMFA analysis
    if isempty(options.tepName)  
        if ~isfield(EEG,'GMFA')
            options.tepName = 'R1';
        else
            options.tepName = ['R',num2str(size(fieldnames(EEG.GMFA),1)+1)];
        end
    end

    %Replace spaces with underscore
    if ~isempty(options.tepName)
        options.tepName = strrep(options.tepName,' ','_');
    end
    
    %Check if GMFA name already exists
    if isfield(EEG,'GMFA')
        if isfield(EEG.GMFA,options.tepName)
            error('tepName ''%s'' already exists. Please choose another name',options.tepName)
        end
    end
   
    %Calculates GMFA (standard deviation across electrodes at each time point)
    EEG.GMFA.(options.tepName).tseries = std(mean(EEG.data,3));
    EEG.GMFA.(options.tepName).time = EEG.times;

    %display message
    if strcmp(options.pairCorrect,'off')
        fprintf('GMFA extracted and saved in EEG.GMFA.%s\n',options.tepName);
    end

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

        if size(elec1,1) ~= size(elec,1)
            error('The channels in the existing file and the file to be subtracted are not the same. Please ensure these match. Script terminated.');
        end

        %Check that sampling rate matches between files
        if EEG.srate ~= EEG1.srate
            error('Sampling rates are not matched between existing file and file for subtraction. Script terminated.');
        end

        %Perform GMFA analysis on file

        %Calculate average over trials
        avgTrials=nanmean(EEG.data,3); %Calculate average over trials
        avgTrials1=nanmean(EEG1.data,3); %Calculate average over trials

        %Shift new time series to align with pulse to be subtracted
        ISIS = (EEG1.srate/1000)*options.ISI; %convert ISI to samples
        avgTrials1(:,1:ISIS) = [];
        temp = zeros(size(avgTrials1,1),ISIS);
        sub = [avgTrials1, temp];
        corrected = avgTrials - sub;

        %Subtract corrected time series from existing time series
        EEG.GMFA.(options.tepName).tseries = std(corrected);
        EEG.GMFA.(options.tepName).corrected = 'yes';
        EEG.GMFA.(options.tepName).time = EEG.times;

        %display message
        fprintf('GMFA extracted and saved in EEG.GMFA.%s\n',options.tepName);
        fprintf('Paired pulse correction applied for ISI of %d ms\n',options.ISI);

    end
end

end
