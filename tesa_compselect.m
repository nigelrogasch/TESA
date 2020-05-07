% tesa_compselect() - applies a series of rules (determined by the
%                   user) to classify components as artifacts, presents the 
%                   components one-by-one for manual
%                   clarification and then removes components selected as
%                   artifacts from the data. During the clarification step,
%                   users can re-classify components using a drop down
%                   menu. When satisfied, press the 'Next' button to 
%                   continue to the next component. Components are classified as one of the
%                   following: keep, reject, reject - TMS-evoked muscle, reject - blink, 
%                   reject - eye movement, reject - muscle, reject - 
%                   electrode noise, or reject - sensory . If a component is selected
%                   as an artifact using a rule, it is plotted as red, if
%                   not it is plotted as blue. Following removal, a before and
%                   after plot is generated and users can select whether to
%                   continue or repeat the process. Component numbers and
%                   variance represented by each component selected for
%                   removal is stored in the EEG structure =>
%                   EEG.icaCompClassTesaX.
% 
%                   Users are encouraged to report the settings used with
%                   this function in any publications.
% 
%                   Note that the component selection rules are intended
%                   to work as a guide only, and each component must be
%                   manually checked and re-classified if necessary. The rule 
%                   thresholds can be adapted by the user and may differ 
%                   between datasets. There is no guarantee for the success
%                   of these rules in accurately classifying components.
%                   For further information on selecting components, please
%                   read:
%                   
%                   Rogasch NC et al. (2014) Removing artefacts from TMS-EEG 
%                   recordings using independent component analysis: Importance 
%                   for assessing prefrontal and motor cortex network
%                   properties. NeuroImage, 101:425-439
%                   
%                   and
% 
%                   Rogasch NC et al. (2014) Analysing concurrent transcranial magnetic stimulation and
%                   electroencephalographic data: A review and introduction to the open-source
%                   TESA software. NeuroImage, 147:934–951               
%
% Usage:
%   >>  EEG = tesa_compselect( EEG );
%   >>  EEG = tesa_compselect( EEG, 'key1',value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs (varargin):
%   'compCheck','str' - 'on' | 'off'. Turns on or off the component plots
%                   for checking automated component selection. On is highly recommended. 
%                   Default: 'on'
%   'remove','str'  - 'on' | 'off'. Turns on or off removal of components
%                   after classification
%                   Default: 'on'
%   'saveWeights','str' - 'on' | 'off'. Saves the winv (topoplot weights)
%                   and icaact (time course weights) that the
%                   classification was performed on. Note this will
%                   increase the storage size of the files.
%                   Default: 'off'
%   'figSize','str' - 'small' | 'medium' | 'large'; Determines the size of
%                   the figures that are plotted displaying the information
%                   on the components.
%                   Default: 'small'
%   'plotTimeX',[start,end] - vector with integers for plotting the
%                   component time course (in ms).
%                   Default: [-200,500]
%   'plotFreqX',[low,high] - vector with integers for plotting the
%                   frequency distribution of components (in Hz)
%                   Default: [1,100]
%   'freqScale','str' - 'raw' | 'log' | 'log10' | 'db'
%                   y-axis scaling for the frequency plots
%                   Default: 'log'
% 
% Optional input pairs for detecting artifact components (varargin):
%   
%   TMS-evoked muscle activity
%   This type of artifact is detected by comparing the mean absolute amplitude
%   of the component time course within a target window ('tmsMuscleWin') 
%   and the mean absolute amplitude across the entire component time course. A threshold is
%   set by the user for detection (e.g. 8 means the mean absolute amplitude in the target window
%   is 8 times larger than the mean absolute amplitude across the entire time course). 
%   Components are stored under tmsMuscle.
%   'tmsMuscle','str' - 'on' | 'off'. A string which turns on TMS-evoked muscle 
%                   activity detection. 
%                   Default: 'on'
%   'tmsMuscleThresh',int - An integer determining the threshold for
%                   detecting components representing TMS-evoked muscle activity.
%                   Default: 8
%   'tmsMuscleWin',[start,end] - a vector describing the target window for
%                   TMS-evoked muscle activity (in ms).
%                   Default: [11,30]
%   'tmsMuscleFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of TMS-evoked muscle threshold value for each component
%                   in a figure.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
% 
%   Eye blinks
%   This type of artifact is detected by comparing the mean absolute z score
%   (calculated on the component topography weights) of two electrodes close 
%   to the eyes (e.g. 'Fp1' and 'Fp2'). A threhold is set by the user for detection 
%   (e.g. 2.5 - the mean z score of the two electrodes needs to be larger than
%   2.5 for detection). Components are stored under eyes.
%   'blink','str' - 'on' | 'off'. A string which turns on blink detection. 
%                   Default: 'on'
%   'blinkThresh',int - An integer determining the threshold for
%                   detecting components representing eye blinks.
%                   Default: 2.5
%   'blinkElecs',{'str','str'} - a cell with the two electrodes for blink
%                   detection. The two electrodes closest to the eyes are
%                   best for this purpose.
%                   Default: {'Fp1','Fp2'}
%   'blinkFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of blink detection threshold value for each component in a figure.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
% 
%   Lateral eye movement
%   This type of artifact is detected by comparing the z scores
%   (calculated on the component topography weights) of two electrodes on either side
%   of the forehead (e.g. 'F7','F8'). The z score must be positive for one and 
%   negative for the other. A threhold is set by the user for detection 
%   (e.g. 2 means the z score in one electrode must be greater than 2 and less than
%   -2 in the other electrode). Components are stored under eyes.
%   'move','str' - 'on' | 'off'. A string which turns on movement detection. 
%                   Default: 'on'
%   'moveThresh',int - An integer determining the threshold for
%                   detecting components representing lateral eye movement.
%                   Default: 2
%   'moveElecs',{'str','str'} - a cell with the two electrodes for movement
%                   detection. The two electrodes on either side of the
%                   forehead are best for this purpose.
%                   Default: {'F7','F8'}
%   'moveFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of movement detection threshold value for each
%                   component in a figure.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
%
%   Persistent muscle activity
%   This type of artifact is detected by calculating the slope of the line 
%   of best fit between the log EEG power and the log frequency. See the 
%   following manuscript for details:
%
%   Fitzgibbon SP et al. (2016) Automatic determination of EMG-contaminated 
%   components and validation of independent component analysis using EEG during
%   pharmacologic paralysis. Clinical Neurophysiology 127:1781-1793
% 
%   Components are stored under muscle.
%   'muscle','str' - 'on' | 'off'. A string which turns on muscle detection. 
%                   Default: 'on'
%   'muscleThresh',int - An integer determining the threshold for
%                   detecting components representing muscle activity.
%                   Default: -0.31
%   'muscleFreqIn',[int,int] - a vector with the frequency limits for
%                   including in the slope analysis. Leave empty to include
%                   all of the frequency spectra. [7,70] is recommended.
%                   Default: []
%   'muscleFreqEx',[int,int] - a vector with frequencies to exclude (e.g
%                   due to line noise. Leave empty to include all of the
%                   frequency spectra.
%                   Example: [48,52] % Exclude 50 Hz line noise
%                   Default: [];
%   'muslceFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of muscle detection threshold value for each component
%                   in a figure.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
%
%   Note that the legacy TESA heuristic for detecting persitent muscle
%   activity can still be called from the command line. The two following
%   arguments must be called with the tesa_compselect function:
%   'muscleLegacy','str' - 'on' | 'off'. Call this key/value pair to turn
%                   on TESA legacy heuristic (e.g. comparing the mean
%                   values between a target window and the remaining
%                   window using a ratio)
%                   Default: 'off'
%   'muscleFreqWin',[int,int] - a vector with the frequencies for the target w
%                   window(in Hz).
%                   Example: [30,100]
%                   Default: [] 
% 
%   Electrode noise
%   This type of artifact is detected by comparing z scores in individual 
%   electrodes (calculated on the component topography weights). A threhold is 
%   set by the user for detection (e.g. 4 means one or more electrodes has
%   an absolute z score of at least 4). Components are stored in elecNoise.
%   'elecNoise','str' - 'on' | 'off'. A string which turns on electrode noise detection. 
%                   Default: 'on'
%   'elecNoiseThresh',int - An integer determining the threshold for
%                   detecting components representing electrode noise.
%                   Default: 4
%   'elecNoiseFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of electrode noise detection threshold value for each component
%                   in a figure.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
% 
%   Note that a sensory category is also available in the dropdown menu.
%   There are currently no specific rules suitable for accurately detecting
%   these components and are therefore classified entirely at the user's discretion.
% 
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples
%   EEG = tesa_compselect( EEG ); % default use
%   EEG = tesa_compselect( EEG, 'figSize', 'medium','plotTimeX',[-600,600], 'plotFreqX', [2,45] ); % changes the number of components to select, the size of the figure and the x axis' of the time course and frequency plots.
%   EEG = tesa_compselect( EEG, 'elecNoise','off','blinkThresh',3,'blinkElecs',{'AF3','AF4'},'blinkFeedback','on'); % turn off electrode noise detection, change threshold for blinks to 3, change electrodes used to AF3 and AF4 and turn on the feedback of blink threhsolds for  individual components in the command window.
% 
% See also:
%   tesa_fastica, tesa_sortcomps 

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

% Change log
% 13.6.2018 - removed sort components step (added directly to tesa_fastica)
% 13.6.2018 - removed plotting function and replaced with tesa_compplot
% 13.6.2018 - changed the way component classifications are stored. Now stored in EEG.icaCompClassTesaX
% 13.6.2018 - feedback on classification thresholds is now provided as a
%               figure as opposed to on the command line
% 13.6.2018 - changed the classification categories to 'keep' and 'reject' categories
% 14.6.2018 - included 'remove' input to determine whether to classify
%               components ('off'), or classify and remove components ('on')
% 22.3.2019 - Changed the method for calculating the power spectrum to
%               pwelch
% 10.4.2019 - Added an option to choose y-axis frequency scale
% 11.4.2019 - Changed the heuristic rule for detecting persistent muscle
%               activity to the slope of the line of best fit on the
%               log-log of the frequency spectra data (see Fitzgibbon et al
%               2016 Clin Neurophysiology). Note that the legacy TESA
%               muscle heuristic can still be called from the command line

function EEG = tesa_compselect( EEG , varargin )
   
    if nargin < 1
        error('Not enough input arguments.');
    end

    %define defaults
    options = struct('compCheck','on','comps',[], 'remove','on','saveWeights','off','figSize','medium','plotTimeX',[],'plotFreqX',[],'freqScale',[],...
        'tmsMuscle','on','tmsMuscleThresh',[],'tmsMuscleWin',[],'tmsMuscleFeedback','off','blink','on','blinkThresh',[],'blinkElecs',[],...
        'blinkFeedback','off','move','on','moveThresh',[],'moveElecs',[],'moveFeedback','off','muscle','on','muscleThresh',[],'muscleLegacy','off',...
        'muscleFreqWin',[],'muscleFreqIn',[],'muscleFreqEx',[],'muscleFeedback','off','elecNoise','on','elecNoiseThresh',[],'elecNoiseFeedback','off');

    % read the acceptable names
    optionNames = fieldnames(options);

    % count arguments
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('tesa_compselect needs key/value pairs')
    end

    for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
       inpName = pair{1}; % make case insensitive

       if any(strcmpi(inpName,optionNames))%looks for known options and replaces these in options
          options.(inpName) = pair{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end
    
    %Set defaults
    if isempty(options.plotTimeX)
        options.plotTimeX = [-200,500];
    end
    if isempty(options.plotFreqX)
        options.plotFreqX = [1,100];
    end
    if isempty(options.freqScale)
        options.freqScale = 'log';
    end
    if isempty(options.tmsMuscleThresh)
        options.tmsMuscleThresh = 8;
    end
    if isempty(options.tmsMuscleWin)
        options.tmsMuscleWin = [11,30];
    end
    if isempty(options.blinkThresh)
        options.blinkThresh = 2.5;
    end
    if isempty(options.blinkElecs)
        options.blinkElecs = {'Fp1','Fp2'};
    end
    if isempty(options.moveThresh)
        options.moveThresh = 2;
    end
    if isempty(options.moveElecs)
        options.moveElecs = {'F7','F8'};
    end
    if isempty(options.muscleThresh)
        options.muscleThresh = -0.31;
    end
    if isempty(options.elecNoiseThresh)
        options.elecNoiseThresh = 4;
    end
    
    %Check for ICA
    if isempty(EEG.icaweights)
        error('There are no stored components in the data. Please run ICA before running this script. See tesa_fastica.')
    end
    
    %Check for channel locations
    if isempty(EEG.chanlocs(1).theta)
        error('Channel locations are required for this function. Please load channel locations: Edit -> Channel locations.')
    end
    
%     %Check comps input
%     if ~isempty(options.comps)
%         if size(EEG.icaweights,1) < options.comps
%             error('The number of components to choose from (%d) is more than the number of independent components in the data (%d).',options.comps,size(EEG.icaweights,1));
%         end
%     end
    
    %Check figure inputs
    if ~(strcmp(options.figSize,'small') | strcmp(options.figSize,'medium') | strcmp(options.figSize,'large'))
        error ('Input for ''figSize'' needs to be either ''small'', ''medium'' or ''large''.');
    end
    
    %Check plotTimeX inputs
    if size(options.plotTimeX,2) ~= 2
        error('Input for ''plotTimeX'' must be in the following format: [start, end]. e.g. [-200,500].');
    elseif options.plotTimeX(1,1) < EEG.times(1,1) || options.plotTimeX(1,2) > EEG.times(1,end)
        error('Input for ''plotTimeX'' must be within the limits of the data [%d to %d].',EEG.times(1,1),EEG.times(1,end));
    end
    
    %Check plotFreqX inputs
    if size(options.plotFreqX,2) ~= 2
        error('Input for ''plotFreqX'' must be in the following format: [low, high]. e.g. [1,100].');
    elseif options.plotFreqX(1,1)<0
        error('Input for ''plotFreqX'' must be larger than 0.');
    end
    
%     %Set muscle windows based on plotFreqX
%     if isempty(options.muscleFreqWin)
%         options.muscleFreqWin = [30,options.plotFreqX(1,2)];
%     end
    
    %Check frequency scaling inputs
    if ~(strcmp(options.freqScale,'raw') | strcmp(options.freqScale,'log') | strcmp(options.freqScale,'log10') | strcmp(options.freqScale,'db'))
        error ('Input for ''freqScale'' needs to be either ''raw'', ''log'', ''log10'' or ''db''.');
    end
    
    %Checks eye movement input, disables if both electrodes are not present.
    e = struct2cell(EEG.chanlocs);
    elec = squeeze(e(1,1,:));
    eNum = [];
    for a = 1:size(options.moveElecs,2)
        if isempty(find(strcmpi(options.moveElecs{1,a},elec)))
        else
            eNum (1,size(eNum,2)+1) = find(strcmpi(options.moveElecs{1,a},elec));
        end
    end
    
    eyeMove = 'on';
    if size(eNum,2) ~= 2
        warning('One or more of the electrodes required for detecting lateral eye movements was not present. This function has been disabled.');
        eyeMove = 'off';
    end
    
   for a = 1:size(options.blinkElecs,2)
        if isempty(find(strcmpi(options.blinkElecs{1,a},elec)))
            warning('%s is not present in the current file. Electrode not included in blink detection.',options.blinkElecs{1,a});
        else
            eNum (1,size(eNum,2)+1) = find(strcmpi(options.blinkElecs{1,a},elec));
        end

    end

    %Creates output storage name
    varsName = fieldnames(EEG);
    varsNameLog = strcmp('icaCompClass',varsName);
    if sum(varsNameLog) == 0
        EEG.icaCompClass = [];
        nameIn = 'TESA1';
    else
        if isempty(EEG.icaCompClass)
            EEG.icaCompClass = [];
            nameIn = 'TESA1';
        else
            varsNameClass = fieldnames(EEG.icaCompClass);
            varsNameClassLog = strncmp('TESA',varsNameClass,4);        
            if sum(varsNameClassLog) == 0
                nameIn = 'TESA1';
            else 
                nameIn = ['TESA', num2str(sum(varsNameClassLog)+1)];
            end
        end
    end

    %Calculate component time series
    if isempty(EEG.icaact)
        fprintf('Calculating component time series\n')
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
        EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
    end

    %Calculate component variance
    fprintf('Calculating component variance\n')
    for x = 1:size(EEG.icaact,1)
        vars(x) = var(mean(EEG.icaact(x,:,:),3));
    end
    varsPerc = vars/sum(vars)*100;
    
    fprintf('Calculating component power spectrum\n')
    
    % Calculate pwelch
    %Resize EEG.icaact if required
    if size(EEG.icaact,3) > 0
        eegData = reshape(EEG.icaact,size(EEG.icaact,1),[]);
    else
        eegData = EEG.icaact;
    end
    [pxx,fp] = pwelch(eegData',size(eegData,2),[],size(eegData,2),EEG.srate);
    FFTout = pxx';
    fp = fp';
    
    % Calculate FFT bins
    [c index]=min(abs(fp-(0.5)));
    freq=options.plotFreqX(1,1):0.5:options.plotFreqX(1,2);
    fftBins = zeros(size(FFTout,1),size(freq,2)); %preallocate
    for a=1:size(freq,2)
        [c, index1]=min(abs(fp-((freq(1,a)-0.25))));
        [c, index2]=min(abs(fp-((freq(1,a)+0.25))));
        fftBins(:,a)=mean(FFTout(:,index1:index2),2); %creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
    end

%     % LEGACY APPROACH USING FFT
%     %Calculates component power spectrum
%     fprintf('Calculating component power spectrum\n')
%     compNumMat = 1:size(EEG.icaact,1);
%     for compNum = compNumMat
%         %##########
%         %Output the data
% 
%         %Calculates time course
%         temp = EEG.icaact(compNum,:,:);
%         temp = squeeze(temp(1,:,:));
% 
%         %Calculates FFT (uV/Hz)
%         T = 1/EEG.srate;             % Sample time
%         L = size(EEG.times,2);       % Length of signal
%         NFFT = 2^nextpow2(L);        % Next power of 2 from length of y
%         f = EEG.srate/2*linspace(0,1,NFFT/2+1); % Frequencies
% 
%         y = reshape(temp,1,[]);
%         Y = fft(zscore(y),NFFT)/L;
% 
%         FFTout(compNum,:) = (abs(Y).^2); 
%         [c index]=min(abs(f-(0.5)));
% 
%         freq=options.plotFreqX(1,1):0.5:options.plotFreqX(1,2);
% 
%         fftBins = zeros(size(FFTout,1),size(freq,2)); %preallocate
% 
%         for x = 1:size(FFTout,1)
%             for a=1:size(freq,2)
%                 [c, index1]=min(abs(f-((freq(1,a)-0.25))));
%                 [c, index2]=min(abs(f-((freq(1,a)+0.25))));
%                 fftBins(x,a)=mean(FFTout(x,index1:index2),2); %creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
%             end
%         end
% 
%     end
    
    % Generate array for component classification: 1 = keep; >1 = reject
    EEG.icaCompClass.(nameIn).name = nameIn;
    EEG.icaCompClass.(nameIn).compClass = ones(1,size(EEG.icawinv,2));
    
    %Creates output storage for components
    EEG.icaCompClass.(nameIn).reject = []; %2
    EEG.icaCompClass.(nameIn).rejectTmsMuscle = []; %3
    EEG.icaCompClass.(nameIn).rejectBlink = []; %4
    EEG.icaCompClass.(nameIn).rejectEyeMove = []; %5
    EEG.icaCompClass.(nameIn).rejectMuscle = []; %6
    EEG.icaCompClass.(nameIn).rejectElectrodeNoise = []; %7
    EEG.icaCompClass.(nameIn).rejectSensory = []; %8

    EEG.icaCompClass.(nameIn).rejectVars = [];
    EEG.icaCompClass.(nameIn).rejectTmsMuscleVars = [];
    EEG.icaCompClass.(nameIn).rejectBlinkVars = [];
    EEG.icaCompClass.(nameIn).rejectEyeMoveVars = [];
    EEG.icaCompClass.(nameIn).rejectMuscleVars = [];
    EEG.icaCompClass.(nameIn).rejectElectrodeNoiseVars = [];
    EEG.icaCompClass.(nameIn).rejectSensoryVars= [];

    %Number of components to consider
    if isempty(options.comps)
        comps = size(EEG.icaweights,1);
    else
        warning('Input ''comps'' has been removed. Select ''Finish'' and hit ''Next'' to terminate IC selection');
        comps = size(EEG.icaweights,1);
    end
    
    fprintf('Classifying components\n');
    
    for compNum =1:comps

        %##########
        %Checks which components fit profile of artifact

        %Calculates the zscores across electrodes for each component
        tempCompZ = zscore(EEG.icawinv(:,compNum));

        %TMS-evoked muscle
        if strcmpi(options.tmsMuscle,'on')

            %Check TMS muscle inputs
            if options.tmsMuscleThresh < 0
                error('Input for ''tmsMuscleThresh'' must be greater than 0.');
            elseif size(options.tmsMuscleWin,2)~=2
                error('Input for ''tmsMusclesWin'' must be in the following format: [start,end]. e.g. [11,50].');
            elseif options.tmsMuscleWin(1,1) < EEG.times(1,1) || options.tmsMuscleWin(1,2) > EEG.times(1,end)
                error('Input for ''tmsMuscleWin'' must be within the limits of the data [%d to %d].',EEG.times(1,1),EEG.times(1,end));
            end
        end

        [val1,mt1] = min(abs(EEG.times-options.tmsMuscleWin(1,1)));
        [val2,mt2] = min(abs(EEG.times-options.tmsMuscleWin(1,2)));
        muscleScore = abs(mean(EEG.icaact(compNum,:,:),3));
        winScore = mean(muscleScore(:,mt1:mt2),2);
        tmsMuscleRatio(compNum) = winScore./mean(muscleScore);


        %Blink 
        e = struct2cell(EEG.chanlocs);
        elec = squeeze(e(1,1,:));

        eNum = [];
        for a = 1:size(options.blinkElecs,2)
            if isempty(find(strcmpi(options.blinkElecs{1,a},elec)))
%                 warning('%s is not present in the current file. Electrode not included in blink detection.',options.blinkElecs{1,a});
            else
                eNum (1,size(eNum,2)+1) = find(strcmpi(options.blinkElecs{1,a},elec));
            end

        end

        if strcmpi(options.blink,'on')
            if isempty(eNum)
                error('None of the electrodes selected to detect eye blinks are present in the data. Please enter alternative electrodes using ''blinkElecs'',{''elec1'',''elec2''}.');
            end
        end

        blinkRatio(compNum) = mean(tempCompZ(eNum,:));

        %Lateral eye movement
        if strcmpi(options.move,'on')

            %Check TMS muscle inputs
            if options.moveThresh < 0
                error('Input for ''moveThresh'' must be greater than 0.');   
            end

        end

        e = struct2cell(EEG.chanlocs);
        elec = squeeze(e(1,1,:));

        eNum = [];
        if strcmpi(options.move,'on')
            for a = 1:size(options.moveElecs,2)
                if isempty(find(strcmpi(options.moveElecs{1,a},elec)))
%                     error('%s is not present in the current file. Eye movement detection not possible. Please enter alternative electrodes using ''moveElecs'',{''elec1'',''elec2''}.',options.moveElecs{1,a});
                else
                    eNum (1,size(eNum,2)+1) = find(strcmpi(options.moveElecs{1,a},elec));
                end

            end

            if isempty(eNum)
                error('None of the electrodes selected to detect eye movements are present in the data.Please enter alternative electrodes using ''moveElecs'',{''elec1'',''elec2''}.');
            end
        end
        moveVal(compNum,:) = tempCompZ(eNum,:);

        %Muscle
        if strcmpi(options.muscle,'on')
            
            if strcmpi(options.muscleLegacy,'off') && ~isempty(options.muscleFreqWin)
                warning('The heuristic TESA uses to detect persistent muscle activity has changed. Please see the user manual.');
            end
            
            if strcmpi(options.muscleLegacy,'on')
                if isempty(muscleFreqWin)
                    error('To use the legacy muscle heuristic, please also define ''muscleFreqWin''.');
                end
                %Checks input for muscle
                if size(options.muscleFreqWin,2) ~=2
                    error('Inputs for ''muscleFreqWin'' must be in the following format: [low, high]. e.g. [31,100].');
                elseif (options.muscleFreqWin(1,1) < options.plotFreqX(1,1) | options.muscleFreqWin(1,2) > options.plotFreqX(1,2))
                    %                 error('Inputs for ''muscleFreqWin'' (%d to %d) are outside of the frequency range set by input ''plotFreqX'' (%d to %d). Please adjust.',options.muscleFreqWin(1,1),options.muscleFreqWin(1,2),options.plotFreqX(1,1),options.plotFreqX(1,2));
                end
                
                [val1,winF1] = min(abs(freq-options.muscleFreqWin(1,1)));
                [val2,winF2] = min(abs(freq-options.muscleFreqWin(1,2)));
                winFreq = mean(fftBins(compNum,winF1:winF2),2);
                allFreq = mean(fftBins(compNum,:),2);
                muscleRatio(compNum) = winFreq./allFreq;
            else
                
                %Checks input for muscle
                if size(options.muscleFreqIn,2) ~=2
                    error('Inputs for ''muscleFreqIn'' must be in the following format: [low, high]. e.g. [7,75].');
%                 elseif (options.muscleFreqIn(1,1) < freq(1,1) | options.muscleFreqIn(1,2) > freq(1,end))
%                     error('Inputs for ''muscleFreqIn'' (%d to %d) are outside of the frequency range (%d to %d). Please adjust.',options.muscleFreqIn(1,1),options.muscleFreqIn(1,2),freq(1,1),freq(1,end));
                end
%                 if size(options.muscleFreqEx,2) ~=2
%                     error('Inputs for ''muscleFreqEx'' must be in the following format: [low, high]. e.g. [48,52].');
%                 elseif (options.muscleFreqEx(1,1) < freq(1,1) | options.muscleFreqEx(1,2) > freq(1,end))
%                     error('Inputs for ''muscleFreqEx'' (%d to %d) are outside of the frequency range (%d to %d). Please adjust.',options.muscleFreqEx(1,1),options.muscleFreqEx(1,2),freq(1,1),freq(1,end));
%                 end
                
                % Define frequencies to include in the analysis
                if ~isempty(options.muscleFreqIn)
                    [~,fin1] = min(abs(options.muscleFreqIn(1) - freq));
                    [~,fin2] = min(abs(options.muscleFreqIn(2) - freq));
                    
                    freqHz = freq(1,fin1:fin2);
                    freqPow = fftBins(compNum,fin1:fin2);
                else
                    freqHz = freq;
                    freqPow = fftBins(compNum,:);
                end
                
                % Define frequencies to exclude from fit
                if ~isempty(options.muscleFreqEx)
                    [~,fex1] = min(abs(options.muscleFreqEx(1) - freqHz));
                    [~,fex2] = min(abs(options.muscleFreqEx(2) - freqHz));
                    freqHz(fex1:fex2) = [];
                    freqPow(fex1:fex2) = [];
                end
                
                % Fit linear regression to log-log data
                p = polyfit(log(freqHz),log(freqPow),1);
                
                % Store the slope
                muscleRatio(compNum) = p(1);
            end
        end

        %Electrode noise
        elecNoise(compNum) = max(abs(round(tempCompZ,2)));

        %Select if component is artifact
        if strcmpi(options.tmsMuscle,'on') && tmsMuscleRatio(compNum) >= options.tmsMuscleThresh
            EEG.icaCompClass.(nameIn).compClass(compNum)= 3;
        elseif strcmpi(options.blink,'on') && abs(blinkRatio(compNum)) >= options.blinkThresh
            EEG.icaCompClass.(nameIn).compClass(compNum) = 4;
        elseif strcmpi(eyeMove,'on') && strcmpi(options.move,'on') && (moveVal(compNum,1) >= options.moveThresh && moveVal(compNum,2) <= -options.moveThresh) | (moveVal(compNum,2) >= options.moveThresh && moveVal(compNum,1) <= -options.moveThresh)
            EEG.icaCompClass.(nameIn).compClass(compNum) = 5;    
        elseif strcmpi(options.muscle,'on') && muscleRatio(compNum) >= options.muscleThresh
            EEG.icaCompClass.(nameIn).compClass(compNum) = 6;
        elseif strcmpi(options.elecNoise,'on') && elecNoise(compNum) >= abs(options.elecNoiseThresh);
            EEG.icaCompClass.(nameIn).compClass(compNum) = 7;
        else
            EEG.icaCompClass.(nameIn).compClass(compNum) = 1;
        end

        %Save components for removal
        if EEG.icaCompClass.(nameIn).compClass(compNum) == 3
            EEG.icaCompClass.(nameIn).rejectTmsMuscle(1,size(EEG.icaCompClass.(nameIn).rejectTmsMuscle,2)+1) = compNum;
            EEG.icaCompClass.(nameIn).rejectTmsMuscleVars(1,size(EEG.icaCompClass.(nameIn).rejectTmsMuscleVars,2)+1) = varsPerc(1,compNum);
        elseif EEG.icaCompClass.(nameIn).compClass(compNum) == 4
            EEG.icaCompClass.(nameIn).rejectBlink(1,size(EEG.icaCompClass.(nameIn).rejectBlink,2)+1) = compNum;
            EEG.icaCompClass.(nameIn).rejectBlinkVars(1,size(EEG.icaCompClass.(nameIn).rejectBlinkVars,2)+1) = varsPerc(1,compNum);
        elseif EEG.icaCompClass.(nameIn).compClass(compNum) == 5
            EEG.icaCompClass.(nameIn).rejectEyeMove(1,size(EEG.icaCompClass.(nameIn).rejectEyeMove,2)+1) = compNum;
            EEG.icaCompClass.(nameIn).rejectEyeMoveVars(1,size(EEG.icaCompClass.(nameIn).rejectEyeMoveVars,2)+1) = varsPerc(1,compNum);
        elseif EEG.icaCompClass.(nameIn).compClass(compNum) == 6
            EEG.icaCompClass.(nameIn).rejectMuscle(1,size(EEG.icaCompClass.(nameIn).rejectMuscle,2)+1) = compNum;
            EEG.icaCompClass.(nameIn).rejectMuscleVars(1,size(EEG.icaCompClass.(nameIn).rejectMuscleVars,2)+1) = varsPerc(1,compNum);
        elseif EEG.icaCompClass.(nameIn).compClass(compNum) == 7
            EEG.icaCompClass.(nameIn).rejectElectrodeNoise(1,size(EEG.icaCompClass.(nameIn).rejectElectrodeNoise,2)+1) = compNum;
            EEG.icaCompClass.(nameIn).rejectElectrodeNoiseVars(1,size(EEG.icaCompClass.(nameIn).rejectElectrodeNoiseVars,2)+1) = varsPerc(1,compNum);           
        end

    end

    % Save component variance
    EEG.icaCompClass.(nameIn).compVars = varsPerc;
    
    % Save ICA weights and time course
    if strcmp(options.saveWeights,'on')
        EEG.icaCompClass.(nameIn).icawinv = EEG.icawinv;
        EEG.icaCompClass.(nameIn).icaact = EEG.icaact;
    end
    
    feedBack = [];
    % Print feedback figures
    if strcmpi(options.tmsMuscleFeedback,'on')
        
        % Component details
        artRatio = tmsMuscleRatio;
        artThresh = options.tmsMuscleThresh;
        artName = 'TMS-evoked muscle';
        
        timeX = 1:length(artRatio);
        artLog = artRatio >= artThresh;
        
        % Plot figure
        feedBack.t1=figure;
        h1 = stem(timeX,artRatio,'lineWidth',2); hold on;
        
        % Check if any components are above threshold and plot as red
        if sum(artLog)>= 1
            artRatio2 = artRatio(artLog);
            timeX2 = timeX(artLog);
            h2 = stem(timeX2,artRatio2,'lineWidth',2);
            h2(1).Color = 'r';
        end
        
        % Plot threshold line
        plot([1,length(artRatio)],[artThresh artThresh],'r--','lineWidth',2);
        set(gca,'box','off','tickDir','out','lineWidth',2,'fontsize',14,'xlim',[0,length(artRatio)+1])
        title([artName,': ',num2str(sum(artLog)),' of ',num2str(length(artRatio))]);
        xlabel('Components');
        ylabel('Classification value');
    end
    if strcmpi(options.blinkFeedback,'on')
        
        % Component details
        artRatio = abs(blinkRatio);
        artThresh = options.blinkThresh;
        artName = 'Blink';
        
        timeX = 1:length(artRatio);
        artLog = artRatio >= artThresh;
        
        % Plot figure
        feedBack.t2=figure;
        h1 = stem(timeX,artRatio,'lineWidth',2); hold on;
        
        % Check if any components are above threshold and plot as red
        if sum(artLog)>= 1
            artRatio2 = artRatio(artLog);
            timeX2 = timeX(artLog);
            h2 = stem(timeX2,artRatio2,'lineWidth',2);
            h2(1).Color = 'r';
        end
        
        % Plot threshold line
        plot([1,length(artRatio)],[artThresh artThresh],'r--','lineWidth',2);
        set(gca,'box','off','tickDir','out','lineWidth',2,'fontsize',14,'xlim',[0,length(artRatio)+1])
        title([artName,': ',num2str(sum(artLog)),' of ',num2str(length(artRatio))]);
        xlabel('Components');
        ylabel('Classification value');
    end
    if strcmpi(options.moveFeedback,'on')
        
        % Component details
        artRatioA = abs(moveVal(:,1));
        artRatioB = abs(moveVal(:,2));
        artThresh = options.moveThresh;
        artName = 'Eye movement';
        
        timeX = 1:length(artRatioA);
        artLog = (artRatioA >= artThresh) & (artRatioB >= artThresh);
        
        % Plot figure
        feedBack.t3=figure;
        h1a = stem(timeX,artRatioA,'lineWidth',2); hold on;
        h1b = stem(timeX,artRatioB,'lineWidth',2);
        h1b(1).Color = h1a(1).Color;
        
        % Check if any components are above threshold and plot as red
        if sum(artLog)>= 1
            artRatio2a = artRatioA(artLog);
            artRatio2b = artRatioB(artLog);
            timeX2 = timeX(artLog);
            h2a = stem(timeX2,artRatio2a,'lineWidth',2);
            h2b = stem(timeX2,artRatio2b,'lineWidth',2);
            h2a(1).Color = 'r';
            h2b(1).Color = 'r';
        end
        
        % Plot threshold line
        plot([1,length(artRatioA)],[artThresh artThresh],'r--','lineWidth',2);
        set(gca,'box','off','tickDir','out','lineWidth',2,'fontsize',14,'xlim',[0,length(artRatioA)+1])
        title([artName,': ',num2str(sum(artLog)),' of ',num2str(length(artRatioA))]);
        xlabel('Components');
        ylabel('Classification value');    
    end
    if strcmpi(options.muscleFeedback,'on')
        % Component details
        artRatio = muscleRatio;
        artThresh = options.muscleThresh;
        artName = 'Muscle';
        
        timeX = 1:length(artRatio);
        artLog = artRatio >= artThresh;
        
        % Plot figure
        feedBack.t4=figure;
        h1 = stem(timeX,artRatio,'lineWidth',2); hold on;
        
        % Check if any components are above threshold and plot as red
        if sum(artLog)>= 1
            artRatio2 = artRatio(artLog);
            timeX2 = timeX(artLog);
            h2 = stem(timeX2,artRatio2,'lineWidth',2);
            h2(1).Color = 'r';
        end
        
        % Plot threshold line
        plot([1,length(artRatio)],[artThresh artThresh],'r--','lineWidth',2);
        set(gca,'box','off','tickDir','out','lineWidth',2,'fontsize',14,'xlim',[0,length(artRatio)+1])
        title([artName,': ',num2str(sum(artLog)),' of ',num2str(length(artRatio))]);
        xlabel('Components');
        ylabel('Classification value');           
    end
    if strcmpi(options.elecNoiseFeedback,'on')
         % Component details
        artRatio = abs(elecNoise);
        artThresh = options.elecNoiseThresh;
        artName = 'Electrode noise';
        
        timeX = 1:length(artRatio);
        artLog = artRatio >= artThresh;
        
        % Plot figure
        feedBack.t5=figure;
        h1 = stem(timeX,artRatio,'lineWidth',2); hold on;
        
        % Check if any components are above threshold and plot as red
        if sum(artLog)>= 1
            artRatio2 = artRatio(artLog);
            timeX2 = timeX(artLog);
            h2 = stem(timeX2,artRatio2,'lineWidth',2);
            h2(1).Color = 'r';
        end
        
        % Plot threshold line
        plot([1,length(artRatio)],[artThresh artThresh],'r--','lineWidth',2);
        set(gca,'box','off','tickDir','out','lineWidth',2,'fontsize',14,'xlim',[0,length(artRatio)+1])
        title([artName,': ',num2str(sum(artLog)),' of ',num2str(length(artRatio))]);
        xlabel('Components');
        ylabel('Classification value');    
    end
    
    if strcmp(options.remove,'on')
        % Check if component need to be plotted, if not remove bad components
        if strcmp(options.compCheck,'on')
            EEG = tesa_compplot(EEG,'compClassInput',EEG.icaCompClass.(nameIn),'saveWeights',options.saveWeights,'figSize',options.figSize,'plotTimeX',options.plotTimeX,'plotFreqX',options.plotFreqX,'varsPerc',varsPerc,'fftBins',fftBins,'freqBins',freq,'freqScale',options.freqScale);
        else
            %Remove bad components
            allComps = 1:length(EEG.icaCompClass.(nameIn).compClass);
            tempCompBad = EEG.icaCompClass.(nameIn).compClass > 1;
            removeComps = allComps(tempCompBad);
            if ~isempty(removeComps)
                EEG = pop_subcomp( EEG,removeComps, 0);
            end
        end
    end
    
    % Close feedback figure
    if ~isempty(feedBack)
        feedBackAll = fieldnames(feedBack);
        for fx = 1:length(feedBackAll)
            close(feedBack.(feedBackAll{fx}));
        end
    end
end
