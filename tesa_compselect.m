% tesa_compselect() - sorts components following ICA by time course
%                   variance, applies a series of rules (determined by the
%                   user) to classify components as artefacts or neural
%                   activity, presents the components one-by-one for manual
%                   clarification and then removes components selected as
%                   artifacts from the data. During the clarification step,
%                   users can re-classify components using a drop down
%                   menu. When satisfied, press enter (or any other button) to 
%                   continue to the next component. Components are classified as one of the
%                   following: neural, TMS-evoked muscle, eye, muscle,
%                   electrode noise or sensory. If a component is selected
%                   as an artifact using a rule, it is plotted as red, if
%                   neural it is plotted as blue. Following removal, a before and
%                   after plot is generated and users can select whether to
%                   continue or repeat the process. Component numbers and
%                   variance represented by each component selected for
%                   removal is stored in the EEG structure =>
%                   EEG.icaBadComp.
% 
%                   Users are encouraged to report the settings used with
%                   this function in any publications.
% 
%                   Note that the component selection rules are intedended
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
% Usage:
%   >>  EEG = tesa_compselect( EEG );
%   >>  EEG = tesa_compselect( EEG, 'key1',value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs (varargin):
%   'comps',int     - int is an integer describing the number of components
%                   to perform selection on (e.g. first 10 components).
%                   Leave empty for all components. 
%                   Default: []
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
% 
% Optional input pairs for detecting artifact components (varargin):
%   
%   TMS-evoked muscle activity
%   This type of artifact is detected by comparing the standard deviation
%   of the mean time course of the component within a target window
%   ('tmsMuscleWin') and a baseline window ('tmsMuscleBase'). A threhold is
%   set by the user for detection (e.g. 6 means the SD in the target window
%   is 6 times larger than in the base window). Components are stored under
%   tmsMuscle.
%   'tmsMuscle','str' - 'on' | 'off'. A string which turns on TMS-evoked muscle 
%                   activity detection. 
%                   Default: 'on'
%   'tmsMuscleThresh',int - An integer determining the threshold for
%                   detecting components representing TMS-evoked muscle activity.
%                   Default: 6
%   'tmsMuscleWin',[start,end] - a vector describing the target window for
%                   TMS-evoked muscle activity (in ms).
%                   Default: [11,50]
%   'tmsMuscleBase',[start,end] - a vector describing a comparison window
%                   for comparing the target window with (in ms).
%                   Default: [51,100]
%   'tmsMuscleFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of TMS-evoked muscle threshold value for each component
%                   in the command window.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
% 
%   Eye blinks
%   This type of artifact is detected by comparing the mean z score
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
%                   of blink detection threshold value for each component in the command window.
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
%                   component in the command window.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
%
%   Persistent muscle activity
%   This type of artifact is detected by comparing the mean power of the 
%   component time course frequency distribution between a low and high frequency window
%   (calculated using an FFT across all trials). A threhold is set by the user for detection 
%   (e.g. 0.5 means the high frequency power is 50% of the low frequency power). 
%   Components are stored under muscle.
%   'muscle','str' - 'on' | 'off'. A string which turns on muscle detection. 
%                   Default: 'on'
%   'muscleThresh',int - An integer determining the threshold for
%                   detecting components representing muscle activity.
%                   Default: 0.5
%   'muscleFreqLow',[int,int] - a vector with the frequencies for the low
%                   window (in Hz).
%                   Default: [1,30] - Note the default low value will
%                   adjust to match the minimum plot frequency defined by
%                   'plotFreqX'.
%   'muscleFreqHigh',[int,int] - a vector with the frequencies for the high
%                   window (in Hz).
%                   Default: [31,100] - Note the default high value will
%                   adjust to match the maximum plot frequency defined by
%                   'plotFreqX'.
%   'muslceFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of muscle detection threshold value for each component
%                   in the command window.
%                   (Useful for determining a suitable threshold).
%                   Default: 'off'
% 
%   Electrode noise
%   This type of artifact is detected by comparing z scores in individual 
%   electrodes (calculated on the component topography weights). A threhold is 
%   set by the user for detection (e.g. 4 means one or more electrodes has
%   a z score of at least 4). Components are stored in elecNoise.
%   'elecNoise','str' - 'on' | 'off'. A string which turns on electrode noise detection. 
%                   Default: 'on'
%   'elecNoiseThresh',int - An integer determining the threshold for
%                   detecting components representing electrode noise.
%                   Default: 4
%   'elecNoiseFeedback','str' - 'on' | 'off'. String turning on feedback
%                   of electrode noise detection threshold value for each component
%                   in the command window.
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
%   EEG = tesa_compselect( EEG, 'comps', 10, 'figSize', 'medium','plotTimeX',[-600,600], 'plotFreqX', [2,45] ); % changes the number of components to select, the size of the figure and the x axis' of the time course and frequency plots.
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

function EEG = tesa_compselect( EEG , varargin )
   
    if nargin < 1
        error('Not enough input arguments.');
    end

    %define defaults
    options = struct('comps',[], 'figSize','small','plotTimeX',[-200,500],'plotFreqX',[1,100],'tmsMuscle','on','tmsMuscleThresh',6,...
        'tmsMuscleWin',[11,50],'tmsMuscleBase',[51,100],'tmsMuscleFeedback','off','blink','on','blinkThresh',2.5,'blinkElecs',[],...
        'blinkFeedback','off','move','on','moveThresh',2,'moveElecs',[],'moveFeedback','off','muscle','on','muscleThresh',0.5,...
        'muscleFreqLow',[],'muscleFreqHigh',[],'muscleFeedback','off','elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off');

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
    
    %Set electrode defaults
    if isempty(options.blinkElecs)
        options.blinkElecs = {'Fp1','Fp2'};
    end
    if isempty(options.moveElecs)
        options.moveElecs = {'F7','F8'};
    end
    
    %Check for ICA
    if isempty(EEG.icaweights)
        error('There are no stored components in the data. Please run ICA before running this script. See tesa_fastica.')
    end
    
    %Check comps input
    if ~isempty(options.comps)
        if size(EEG.icaweights,1) < options.comps
            error('The number of components to choose from (%d) is more than the number of independent components in the data (%d).',options.comps,size(EEG.icaweights,1));
        end
    end
    
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
    
    %Set muscle windows based on plotFreqX
    if isempty(options.muscleFreqLow)
        options.muscleFreqLow = [options.plotFreqX(1,1),30];
    end
    if isempty(options.muscleFreqHigh)
        options.muscleFreqHigh = [31,options.plotFreqX(1,2)];
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
    varsNameLog = strncmp('icaBadComp',varsName,10);
    if sum(varsNameLog) == 0
        nameIn = 'icaBadComp1';
    else
        nameIn = ['icaBadComp', num2str(sum(varsNameLog)+1)];
    end

    %Ranks components, outputs variance
    [EEG, varsPerc] = tesa_sortcomps(EEG);

    %Calculates time course and variance
    compTimeCourse = arrayfun(@(x)eeg_getdatact(EEG, 'component', [x], 'projchan', []),1:size(EEG.icawinv,2),'UniformOutput', false);

    global redo
    redo = true;

    while redo == true;

        %Creates output storage for components
        EEG.(nameIn).tmsMuscle = [];
        EEG.(nameIn).eye = [];
        EEG.(nameIn).muscle = [];
        EEG.(nameIn).electrodeNoise = [];
        EEG.(nameIn).sensory = [];

        EEG.(nameIn).tmsMuscleVars = [];
        EEG.(nameIn).eyeVars = [];
        EEG.(nameIn).muscleVars = [];
        EEG.(nameIn).electrodeNoiseVars = [];
        EEG.(nameIn).sensoryVars= [];

        %Number of components to consider
        if isempty(options.comps)
            comps = size(EEG.icaweights,1);
        else
            comps = options.comps;
        end

    for compNum =1:comps

        %##########
        %Output the data

        %Calculates time course
        temp = compTimeCourse{1,compNum};
        temp = squeeze(temp(1,:,:));

        %Calculates FFT (uV/Hz)
        T = 1/EEG.srate;             % Sample time
        L = size(EEG.times,2);       % Length of signal
        NFFT = 2^nextpow2(L);        % Next power of 2 from length of y
        f = EEG.srate/2*linspace(0,1,NFFT/2+1); % Frequencies

        y = reshape(temp,1,[]);
        Y = fft(zscore(y),NFFT)/L;

        Yout = (abs(Y).^2); 
        [c index]=min(abs(f-(0.5)));

        for x = 1:size(Yout,1);
            freq=options.plotFreqX(1,1):0.5:options.plotFreqX(1,2);
            for a=1:size(freq,2);
                [c index1]=min(abs(f-((freq(1,a)-0.25))));
                [c index2]=min(abs(f-((freq(1,a)+0.25))));
                Y2(x,a)=mean(Yout(x,index1:index2),2); %creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
            end;
        end;

        %##########
        %Checks which components fit profile of artefact

        %Calculates the zscores across electrodes for each component
        tempCompZ = zscore(EEG.icawinv(:,compNum));

        %TMS-evoked muscle
        if strcmpi(options.tmsMuscle,'on')
            
            %Check TMS muscle inputs
            if options.tmsMuscleThresh < 0
                error('Input for ''tmsMuscleThresh'' must be greater than 0.');
            elseif size(options.tmsMuscleWin,2)~=2
                error('Input for ''tmsMusclesWin'' must be in the following format: [start,end]. e.g. [11,51].');
            elseif size(options.tmsMuscleBase,2)~=2
                error('Input for ''tmsMusclesBase'' must be in the following format: [start,end]. e.g. [51,100].');
            elseif options.tmsMuscleWin(1,1) < EEG.times(1,1) || options.tmsMuscleWin(1,2) > EEG.times(1,end)
                error('Input for ''tmsMuscleWin'' must be within the limits of the data [%d to %d].',EEG.times(1,1),EEG.times(1,end));
            elseif options.tmsMuscleBase(1,1) < EEG.times(1,1) || options.tmsMuscleBase(1,2) > EEG.times(1,end)
                error('Input for ''tmsMuscleBase'' must be within the limits of the data [%d to %d].',EEG.times(1,1),EEG.times(1,end));
            end
        end
                    
        [val1,mt1] = min(abs(EEG.times-options.tmsMuscleWin(1,1)));
        [val2,mt2] = min(abs(EEG.times-options.tmsMuscleWin(1,2)));
        tmsMuscleSD = std(mean(temp(mt1:mt2,:),2));
        [val1,bt1] = min(abs(EEG.times-options.tmsMuscleBase(1,1)));
        [val2,bt2] = min(abs(EEG.times-options.tmsMuscleBase(1,2)));    
        baseSD = std(mean(temp(bt1:bt2,:),2));
        tmsMuscleRatio = tmsMuscleSD./baseSD;
        if strcmpi(options.tmsMuscleFeedback,'on')
            fprintf('Comp. %d TMS-evoked muscle ratio is %s.\n', compNum,num2str(round(tmsMuscleRatio,2)));
        end

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

        blinkRatio = mean(tempCompZ(eNum,:));
        if strcmpi(options.blinkFeedback,'on')
            fprintf('Comp. %d blink ratio is %s.\n', compNum,num2str(round(blinkRatio,2)));
        end

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
        moveVal = tempCompZ(eNum,:);
        if strcmpi(options.moveFeedback,'on')
            fprintf('Comp. %d eye movement values are %s and %s.\n', compNum,num2str(round(moveVal(1,1),2)),num2str(round(moveVal(2,1),2)));
        end

        %Muscle
        if strcmpi(options.muscle,'on')
            
            %Checks input for muscle
            if size(options.muscleFreqLow,2) ~=2
                error('Inputs for ''muscleFreqLow'' must be in the following format: [low, high]. e.g. [1,30].');
            elseif size(options.muscleFreqHigh,2) ~=2
                error('Inputs for ''muscleFreqHigh'' must be in the following format: [low, high]. e.g. [31,100].');
            elseif (options.muscleFreqLow(1,1) < options.plotFreqX(1,1) | options.muscleFreqLow(1,2) > options.plotFreqX(1,2))
                error('Inputs for ''muscleFreqLow'' (%d to %d) are outside of the frequency range set by input ''plotFreqX'' (%d to %d). Please adjust.',options.muscleFreqLow(1,1),options.muscleFreqLow(1,2),options.plotFreqX(1,1),options.plotFreqX(1,2));
            elseif (options.muscleFreqHigh(1,1) < options.plotFreqX(1,1) | options.muscleFreqHigh(1,2) > options.plotFreqX(1,2))
                error('Inputs for ''muscleFreqHigh'' (%d to %d) are outside of the frequency range set by input ''plotFreqX'' (%d to %d). Please adjust.',options.muscleFreqHigh(1,1),options.muscleFreqHigh(1,2),options.plotFreqX(1,1),options.plotFreqX(1,2));
            end
        end
            
        [val1,lowF1] = min(abs(freq-options.muscleFreqLow(1,1)));
        [val2,lowF2] = min(abs(freq-options.muscleFreqLow(1,2)));
        lowFreq = mean(Y2(:,lowF1:lowF2),2);
        [val1,highF1] = min(abs(freq-options.muscleFreqHigh(1,1)));
        [val2,highF2] = min(abs(freq-options.muscleFreqHigh(1,2)));    
        highFreq = mean(Y2(:,highF1:highF2),2);
        muscleRatio = highFreq./lowFreq;
        if strcmpi(options.muscleFeedback,'on')
            fprintf('Comp. %d muscle ratio is %s.\n', compNum,num2str(round(muscleRatio,2)));
        end

        %Electrode noise
           
        elecNoise = abs(tempCompZ) > abs(options.elecNoiseThresh);
        if strcmpi(options.elecNoiseFeedback,'on')
            fprintf('Comp. %d maximum electrode z score is %s.\n', compNum,num2str(round(max(tempCompZ),2)));
        end
        
        %Select if component is artefact
        if strcmpi(options.tmsMuscle,'on') && tmsMuscleRatio >= options.tmsMuscleThresh
            compVal = 2;
        elseif strcmpi(options.blink,'on') && abs(blinkRatio) >= options.blinkThresh
            compVal = 3;
        elseif strcmpi(eyeMove,'on') && strcmpi(options.move,'on') && (moveVal(1,1) >= options.moveThresh && moveVal(2,1) <= -options.moveThresh) | (moveVal(2,1) >= options.moveThresh && moveVal(1,1) <= -options.moveThresh)
            compVal = 3;    
        elseif strcmpi(options.muscle,'on') && abs(muscleRatio) >= options.muscleThresh
            compVal = 4;
        elseif strcmpi(options.elecNoise,'on') && sum(elecNoise) >= 1
            compVal = 5;
        else
            compVal = 1;
        end
        
        %##########
        %Plots the figure

        %Decide colour of figure
        if compVal ~= 1
            colour = 'r'; % if artefact
        else
            colour = 'b'; % if not artefac
        end

        %Figure sizing
        if strcmpi(options.figSize,'small')
            sz = [560, 420]; % figure size
            popPos = [345, 185, 140, 50];
            popFont = 9;
            compPos = [115, 400, 100, 20];
            compFont = 12;
            varPos =  [340, 397, 150, 20];
            varFont = 9;
        elseif strcmpi(options.figSize,'medium')
            sz = [900, 600]; % figure size
            popPos = [560, 260, 210, 75];
            popFont = 14;
            compPos = [190, 565, 150, 30];
            compFont = 18;
            varPos =  [550, 560, 240, 30];
            varFont = 14;
        elseif strcmpi(options.figSize,'large')
            sz = [1200, 900]; % figure size
            popPos = [750, 400, 280, 100];
            popFont = 18;
            compPos = [250, 855, 200, 40];
            compFont = 24;
            varPos =  [740, 850, 300, 40];
            varFont = 18;
        end

        screensize = get(0,'ScreenSize');
        xpos = ceil((screensize(3)-sz(1))/2); % center figure horizontally
        ypos = ceil((screensize(4)-sz(2))/2); % center figure vertically

        f = figure('KeyPressFcn',@(obj,evt) 0);
        f.Position = [xpos ypos sz(1) sz(2)];

        %Plot time course
        subplot(2,2,1);
        plot(EEG.times,mean(temp,2),colour); grid on; hold on;
        plot([0 0], get(gca,'ylim'),'r--');
        set(gca,'Xlim', [options.plotTimeX(1,1), options.plotTimeX(1,2)]);
        xlabel('Time (ms)');
        ylabel('Amplitude (a.u.)');

        %Plot topoplot
        subplot(2,2,2);
        topoplot(EEG.icawinv(:,compNum),EEG.chanlocs,'electrodes','off');

        %Plot time course matrix
        [val1,tp1] = min(abs(EEG.times-options.plotTimeX(1,1)));
        [val2,tp2] = min(abs(EEG.times-options.plotTimeX(1,2)));
        temp1 = temp(tp1:tp2,:);
        subplot(2,2,3);
        imagesc(temp1','XData', options.plotTimeX);
        caxis([-max(abs(temp1(:))), max(abs(temp1(:)))]);
        xlabel('Time (ms)');
        ylabel('Trials');

        subplot(2,2,4);
        plot(freq,Y2,colour);grid on;
        set(gca,'Xlim', options.plotFreqX);
        xlabel('Time (ms)');
        ylabel('Power (\muV/Hz)');

        %Plot popup window
        popup = uicontrol('Style', 'popup',...
            'String', {'Neural','TMS-evoked muscle','Eye','Muscle','Electrode noise','Sensory'},...
            'Position', popPos,...
            'Value',compVal,...
            'fontSize',popFont); 

        %Plot component number
        hT = uicontrol('style', 'text',... 
            'string', ['IC ', num2str(compNum), ' of ', num2str(size(EEG.icaweights,1))],... 
            'position', compPos,...
            'BackgroundColor',f.Color,...
            'fontWeight','bold',...
            'fontSize',compFont);

        %Variance info
        hT2 = uicontrol('style', 'text',... 
            'string', ['Var. accounted for: ', num2str(round(varsPerc(1,compNum),1)),' %'],... 
            'position', varPos,...
            'BackgroundColor',f.Color,...
            'fontSize',varFont);

        waitfor(gcf,'CurrentCharacter');
        curChar=uint8(get(gcf,'CurrentCharacter'));

        output = popup.Value;

        %Save components for removal
        if output == 2
            EEG.(nameIn).tmsMuscle(1,size(EEG.(nameIn).tmsMuscle,2)+1) = compNum;
            EEG.(nameIn).tmsMuscleVars(1,size(EEG.(nameIn).tmsMuscleVars,2)+1) = varsPerc(1,compNum);
        elseif output == 3
            EEG.(nameIn).eye(1,size(EEG.(nameIn).eye,2)+1) = compNum;
            EEG.(nameIn).eyeVars(1,size(EEG.(nameIn).eyeVars,2)+1) = varsPerc(1,compNum);
        elseif output == 4
            EEG.(nameIn).muscle(1,size(EEG.(nameIn).muscle,2)+1) = compNum;
            EEG.(nameIn).muscleVars(1,size(EEG.(nameIn).muscleVars,2)+1) = varsPerc(1,compNum);
        elseif output == 5
            EEG.(nameIn).electrodeNoise(1,size(EEG.(nameIn).electrodeNoise,2)+1) = compNum;
            EEG.(nameIn).electrodeNoiseVars(1,size(EEG.(nameIn).electrodeNoiseVars,2)+1) = varsPerc(1,compNum);
        elseif output == 6
            EEG.(nameIn).sensory(1,size(EEG.(nameIn).sensory,2)+1) = compNum;
            EEG.(nameIn).sensoryVars(1,size(EEG.(nameIn).sensoryVars,2)+1) = varsPerc(1,compNum);
        end

        close;
    end

    pre = mean(EEG.data,3);
    EEG1 = EEG;

    %Remove bad components
    removeComps = [EEG.(nameIn).tmsMuscle, EEG.(nameIn).eye, EEG.(nameIn).muscle, EEG.(nameIn).electrodeNoise, EEG.(nameIn).sensory];
    EEG = pop_subcomp( EEG,removeComps, 0);

    %Plot check
    f1 = figure;
    subplot(1,2,1)
    plot(EEG.times,pre,'b'); grid on; hold on;
    plot([0 0], get(gca,'ylim'),'r--');
    set(gca,'Xlim', options.plotTimeX);
    xlabel('Time (ms)');
    ylabel('Amplitude (\muV)');
    title('Before correction','fontweight','bold');
    yScale = get(gca,'ylim');

    subplot(1,2,2)
    plot(EEG.times,mean(EEG.data,3),'b'); grid on; hold on;
    plot([0 0], yScale,'r--');
    set(gca,'Xlim', options.plotTimeX,'ylim',yScale);
    xlabel('Time (ms)');
    ylabel('Amplitude (\muV)');
    title('After correction','fontweight','bold');

    %Check if happy to continue
    ButtonName = questdlg('Are you happy with component removal?', 'Verify component removal', 'No-Redo', 'Yes', 'Cancel', 'Yes');
    switch ButtonName,
        case 'No-Redo',
            redo = true; EEG = EEG1;
        case 'Yes',
            redo = false; 
        case 'Cancel',
            redo = false; error('Script termninated by user.');
    end 

    close;

    end
end
