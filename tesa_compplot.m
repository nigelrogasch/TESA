% tesa_compplot() - Plots independent components following ICA and allows users to 
%                   classify components for removal using a dropdown menu.
%                   Each plot includes a butterfly plot of the EEG data with
%                   and without the current IC, the average time course of
%                   the component, the trial-by-trial amplitude of the
%                   component, the topography of the component, and the
%                   power spectrum of the component. Infomration on the %
%                   mean variance (of total variance) accounted for by each
%                   component is also provided. Users can interactively
%                   alter the x and y axes of the plots by selecting
%                   'Update plots'. Note that three figures sizes are
%                   available: 'small','medium' and 'large'.
%                   When satisfied with classification, press the 'Next' button to 
%                   continue to the next component. Users can scroll backwards and forward
%                   through components by selecting the 'Next IC' and then pressing
%                   'Next'. Users can also skip to the 'Finish' from this dropdown.
%                   Components are classified as one of the following: 
%                   keep, reject, reject - TMS-evoked muscle, reject - blink, 
%                   reject - eye movement, reject - muscle, reject - 
%                   electrode noise, or reject - sensory . If a component classification
%                   algorithm has been used to classify components in TESA, this
%                   can be used as an input to tesa_compplot to automatically
%                   mark components for rejection. If a component is selected
%                   for rejection, it is plotted as red on subsequent plots, if
%                   not it is plotted as blue. Following removal, a before and
%                   after plot is generated and users can select whether to
%                   continue or repeat the process. Component numbers and
%                   variance represented by each component selected for
%                   removal is stored in the EEG structure =>
%                   EEG.icaCompClass.xxx
% 
%                   Note that tesa_compplot is automatically called by
%                   tesa_compselect if 'compCheck' is 'on'.
% 
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
%                   Rogasch NC et al. (2017) Analysing concurrent transcranial magnetic stimulation and
%                   electroencephalographic data: A review and introduction to the open-source
%                   TESA software NeuroImage, 147:934–951               
%
% Usage:
%   >>  EEG = tesa_compplot( EEG );
%   >>  EEG = tesa_compplot( EEG, 'key1',value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs (varargin):
%   'compClassInput',EEGLAB structure field - location in EEGLAB structure
%                   where output from TESA components classification
%                   algorithms are stored (e.g. tesa_compselect). Note that the number
%                   of components that have been classified must match the
%                   the number of components currently stored in the data
%                   (i.e. this won't work if you have already removed the
%                   components in tesa_compselect).
%                   Example: EEG.icaCompClass.TESA1
%                   Default: empty
%   'saveWeights','str' - 'on' | 'off'. Saves the winv (topoplot weights)
%                   and icaact (time course weights) that the
%                   classification was performed on. Note this will
%                   increase the storage size of the files.
%                   Default: 'off'
%   'figSize','str' - 'small' | 'medium' | 'large'; Determines the size of
%                   the figures that are plotted displaying the information
%                   on the components.
%                   Default: 'medium'
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
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples
%   EEG = tesa_compplot( EEG ); % default use
%   EEG = tesa_compplot( EEG, 'compClassInput', EEG.icaCompClass.TESA1 ); % Loads automated component classifications performed by tesa_compselect
%   EEG = tesa_compplot( EEG, 'figSize', 'small','plotTimeX',[-600,600], 'plotFreqX', [2,45], 'freqScale', 'db' ); % changes the number of components to select, the size of the figure and the x axis' of the time course and frequency plots, and the scale of the y-axis on the frequency plots.
% 
% See also:
%   tesa_compselect, tesa_sortcomps 

% Copyright (C) 2018  Nigel Rogasch, Monash University,
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
% NR 16-6-2020
% - Changed the position of the 'Next' and 'Back' callback functions to
% inside the main function to avoid use of global variables
% - Included a 'check' variable which disables Next/Back buttons while the
% main function is updating after a button press and avoids the need for a
% delay in the callback functions.
% - Included a random list of name variants in Next/Back callback functions
% to prevent function from stalling. The waitfor function was occasionally
% missing the change in figure name which prevented the Next/Back buttons
% from working on subsequent presses, giving the impression that the
% function had frozen. By randomising the names which change on each press
% (i.e. by including different amounts of spacing in string), this issue is
% circumvented.

function EEG = tesa_compplot(EEG,varargin)

if nargin < 1
    error('Not enough input arguments.');
end

%define options
options = struct('compClassInput',[],'saveWeights','off','figSize',[],'plotTimeX',[],'plotFreqX',[],'varsPerc',[],'fftBins',[],'freqBins',[],'freqScale',[]);

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('tesa_compplot needs key/value pairs')
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
if isempty(options.compClassInput)
    options.compClassInput = [];
end
if isempty(options.figSize)
    options.figSize = 'medium';
end
if isempty(options.plotTimeX)
    options.plotTimeX = [-200,500];
end
if isempty(options.plotFreqX)
    options.plotFreqX = [1,100];
end
if isempty(options.varsPerc)
    varsPerc = [];
else
    varsPerc = options.varsPerc;
end
if isempty(options.fftBins)
    fftBins = [];
else
    fftBins = options.fftBins;
end
if isempty(options.freqBins)
    freq = [];
else
    freq = options.freqBins;
end
if isempty(options.freqScale)
    options.freqScale = 'log';
end

%Check for ICA
if isempty(EEG.icaweights)
    error('There are no stored components in the data. Please run ICA before running this script. See tesa_fastica.')
end

%Check for channel locations
if isempty(EEG.chanlocs(1).theta)
    error('Channel locations are required for this function. Please load channel locations: Edit -> Channel locations.')
end

%Check figure inputs
if ~(strcmp(options.figSize,'small') | strcmp(options.figSize,'medium') | strcmp(options.figSize,'large'))
    error ('Input for ''figSize'' needs to be either ''small'', ''medium'' or ''large''.');
end

%Check plotTimeX inputs
if size(options.plotTimeX,2) ~= 2
    error('Input for ''plotTimeX'' must be in the following format: [start, end]. e.g. [-200,500].');
elseif options.plotTimeX(1,1) < EEG.times(1,1) || options.plotTimeX(1,2) > EEG.times(1,end)
    error('Input for ''plotTimeX'' must be within the limits of the data [%d to %d].',round(EEG.times(1,1)),round(EEG.times(1,end)));
end

%Check plotFreqX inputs
if size(options.plotFreqX,2) ~= 2
    error('Input for ''plotFreqX'' must be in the following format: [low, high]. e.g. [1,100].');
elseif options.plotFreqX(1,1)<0
    error('Input for ''plotFreqX'' must be larger than 0.');
end

%Check frequency scaling inputs
if ~(strcmp(options.freqScale,'raw') | strcmp(options.freqScale,'log') | strcmp(options.freqScale,'log10') | strcmp(options.freqScale,'db'))
    error ('Input for ''freqScale'' needs to be either ''raw'', ''log'', ''log10'' or ''db''.');
end
    
%Calculate component time series
if isempty(EEG.icaact)
    fprintf('Calculating component time series\n')
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
end

%Calculate component variance
if isempty(varsPerc)
    fprintf('Calculating component variance\n')
    for x = 1:size(EEG.icaact,1)
        vars(x) = var(mean(EEG.icaact(x,:,:),3));
    end
    varsPerc = vars/sum(vars)*100;
end

% Number of components
compNumMat = 1:size(EEG.icawinv,2);

%Calculates component power spectrum
if isempty(fftBins)
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
end

% Check that both the binned fft outputs and the corresponding frequency
% bins are present
if ~isempty(fftBins)
    if isempty(freq)
        error('Please provide frequency bins for fftBins input with: ''freqBins'', freq - where freq is a 1xn matrix of same size as fftBins');
    end
end

% Generate cell array with number of components + finish
compNumCell = [];
for iComp = 1:length(compNumMat)
    compNumCell{iComp} = num2str(compNumMat(iComp));
end
compNumCell{end+1} = 'Finish';

% Generate array for component classification: 1 = keep; >1 = reject
if isempty(options.compClassInput)
    compClass = ones(1,size(EEG.icawinv,2));
else
    compClass = options.compClassInput.compClass;
end

% Check the compClass input against the stored ICs
if length(compClass) ~= size(EEG.icaweights,1)
    error('The number of components in the classification array (%d) does not match the number of components in the data (%d)\nMake sure the correct component classification was loaded for this data',length(compClass), size(EEG.icaweights,1));
end

% Set up figure inputs
global nextPlot
nextPlot = true;
nComp = 1;

global backPlot
backPlot = 0;

% Launch the figure
fprintf('Launching figure\n')

%##########
%Plots the figure

%Figure sizing
if strcmpi(options.figSize,'small')
    sz = [560, 420]; % figure size
    sf = 0.7;
    compFont = 12; % IC title font
    popFont = 9; % Other title fonts
    varFont = 9; % variance font
    compPos = [sz(1,1)*0.72, sz(1,2)*0.93, 150*sf, 30*sf]; % IC title position
    varPos =  [sz(1,1)*0.64, sz(1,2)*0.88, 300*sf, 30*sf]; % variance position
    figOptPos =  [sz(1,1)*0.64, sz(1,2)*0.82, 300*sf, 30*sf]; % figOpt position
    timePos = [sz(1,1)*0.8, sz(1,2)*0.77, 150*sf, 30*sf]; % time position
    timeTitPos = [sz(1,1)*0.63, sz(1,2)*0.77, 130*sf, 30*sf]; % time title position
    freqPos = [sz(1,1)*0.8, sz(1,2)*0.72, 150*sf, 30*sf]; % freq position
    freqTitPos = [sz(1,1)*0.63, sz(1,2)*0.72, 130*sf, 30*sf]; % freq title position
    ampPos = [sz(1,1)*0.8, sz(1,2)*0.67, 150*sf, 30*sf]; % amp position
    ampTitPos = [sz(1,1)*0.63, sz(1,2)*0.67, 130*sf, 30*sf]; % amp title position
    auPos = [sz(1,1)*0.8, sz(1,2)*0.62, 150*sf, 30*sf]; % au position
    auTitPos = [sz(1,1)*0.63, sz(1,2)*0.62, 130*sf, 30*sf]; % au title position
    figButPos = [sz(1,1)*0.7, sz(1,2)*0.53, 200*sf, 40*sf]; % update figure buttonposition
    scrollTitPos = [sz(1,1)*0.68, sz(1,2)*0.45, 240*sf, 30*sf]; % scrolltitle position
    scrollButPos = [sz(1,1)*0.7, sz(1,2)*0.38, 200*sf, 40*sf]; % scroll buttonposition
    catPos = [sz(1,1)*0.67, sz(1,2)*0.3, 240*sf, 30*sf]; % category title position
    popPos = [sz(1,1)*0.67, sz(1,2)*0.25, 250*sf, 30*sf]; % category dropdown position
    nextPos = [sz(1,1)*0.7, sz(1,2)*0.14, 80*sf, 30*sf]; % next comp position
    popPosComp = [sz(1,1)*0.81, sz(1,2)*0.15, 80*sf, 30*sf]; % next comp position
    NextButPos = [sz(1,1)*0.83, sz(1,2)*0.05, 120*sf, 50*sf]; % next comp buttonposition
    BackButPos = [sz(1,1)*0.67, sz(1,2)*0.05, 120*sf, 50*sf]; % back comp buttonposition
elseif strcmpi(options.figSize,'medium')
    sz = [900, 600]; % figure size       
    sf = 1;
    compFont = 18; % IC title font
    popFont = 14; % Other title fonts
    varFont = 14; % variance font
    compPos = [sz(1,1)*0.72, sz(1,2)*0.93, 150*sf, 30*sf]; % IC title position
    varPos =  [sz(1,1)*0.64, sz(1,2)*0.88, 300*sf, 30*sf]; % variance position
    figOptPos =  [sz(1,1)*0.64, sz(1,2)*0.82, 300*sf, 30*sf]; % figOpt position
    timePos = [sz(1,1)*0.8, sz(1,2)*0.77, 150*sf, 30*sf]; % time position
    timeTitPos = [sz(1,1)*0.63, sz(1,2)*0.77, 150*sf, 30*sf]; % time title position
    freqPos = [sz(1,1)*0.8, sz(1,2)*0.72, 150*sf, 30*sf]; % freq position
    freqTitPos = [sz(1,1)*0.63, sz(1,2)*0.72, 150*sf, 30*sf]; % freq title position
    ampPos = [sz(1,1)*0.8, sz(1,2)*0.67, 150*sf, 30*sf]; % amp position
    ampTitPos = [sz(1,1)*0.63, sz(1,2)*0.67, 150*sf, 30*sf]; % amp title position
    auPos = [sz(1,1)*0.8, sz(1,2)*0.62, 150*sf, 30*sf]; % au position
    auTitPos = [sz(1,1)*0.63, sz(1,2)*0.62, 150*sf, 30*sf]; % au title position
    figButPos = [sz(1,1)*0.7, sz(1,2)*0.53, 200*sf, 40*sf]; % update figure buttonposition
    scrollTitPos = [sz(1,1)*0.68, sz(1,2)*0.45, 240*sf, 30*sf]; % scrolltitle position
    scrollButPos = [sz(1,1)*0.7, sz(1,2)*0.38, 200*sf, 40*sf]; % scroll buttonposition
    catPos = [sz(1,1)*0.67, sz(1,2)*0.3, 240*sf, 30*sf]; % category title position
    popPos = [sz(1,1)*0.67, sz(1,2)*0.25, 250*sf, 30*sf]; % category dropdown position
    nextPos = [sz(1,1)*0.7, sz(1,2)*0.14, 80*sf, 30*sf]; % next comp position
    popPosComp = [sz(1,1)*0.81, sz(1,2)*0.15, 80*sf, 30*sf]; % next comp position
    NextButPos = [sz(1,1)*0.83, sz(1,2)*0.05, 120*sf, 50*sf]; % next comp buttonposition
    BackButPos = [sz(1,1)*0.67, sz(1,2)*0.05, 120*sf, 50*sf]; % back comp buttonposition
elseif strcmpi(options.figSize,'large')
    sz = [1200, 900]; % figure size
    sf = 1.5;
    compFont = 20; % IC title font
    popFont = 18; % Other title fonts
    varFont = 18; % variance font
    compPos = [sz(1,1)*0.72, sz(1,2)*0.93, 150*sf, 30*sf]; % IC title position
    varPos =  [sz(1,1)*0.64, sz(1,2)*0.88, 300*sf, 30*sf]; % variance position
    figOptPos =  [sz(1,1)*0.64, sz(1,2)*0.82, 300*sf, 30*sf]; % figOpt position
    timePos = [sz(1,1)*0.8, sz(1,2)*0.77, 150*sf, 30*sf]; % time position
    timeTitPos = [sz(1,1)*0.63, sz(1,2)*0.77, 150*1.3, 30*sf]; % time title position
    freqPos = [sz(1,1)*0.8, sz(1,2)*0.72, 150*sf, 30*sf]; % freq position
    freqTitPos = [sz(1,1)*0.63, sz(1,2)*0.72, 150*1.3, 30*sf]; % freq title position
    ampPos = [sz(1,1)*0.8, sz(1,2)*0.67, 150*sf, 30*sf]; % amp position
    ampTitPos = [sz(1,1)*0.63, sz(1,2)*0.67, 150*1.3, 30*sf]; % amp title position
    auPos = [sz(1,1)*0.8, sz(1,2)*0.62, 150*sf, 30*sf]; % au position
    auTitPos = [sz(1,1)*0.63, sz(1,2)*0.62, 150*1.3, 30*sf]; % au title position
    figButPos = [sz(1,1)*0.7, sz(1,2)*0.53, 200*sf, 40*sf]; % update figure buttonposition
    scrollTitPos = [sz(1,1)*0.68, sz(1,2)*0.45, 240*sf, 30*sf]; % scrolltitle position
    scrollButPos = [sz(1,1)*0.7, sz(1,2)*0.38, 200*sf, 40*sf]; % scroll buttonposition
    catPos = [sz(1,1)*0.67, sz(1,2)*0.3, 240*sf, 30*sf]; % category title position
    popPos = [sz(1,1)*0.67, sz(1,2)*0.25, 250*sf, 30*sf]; % category dropdown position
    nextPos = [sz(1,1)*0.7, sz(1,2)*0.14, 80*sf, 30*sf]; % next comp position
    popPosComp = [sz(1,1)*0.81, sz(1,2)*0.15, 80*sf, 30*sf]; % next comp position
    NextButPos = [sz(1,1)*0.83, sz(1,2)*0.05, 120*sf, 50*sf]; % next comp buttonposition
    BackButPos = [sz(1,1)*0.67, sz(1,2)*0.05, 120*sf, 50*sf]; % back comp buttonposition
end

screensize = get(0,'ScreenSize');
xpos = ceil((screensize(3)-sz(1))/2); % center figure horizontally
ypos = ceil((screensize(4)-sz(2))/2); % center figure vertically

f = figure('KeyPressFcn',@(obj,evt) 0);
f.Position = [xpos ypos sz(1) sz(2)];
f.Name = 'Press next when selection is made.';
f.NumberTitle = 'off';

while nextPlot == true
    
    if exist('f')==0
        f = figure('KeyPressFcn',@(obj,evt) 0);
        f.Position = [xpos ypos sz(1) sz(2)];
        f.Name = 'Press next when selection is made.';
        f.NumberTitle = 'off';
    end
    
    % Determine what colour to make the plots based on component category
    if compClass(nComp) == 1
        colour = 'b';
    else
        colour = 'r';
    end
    
    % Set up component removal plot
    % Calculate manually on average (faster than calculating on individual
    % trials)
    % Time series with component in
    allComps = 1:length(compClass);
    tempCompGood = compClass == 1;
    tempCompGood(nComp) = 1;
    goodC = allComps(tempCompGood);
    dataIn = mean(EEG.data,3); % Calculate on average data - faster
    S = (EEG.icaweights(goodC,:)*EEG.icasphere)*dataIn(EEG.icachansind,:);
    Winv = EEG.icawinv(:, goodC);
    EEGdataIn = Winv*S;

    % Time series with component out
    tempCompGood = compClass == 1;
    tempCompGood(nComp) = 0;
    goodC = allComps(tempCompGood);
    dataIn = mean(EEG.data,3); % Calculate on average data - faster
    S = (EEG.icaweights(goodC,:)*EEG.icasphere)*dataIn(EEG.icachansind,:);
    Winv = EEG.icawinv(:, goodC);
    EEGdataOut = Winv*S;
    
    %Plot with component in
    P.p1 = subplot(3,3,1);
    plot(EEG.times,EEGdataIn,'b'); grid on; hold on;
    tempP1 = P.p1.YLim;
    plot([0 0], [-1000,1000],'r--'); hold off;
    set(gca,'Xlim', [options.plotTimeX(1,1), options.plotTimeX(1,2)]);
    xlabel('Time (ms)');
    ylabel('Amplitude (\muV)');
    title('With IC');
    
    %Plot with component out
    P.p2 = subplot(3,3,2);
    plot(EEG.times,EEGdataOut,'b'); grid on; hold on;
    tempP2 = P.p2.YLim;
    plot([0 0], [-1000,1000],'r--'); hold off;
    set(gca,'Xlim', [options.plotTimeX(1,1), options.plotTimeX(1,2)]);
    xlabel('Time (ms)');
    ylabel('Amplitude (\muV)');
    title('Without IC');
    
    % Y limits
    negMax = max([abs(tempP1(1,1)),abs(tempP2(1,1))]);
    posMax = max([abs(tempP1(1,2)),abs(tempP2(1,2))]);
    P.p1.YLim = [-negMax,posMax];
    P.p2.YLim = [-negMax,posMax];
    
    %Plot time course
    P.p3 = subplot(3,3,4);
    plot(EEG.times,mean(EEG.icaact(nComp,:,:),3),colour); grid on; hold on;
    plot([0 0], get(gca,'ylim'),'r--'); hold off;
    set(gca,'Xlim', [options.plotTimeX(1,1), options.plotTimeX(1,2)]);
    xlabel('Time (ms)');
    ylabel('Amplitude (a.u.)');
    title('IC time series');
    
    %Plot topoplot
    P.p4 = subplot(3,3,5);
    cla(P.p4)
    topoplot(EEG.icawinv(:,nComp),EEG.chanlocs,'electrodes','off');hold off;
%     colorbar;
    title('IC topography');
    
    %Plot time course matrix
    [val1,tp1] = min(abs(EEG.times-options.plotTimeX(1,1)));
    [val2,tp2] = min(abs(EEG.times-options.plotTimeX(1,2)));
    temp1 = squeeze(EEG.icaact(nComp,1:end,:));
    P.p5 = subplot(3,3,7);
    imagesc(temp1','XData', [EEG.times(1),EEG.times(end)]);
    P.p5.XLim = [options.plotTimeX(1,1),options.plotTimeX(1,2)];
    caxis([-max(abs(temp1(:))), max(abs(temp1(:)))]);
    xlabel('Time (ms)');
    ylabel('Trials');
    title('IC time series (trials)');
    
    % Plot frequency spectra
    P.p6 = subplot(3,3,8);
    if strcmp(options.freqScale,'raw')
        plot(freq,fftBins(nComp,:),colour);grid on;
        ylabel('Power (\muV^{2}/Hz)');
    elseif strcmp(options.freqScale,'log')
        plot(freq,log(fftBins(nComp,:)),colour);grid on;
        ylabel('Power log(\muV^{2}/Hz)');
    elseif strcmp(options.freqScale,'log10')
        plot(freq,log10(fftBins(nComp,:)),colour);grid on;
        ylabel('Power log10(\muV^{2}/Hz)');
    elseif strcmp(options.freqScale,'db')
        plot(freq,10*log10(fftBins(nComp,:)),colour);grid on;
        ylabel('Power dB 10*log10(\muV^{2}/Hz)');
    end
    set(gca,'Xlim', options.plotFreqX);
    xlabel('Frequency (Hz)');
    title('IC power spectrum');
       
    %Plot popup window for categorisation
    popup = uicontrol('Style', 'popup',...
        'String', {'Keep','Reject','Reject - TMS-evoked muscle','Reject - Blink','Reject - Eye-movement','Reject - Muscle','Reject - Electrode noise','Reject - Sensory'},...
        'Position', popPos,...
        'Value',compClass(nComp),...
        'fontSize',popFont);
    
%       'Keep' = 1
%       'Reject' = 2
%       'Reject - TMS-evoked muscle' = 3
%       'Reject - Blink' = 4 
%       'Reject - Eye-movement' = 5
%       'Reject - Muscle' = 6
%       'Reject - Electrode noise' = 7
%       'Reject - Sensory' = 8
    
    %Plot popup window from component control
    popupComp = uicontrol('Style', 'popup',...
        'String', compNumCell,...
        'Position', popPosComp,...
        'Value',nComp+1,...
        'fontSize',popFont);
    
    %Plot next button
    nextButton = uicontrol('Style', 'pushbutton',...
        'String', {'Next'},...
        'Position', NextButPos,...
        'fontSize',popFont,...
        'Callback', {@nextCommand});
    
    %Plot back button
    backButton = uicontrol('Style', 'pushbutton',...
        'String', {'Back'},...
        'Position', BackButPos,...
        'fontSize',popFont,...
        'Callback', {@backCommand});
    
    %Plot scroll IC button
    scrollButton = uicontrol('Style', 'pushbutton',...
        'String', {'Scroll ICs'},...
        'Position', scrollButPos,...
        'fontSize',popFont,...
        'Callback', 'pop_eegplot( EEG, 0, 1, 0);');
    
    %Plot edit window for time axis
    U.timeEdit = uicontrol('Style', 'edit',...
        'String', {[num2str(options.plotTimeX(1,1)),' ',num2str(options.plotTimeX(1,2))]},...
        'Position', timePos,...
        'fontSize',popFont);
    
    %Plot edit window for freq axis
    U.freqEdit = uicontrol('Style', 'edit',...
        'String', {[num2str(options.plotFreqX(1,1)),' ',num2str(options.plotFreqX(1,2))]},...
        'Position', freqPos,...
        'fontSize',popFont);
    
    %Plot edit window for amp axis
    U.ampEdit = uicontrol('Style', 'edit',...
        'String', {[num2str(P.p1.YLim(1,1)),' ',num2str(P.p1.YLim(1,2))]},...
        'Position', ampPos,...
        'fontSize',popFont);
    
    %Plot edit window for amp axis
    U.auEdit = uicontrol('Style', 'edit',...
        'String', {[num2str(P.p3.YLim(1,1)),' ',num2str(P.p3.YLim(1,2))]},...
        'Position', auPos,...
        'fontSize',popFont);
       
    %Plot update figure button
    figButton = uicontrol('Style', 'pushbutton',...
        'String', {'Update figures'},...
        'Position', figButPos,...
        'fontSize',popFont,...
        'Callback', {@updateFig,U,P});
    
    %Plot component number
    hT = uicontrol('style', 'text',...
        'string', ['IC ', num2str(nComp), ' of ', num2str(size(EEG.icaweights,1))],...
        'position', compPos,...
        'BackgroundColor',f.Color,...
        'fontWeight','bold',...
        'fontSize',compFont);
    
    %Variance info
    hT2 = uicontrol('style', 'text',...
        'string', ['Var. accounted for: ', num2str(round(varsPerc(1,nComp),1)),' %'],...
        'position', varPos,...
        'BackgroundColor',f.Color,...
        'fontSize',varFont);
    
    %Category title
    hT3 = uicontrol('style', 'text',...
        'string', ['IC classification'],...
        'position', catPos,...
        'BackgroundColor',f.Color,...
        'fontWeight','bold',...
        'fontSize',varFont);
    
    %Next component title
    hT4 = uicontrol('style', 'text',...
        'string', ['Next IC:'],...
        'position', nextPos,...
        'BackgroundColor',f.Color,...
        'fontWeight','bold',...
        'fontSize',varFont);
    
    %Figure options title
    hT5 = uicontrol('style', 'text',...
        'string', ['Figure options'],...
        'position', figOptPos,...
        'BackgroundColor',f.Color,...
        'fontWeight','bold',...
        'fontSize',varFont);
    
    %Time title
    hT6 = uicontrol('style', 'text',...
        'string', ['Time (ms)'],...
        'position', timeTitPos,...
        'BackgroundColor',f.Color,...
        'fontSize',varFont);
    
    %Freq title
    hT7 = uicontrol('style', 'text',...
        'string', ['Freq. (Hz)'],...
        'position', freqTitPos,...
        'BackgroundColor',f.Color,...
        'fontSize',varFont);
    
    %Amp title
    hT8 = uicontrol('style', 'text',...
        'string', ['Amp. (uV)'],...
        'position', ampTitPos,...
        'BackgroundColor',f.Color,...
        'fontSize',varFont);
    
    %AU title
    hT9 = uicontrol('style', 'text',...
        'string', ['Amp. (a.u.)'],...
        'position', auTitPos,...
        'BackgroundColor',f.Color,...
        'fontSize',varFont);
    
    %Scroll IC title
    hT9 = uicontrol('style', 'text',...
        'string', ['Scroll ICs by trial'],...
        'position', scrollTitPos,...
        'BackgroundColor',f.Color,...
        'fontWeight','bold',...
        'fontSize',varFont);
       
    check = [];
    drawnow;
    waitfor(f, 'Name');
    drawnow;
    f.Name = 'Press next when selection is made.';
    check = 1;
    
    % Update categorisation of component
    compClass(nComp) = popup.Value;
    
    % Update next component to plot
    if backPlot == 0
        nComp = popupComp.Value;
    elseif backPlot == 1
        if nComp == 1
            nComp = 1;
        else
            nComp = nComp-1;
        end
        backPlot = 0;
    end
        
    if strcmp(compNumCell(nComp),'Finish')
        
        % close figure
        close;
        clear f;
        
        %Figure sizing
        if strcmpi(options.figSize,'small')
            sz = [560, 420]; % figure size
            sf = 0.7;
            compFont = 12; % IC title font
            popFont = 8; % Other title fonts
            varFont = 14; % variance font
            yesButPos = [sz(1,1)*0.1, sz(1,2)*0.1, 200*sf, 75*sf]; % IC title position
            noButPos =  [sz(1,1)*0.4, sz(1,2)*0.1, 200*sf, 75*sf]; % variance position
            cancelButPos =  [sz(1,1)*0.7, sz(1,2)*0.1, 200*sf, 75*sf]; % figOpt position
            butTitPos = [sz(1,1)*0.2, sz(1,2)*0.25, 500*sf, 40*sf];
            compRemTitPos = [sz(1,1)*0.2, sz(1,2)*0.4, 500*sf, 40*sf];
            varRemTitPos = [sz(1,1)*0.2, sz(1,2)*0.35, 500*sf, 40*sf];
        elseif strcmpi(options.figSize,'medium')
            sz = [900, 600]; % figure size   
            sf = 1;
            compFont = 18; % IC title font
            popFont = 12; % Other title fonts
            varFont = 14; % variance font
            yesButPos = [sz(1,1)*0.1, sz(1,2)*0.1, 200*sf, 75*sf]; % IC title position
            noButPos =  [sz(1,1)*0.4, sz(1,2)*0.1, 200*sf, 75*sf]; % variance position
            cancelButPos =  [sz(1,1)*0.7, sz(1,2)*0.1, 200*sf, 75*sf]; % figOpt position
            butTitPos = [sz(1,1)*0.25, sz(1,2)*0.25, 500*sf, 40*sf];
            compRemTitPos = [sz(1,1)*0.25, sz(1,2)*0.4, 500*sf, 40*sf];
            varRemTitPos = [sz(1,1)*0.25, sz(1,2)*0.35, 500*sf, 40*sf];
        elseif strcmpi(options.figSize,'large')
            sz = [1200, 900]; % figure size
            sf = 1.5;
            compFont = 25; % IC title font
            popFont = 16; % Other title fonts
            varFont = 14; % variance font
            yesButPos = [sz(1,1)*0.1, sz(1,2)*0.1, 200*sf, 75*sf]; % IC title position
            noButPos =  [sz(1,1)*0.4, sz(1,2)*0.1, 200*sf, 75*sf]; % variance position
            cancelButPos =  [sz(1,1)*0.7, sz(1,2)*0.1, 200*sf, 75*sf]; % figOpt position
            butTitPos = [sz(1,1)*0.2, sz(1,2)*0.25, 500*sf, 40*sf];
            compRemTitPos = [sz(1,1)*0.2, sz(1,2)*0.4, 500*sf, 40*sf];
            varRemTitPos = [sz(1,1)*0.2, sz(1,2)*0.35, 500*sf, 40*sf];
        end
        
        %Plot check
        global exitScript
        exitScript = 'no';
        f1 = figure('KeyPressFcn',@(obj,evt) 0);
        f1.Name = 'Make response when ready to continue.';
        f1.Position = [xpos ypos sz(1) sz(2)];
        f1.NumberTitle = 'off';
        
        % Plot output with no ICs removed
        allComps = 1:length(compClass);
        tempCompGood = ones(1,length(compClass)) == 1;
        goodC = allComps(tempCompGood);
        dataIn = mean(EEG.data,3); % Calculate on average data - faster
        S = (EEG.icaweights(goodC,:)*EEG.icasphere)*dataIn(EEG.icachansind,:);
        Winv = EEG.icawinv(:, goodC);
        EEGdataIn = Winv*S;

        % Plot output with selected ICs removed
        tempCompGood = compClass == 1;
        tempCompBad = compClass > 1;
        goodC = allComps(tempCompGood);
        dataIn = mean(EEG.data,3); % Calculate on average data - faster
        S = (EEG.icaweights(goodC,:)*EEG.icasphere)*dataIn(EEG.icachansind,:);
        Winv = EEG.icawinv(:, goodC);
        EEGdataOut = Winv*S;

        subplot(2,2,1)
        plot(EEG.times,EEGdataIn,'b'); grid on; hold on;
        plot([0 0], get(gca,'ylim'),'r--');
        set(gca,'Xlim', options.plotTimeX);
        xlabel('Time (ms)');
        ylabel('Amplitude (\muV)');
        title('Before correction','fontweight','bold');
        yScale = get(gca,'ylim');

        subplot(2,2,2)
        plot(EEG.times,EEGdataOut,'b'); grid on; hold on;
        plot([0 0], yScale,'r--');
        set(gca,'Xlim', options.plotTimeX,'ylim',yScale);
        xlabel('Time (ms)');
        ylabel('Amplitude (\muV)');
        title('After correction','fontweight','bold');
        
        %Plot Yes button
        yesButton = uicontrol('Style', 'pushbutton',...
        'String', {'Yes'},...
        'Position', yesButPos,...
        'fontSize',popFont,...
        'Callback', {@yesCommand});
    
        %Plot No button
        noButton = uicontrol('Style', 'pushbutton',...
        'String', {'No'},...
        'Position', noButPos,...
        'fontSize',popFont,...
        'Callback', {@noCommand});
    
        %Plot Cancel button
        cancelButton = uicontrol('Style', 'pushbutton',...
        'String', {'Cancel'},...
        'Position', cancelButPos,...
        'fontSize',popFont,...
        'Callback', {@cancelCommand});
    
        %Plot component number
        hP1 = uicontrol('style', 'text',...
            'string', ['Are you happy with component removal?'],...
            'position', butTitPos,...
            'BackgroundColor',f1.Color,...
            'fontWeight','bold',...
            'fontSize',compFont);
        
        %Plot components removed
        hP2 = uicontrol('style', 'text',...
            'string', ['Components removed: ',num2str(sum(tempCompBad)),' of ',num2str(length(allComps))],...
            'position', compRemTitPos,...
            'BackgroundColor',f1.Color,...
            'fontSize',compFont);
        
        %Plot variance removed
        hP3 = uicontrol('style', 'text',...
            'string', ['Variance removed: ',num2str(round(sum(varsPerc(tempCompBad)),2)),' %'],...
            'position', varRemTitPos,...
            'BackgroundColor',f1.Color,...
            'fontSize',compFont);
        
        waitfor(gcf, 'Name');
        
        if strcmp(exitScript,'yes')
            error('Plot and remove cancelled. No components removed');
        else
            nComp = 1;
        end
        
    end
         
end

% Make copy of original ICA weights and time course
copyicawinv = EEG.icawinv;
copyicaact = EEG.icaact;

%Remove bad components
allComps = 1:length(compClass);
tempCompBad = compClass > 1;
removeComps = allComps(tempCompBad);
if ~isempty(removeComps)
    EEG = pop_subcomp( EEG,removeComps, 0);
end

% Recalculate EEG.icaact
EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
EEG.icaact = reshape( EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);

% Update EEG.icaCompClass
if isempty(options.compClassInput)
    varsName = fieldnames(EEG);
    varsNameLog = strcmp('icaCompClass',varsName);
    if sum(varsNameLog) == 0
        EEG.icaCompClass = [];
        nameIn = 'Manual1';
    else
        if isempty(EEG.icaCompClass)
            EEG.icaCompClass = [];
            nameIn = 'Manual1';
        else
            varsNameClass = fieldnames(EEG.icaCompClass);
            varsNameClassLog = strncmp('Manual',varsNameClass,6);        
            if sum(varsNameClassLog) == 0
                nameIn = 'Manual1';
            else 
                nameIn = ['Manual', num2str(sum(varsNameClassLog)+1)];
            end
        end
    end
else
    nameIn = options.compClassInput.name;
end

fprintf('Updating EEG.icaCompClass.%s\n',nameIn);
EEG.icaCompClass.(nameIn) = [];
EEG.icaCompClass.(nameIn).name = nameIn;
EEG.icaCompClass.(nameIn).compClass = compClass;
EEG.icaCompClass.(nameIn).compVars = varsPerc;

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

for compNum = 1:length(EEG.icaCompClass.(nameIn).compClass)
    
    %Save components for removal
    if EEG.icaCompClass.(nameIn).compClass(compNum) == 2
        EEG.icaCompClass.(nameIn).reject(1,size(EEG.icaCompClass.(nameIn).reject,2)+1) = compNum;
        EEG.icaCompClass.(nameIn).rejectVars(1,size(EEG.icaCompClass.(nameIn).rejectVars,2)+1) = varsPerc(1,compNum);
    elseif EEG.icaCompClass.(nameIn).compClass(compNum) == 3
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

% Save ICA weights and time course
if strcmp(options.saveWeights,'on')
    EEG.icaCompClass.(nameIn).icawinv = copyicawinv;
    EEG.icaCompClass.(nameIn).icaact = copyicaact;
end

clear('h_gcbf','nextPlot','exitScript');

fprintf('Component removal complete\n');

function nextCommand(src,event)
% persistent check % Shared with all calls of pushbutton1_Callback.
% delay = 0.5; % Delay in seconds.
if isempty(check)
    
    % Hack to ensure waitfor does not miss change in name
    inputName = {'Press next when selection is made. ',...
        'Press next when selection is made.  ',...
        'Press next when selection is made.   ',...
        'Press next when selection is made.    ',...
        'Press next when selection is made.     '};
    
    av = [1 2 3 4 5];
    a_rand = av(randperm(length(av)));
    
    h_gcbf = gcbf;
    h_gcbf.Name = inputName{a_rand(1)};

else
    return
end
end

function back = backCommand(src,event)
% persistent check % Shared with all calls of pushbutton1_Callback.
% delay = 0.5; % Delay in seconds.
if isempty(check)
    
    % Hack to ensure waitfor does not miss change in name
    inputName = {'Press next when selection is made. ',...
        'Press next when selection is made.  ',...
        'Press next when selection is made.   ',...
        'Press next when selection is made.    ',...
        'Press next when selection is made.     '};
    
    av = [1 2 3 4 5];
    a_rand = av(randperm(length(av)));

    h_gcbf = gcbf;
    h_gcbf.Name = inputName{a_rand(1)};
    backPlot = 1;
else
    return
end
end

end


function updateFig(src,event,U,P)

timeIn = strsplit(U.timeEdit.String{1});
freqIn = strsplit(U.freqEdit.String{1});
ampIn = strsplit(U.ampEdit.String{1});
auIn = strsplit(U.auEdit.String{1});

% Update time axes
P.p1.XLim = [str2num(timeIn{1}),str2num(timeIn{2})];
P.p2.XLim = [str2num(timeIn{1}),str2num(timeIn{2})];
P.p3.XLim = [str2num(timeIn{1}),str2num(timeIn{2})];
P.p5.XLim = [str2num(timeIn{1}),str2num(timeIn{2})];

% Update frequency axes
P.p6.XLim = [str2num(freqIn{1}),str2num(freqIn{2})];

% Update amplitude (uV) axes
P.p1.YLim = [str2num(ampIn{1}),str2num(ampIn{2})];
P.p2.YLim = [str2num(ampIn{1}),str2num(ampIn{2})];

% Update amplitude (a.u.) axes
P.p3.YLim = [str2num(auIn{1}),str2num(auIn{2})];

end



function yesCommand(src,event)
close; 
global nextPlot
nextPlot = false;
h_gcbf.Name = 'Press next when selection is made. ';
end

function noCommand(~,event)
close; 
h_gcbf.Name = 'Press next when selection is made. ';
end

function cancelCommand(src,event)
close;
global exitScript
exitScript = 'yes';
h_gcbf.Name = 'Press next when selection is made. ';
end

