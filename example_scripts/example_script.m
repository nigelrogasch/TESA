clear all; close all; clc;

% ##### SETTINGS #####

%Enter path to file
path = 'C:\Users\Nigel\Desktop\TESA_example_data\';

%%
% ##### LOAD FILE AND CHANNEL LOCATIONS #####

% %Load file
% EEG = pop_loadbv(path, 'testData.vhdr');
% 
% %Load in channel data
% EEG = pop_chanedit(EEG,  'lookup', 'standard-10-5-cap385.elp');
% 
% %Remove eye electrodes
% EEG = pop_select(EEG,'nochannel',{'31','32'});
% 
% %Save channels - will need this when interpolating missing channels
% EEG.allchan = EEG.chanlocs;

%%
% ##### FIND TMS PULSES (TESA) #####
% !!!!! Optional - only if required

% %-----
% %Option 1: Find artefacts using first derivative and voltage threshold
%     EEG = pop_tesa_findpulse( EEG, 'Cz');
% 
% %-----
% %Option 2: Find artefacts using peak data
% %   EEG = pop_tesa_findpulsepeak( EEG, 'Cz' );

%%
% ##### EPOCH AND BASELINE CORRECT DATA #####

% %Epoch (-1 s to 1 s)
% EEG = pop_epoch( EEG, {  'TMS'  }, [-1  1], 'epochinfo', 'yes');
% 
% %Baseline correct (-500 ms to -10 ms)
% EEG = pop_rmbase( EEG, [-500  -10]);

%%
% ##### CORRECT MARKING OF TMS PULSES, MARK PAIRED PULSES (TESA) #####
% !!!!! Optional - only if required

% %Find trigger and take new epoch between 0.8 s to 0.8 s
% EEG = pop_tesa_fixtrigger( EEG, 'Cz', [-0.8,0.8] );

%%
% ##### DOWNSAMPLE DATA #####

% %Downsample from 10,000 Hz to 1,000 Hz
% EEG = pop_resample( EEG, 1000);

%%
% ##### LOAD SAMPLE DATA #####

%Call eeglab
eeglab

%Load data
EEG = pop_loadset( 'filename', 'testData_ep_bc_ds.set', 'filepath', path);

%%
% ##### REMOVE THE TMS PULSE ARTEFACT (TESA) #####

%Remove TMS pulse artefact and peak of muscle artefacts taking in to
%account ringing artefacts from downsample (-12 ms to 12 ms)
EEG = pop_tesa_removedata(EEG,[-12,12]);

%%
% ##### REMOVE BAD ELECTRODES #####

%Check for bad electrodes
pop_eegplot( EEG, 1, 1, 1);
pause_script = input('Press enter when you are ready to continue.');

%Remove bad electrodes
answer = inputdlg('Enter bad channels', 'Bad channel removal', 1);
EEG.badChan = strsplit(answer{1});
close all;
EEG = pop_select( EEG,'nochannel',EEG.badChan);

%%
% ##### REMOVE BAD TRIALS #####

%Select bad trials
pop_rejmenu( EEG, 1);
pause_script = input('Highlight bad trials, update marks and then press enter');

%Reject bad trials
EEG.BadTr = unique(find(EEG.reject.rejmanual==1));
EEG = pop_rejepoch(EEG,EEG.BadTr,0);

%%
% ##### MINIMISE TMS-EVOKED MUSCLE ARTEFACT (TESA) #####

%-----
%Option 1: Use fastICA to remove muscle component (EEGLAB, TESA)

%     %Run fastica (EEGLAB)
%     EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', 'tanh');
% 
%     %Sort components by time course variance (TESA)
%     EEG = pop_tesa_sortcomps(EEG);
% 
%     %Select TMS-evoked muscle components
%     pop_selectcomps(EEG, [1:size(EEG.icawinv,2)] );
%     pause_script = input('Highlight bad trials, update marks and then press enter');
% 
%     %Removes TMS-evoked muscle components
%     EEG.tmsComp = unique(find(EEG.reject.gcompreject==1));
%     EEG = pop_subcomp( EEG,EEG.tmsComp, 0);

%-----
%Option 2: Use Enhanced Deflation Method (EDM) to remove muscle component (TESA)


%-----
%Option 3: Use PCA to remove muscle component (TESA)
    EEG = pop_tesa_pcacompress( EEG, 30 ); %OPTIONAL - can compress data to top 25-30 dimensions prior to supression
    EEG = pop_tesa_pcasupress(EEG,[13,50]);

%-----
%Option 4: Use PCA with wavelet filtering to remove muscle component (TESA)

%-----
%Option 5: Use linear or exponential modeling to remove muscle component (TESA)
%     EEG = pop_tesa_detrend(EEG,'double',[13,999]);

%-----
%Option 6: Use median filter to remove muscle component (TESA)
%     EEG = pop_tesa_medianfilt(EEG,[13,50]);

%%
% ##### INTERPOLATE MISSING DATA AROUND TMS ARTEFACT (TESA) #####

%Interpolate missing data
EEG = pop_tesa_interpdata(EEG,'linear');

%%
% ##### FILTER THE DATA #####

%Band pass filter data (1-100 Hz)
EEG = pop_eegfiltnew(EEG, 1, 100, 3300, 0, [], 0);

%Band stop filter data (48-52 Hz) 
EEG = pop_eegfiltnew(EEG, 48, 52, 1650, 1, [], 0);

%%
% ##### RUN FASTICA AND SORT COMPONENTS (TESA) #####

%Run fastica (EEGLAB)
EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', 'tanh');

%Sort components by time course variance (TESA)
EEG = pop_tesa_sortcomps(EEG);

%%
% ##### REMOVE BAD COMPONENTS #####

%Select bad components
pop_selectcomps(EEG, [1:size(EEG.icawinv,2)] );
pause_script = input('Highlight bad trials, update marks and then press enter');

%Removes bad components
EEG.badComps = unique(find(EEG.reject.gcompreject==1));
EEG = pop_subcomp( EEG,EEG.badComps, 0);

%%
% ##### INTERPOLATE MISSING ELECTRODES #####

%Interpolate missing electrodes
EEG = pop_interp(EEG, EEG.allchan, 'spherical');

%%
% ##### RE-REFERENCE TO AVERAGE OF ALL CHANNELS #####

%Re-reference to average
EEG = pop_reref( EEG, []);

%%
% ##### ANALYSE TEPS WITH ROI(TESA) #####

EEG = pop_tesa_tepextract(EEG,'ROI','roi',{'P1','P3','P5','CP1','CP3','CP5'});
EEG = pop_tesa_peakanalysis(EEG,'ROI','positive',[40,80],[35,45;75,85]);
EEG = pop_tesa_peakanalysis(EEG,'ROI','negative',[20,60,100],[15,25;55,65;95,105]);

%%
% ##### ANALYSE TEPS WITH GMFA(TESA) #####

EEG = pop_tesa_tepextract(EEG,'GMFA');
EEG = pop_tesa_peakanalysis(EEG,'GMFA','positive',[20,40,60,80,120],[15,25;35,45;55,65;75,85;115,125]);

%%
% ##### OUTPUT ANALYSIS RESULTS TO WORKSPACE (TESA) #####

output = pop_tesa_peakoutput(EEG,'averageWin',5);

%%

% ##### PLOT DATA #####

tesa_plot(EEG,'input','data');
tesa_plot(EEG,'input','ROI','roiName','R1','CI','on');
tesa_plot(EEG,'input','GMFA','plotPeak','on');



