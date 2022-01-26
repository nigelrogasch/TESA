
clear; close all; clc;

% This pipeline might not fit all purposes or datasets. It is merely a
% practical example showing how to apply the TESA SOUND and TESA SSP-SIR functions
% to TMS-EEG data.

% Tuomas Mutanen, 17-05-2020


%% Step 1: Import data to EEGLAB

filepath = '/your/path/to/tesa/example/data/';
fileprefix = 'example_data';

%Open EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

filename = [fileprefix,'.set'];
EEG = pop_loadset('filename',filename,'filepath',filepath);

%% Step 2: Load channel locations (note that you might need to change the file location here for your computer).

eeglab_path = '/your/path/to/eeglab/';
EEG = pop_chanedit(EEG, 'lookup', fullfile(eeglab_path, 'plugins','dipfit','standard_BESA','standard-10-5-cap385.elp') );

%% Step 3: Remove unused electrodes

EEG = pop_select( EEG,'nochannel',{'31' '32'});

%% Step 4: Epoch data (-1000 to 1000 ms)

EEG = pop_epoch( EEG, {  'R128'  }, [-1  1], 'epochinfo', 'yes');

%% Step 5: Baseline correct the data (-1000 ms to -5 ms)

EEG = pop_rmbase( EEG, [-1000  -2]);

%% Step 6: Remove TMS pulse artifact (-2 to 4 ms)

% NOTE: Here, we remove a much shorter time window than in the original TESA pipeline, to ensure that we
% have enough information about the muscle artifacts that will be removed
% with SSP-SIR later on.

EEG = pop_tesa_removedata( EEG, [-2 4] );

%% Step 7: Interpolate the missing data around the TMS pulse

EEG = pop_tesa_interpdata( EEG, 'cubic', [1,1] );

%% Step 8: Downsample data (5000 Hz to 1000 Hz)

EEG = pop_resample( EEG, 1000);

%% Step 9: Remove bad trials  

EEG = pop_jointprob(EEG,1,[1:size(EEG.data,1)] ,5,5,0,0);
pop_rejmenu(EEG,1);
pause_script = input('Highlight bad trials, update marks and then press enter');
EEG.BadTr = unique([find(EEG.reject.rejjp==1) find(EEG.reject.rejmanual==1)]);
EEG = pop_rejepoch( EEG, EEG.BadTr ,0); 

%% Step 10: Remove blinks and ocular artifacts prior to SOUND with ICA (using FastICA and auto component selection)

EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );

EEG = pop_tesa_compselect( EEG,'compCheck','on','comps',[],'figSize','large',...
    'plotTimeX',[-200 500],'plotFreqX',[1 100],'tmsMuscle','off','tmsMuscleThresh',...
    8,'tmsMuscleWin',[11 30],'tmsMuscleFeedback','off','blink','on','blinkThresh',...
    2.5,'blinkElecs',{'Fp1','Fp2'},'blinkFeedback','off','move','on','moveThresh',...
    2,'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','off','elecNoise','off',...
    'elecNoiseThresh',4,'elecNoiseFeedback','off' );

%% Step 11: Perform the baseline correction again after ICA:

EEG = pop_rmbase( EEG, [-1000 -5] ,[]);

% Visualize the data:

figure; pop_timtopo(EEG, [-50  200], [10  40  60  80], 'ERP data and scalp maps of  resampled');

%% Step 12: Detect and clean the noise components and noisy channels with SOUND (default settings) 

EEG = pop_tesa_sound(EEG, 'lambdaValue', 0.1, 'iter', 5, 'leadfieldInFile', [], 'leadfieldChansFile', [], 'replaceChans', [], 'multipleConds', []);

% Visualize the data:

figure; pop_timtopo(EEG, [-50  200], [10  40  60  80], 'ERP data and scalp maps of  resampled');


%% Step 13: Remove residual muscle artifact with the SSP-SIR (default settings):

EEG = pop_tesa_SSPSIR( EEG, 'artScale','automatic','PC',[]);

% Visualize the data:
figure; pop_timtopo(EEG, [-50  200], [10  40  60  80], 'ERP data and scalp maps of  resampled');


%%  Step 14: Bandpass (1-100 Hz) and bandstop (48-52 Hz) filter data

EEG = pop_tesa_filtbutter( EEG, 1, 100, 4, 'bandpass' ); 
EEG = pop_tesa_filtbutter( EEG, 48, 52, 4, 'bandstop' );

% Visualize the data:
figure; pop_timtopo(EEG, [-50  200], [10  40  60  80], 'ERP data and scalp maps of  resampled');




