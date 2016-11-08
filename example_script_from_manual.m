clear; close all; clc;

filepath = 'H:\\TESA_example_data\\test_TESA\\';
fileprefix = 'example_data';

%Open EEGLAB
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%Step 1: Import data to EEGLAB
% EEG = pop_loadset('filename','example_data.set','filepath','H:\\TESA_example_data\\test_TESA\\');
filename = [fileprefix,'.set'];
EEG = pop_loadset('filename',filename,'filepath',filepath);

%Step 2: Load channel locations (note that you might need to change the file location here for your computer).
EEG = pop_chanedit(EEG, 'lookup','C:\\Program Files\\MATLAB\\eeglab13_5_4b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp');

%Step 3: Remove unused electrodes
EEG = pop_select( EEG,'nochannel',{'31' '32'});

%Save the original EEG locations for use in interpolation later
EEG.allchan = EEG.chanlocs;

%Step 4: Find TMS pulses (if no triggers recorded) - skipped
%EEG = pop_tesa_findpulse( EEG, 'Cz', 'refract', 4, 'rate', 10000, 'tmsLabel', 'TMS', 'plots', 'on');

%Step 5: Remove bad electrodes
EEG = pop_rejchan(EEG, 'elec',[1:size(EEG.data,1)] ,'threshold',5,'norm','on','measure','kurt');

%Step 6: Epoch data (-1000 to 1000 ms)
EEG = pop_epoch( EEG, {  'R128'  }, [-1  1], 'epochinfo', 'yes');

%Step 7: Demean data (-1000 ms to 1000 ms)
EEG = pop_rmbase( EEG, [-1000  1000]);

%Step 8: Remove TMS pulse artifact and peaks of TMS-evoked muscle activity (-2 to 10 ms)
EEG = pop_tesa_removedata( EEG, [-2 10] );

%Step 9: Interpolate missing data around TMS pulse
EEG = pop_tesa_interpdata( EEG, 'cubic', [1,1] );

%Step 10: Downsample data (5000 Hz to 1000 Hz)
EEG = pop_resample( EEG, 1000);

%Step 11: Remove bad trials  
EEG = pop_jointprob(EEG,1,[1:size(EEG.data,1)] ,5,5,0,0);
pop_rejmenu(EEG,1);
pause_script = input('Highlight bad trials, update marks and then press enter');
EEG.BadTr = unique([find(EEG.reject.rejjp==1) find(EEG.reject.rejmanual==1)]);
EEG = pop_rejepoch( EEG, EEG.BadTr ,0); 

%Step 12: Replace interpolated data around TMS pulse with constant amplitude data (-2 to 10 ms)
EEG = pop_tesa_removedata( EEG, [-2 10] ); 

%Step 13: Remove TMS-evoked muscle activity (using FastICA and auto component selection)
EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
EEG = pop_tesa_compselect( EEG,'comps',15,'figSize','small','plotTimeX',[-200 500],'plotFreqX',[1 100],'tmsMuscle','on','tmsMuscleThresh',8,'tmsMuscleWin',[11 30],'tmsMuscleFeedback','off','blink','off','blinkThresh',2.5,'blinkElecs',{'Fp1','Fp2'},'blinkFeedback','off','move','off','moveThresh',2,'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','off','muscleThresh',0.6,'muscleFreqWin',[30 100],'muscleFeedback','off','elecNoise','off','elecNoiseThresh',4,'elecNoiseFeedback','off' );

%Step 14: Extend data removal to 15 ms (-2 to 15 ms)
EEG = pop_tesa_removedata( EEG, [-2 15] );

%Step 15: Interpolate missing data around TMS pulse
EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );

%Step 16: Bandpass (1-100 Hz) and bandstop (48-52 Hz) filter data
EEG = pop_tesa_filtbutter( EEG, 1, 100, 4, 'bandpass' ); 
EEG = pop_tesa_filtbutter( EEG, 48, 52, 4, 'bandstop' );

%Step 17: Replace interpolated data around TMS pulse with constant amplitude data (-2 to 15 ms)
EEG = pop_tesa_removedata( EEG, [-2 15] );

%Step 18: Remove all other artifacts (using FastICA and auto component selection)
EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
EEG = pop_tesa_compselect( EEG,'compCheck','on','comps',[],'figSize','small','plotTimeX',[-200 500],'plotFreqX',[1 100],'tmsMuscle','on','tmsMuscleThresh',8,'tmsMuscleWin',[11 30],'tmsMuscleFeedback','off','blink','on','blinkThresh',2.5,'blinkElecs',{'Fp1','Fp2'},'blinkFeedback','off','move','on','moveThresh',2,'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','on','muscleThresh',0.6,'muscleFreqWin',[30 100],'muscleFeedback','off','elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off' );

%Step 19: Interpolate missing data around TMS pulse
EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );

%Step 20: Interpolate missing channels
EEG = pop_interp(EEG, EEG.allchan, 'spherical');

%Step 21: Re-reference to average
EEG = pop_reref( EEG, []);

%Save point
filename = [fileprefix,'_cleaned.set'];
EEG = pop_saveset( EEG, 'filename',filename,'filepath',filepath);

%Step 22: Extract a region of interest
EEG = pop_tesa_tepextract( EEG, 'ROI', 'elecs', {'P3','CP1','P1','CP3'}, 'tepName', 'Parietal' );

%Step 23: Find peaks of interest
EEG = pop_tesa_peakanalysis( EEG, 'ROI', 'positive', [40 80 200], [30 50;70 90;180 220], 'method' ,'largest', 'samples', 5, 'tepName', 'Parietal' );
EEG = pop_tesa_peakanalysis( EEG, 'ROI', 'negative', [20 60 100], [10 30;50 70;90 110], 'method' ,'largest', 'samples', 5, 'tepName', 'Parietal' );

%Step 24: Output peak analysis in table
output = pop_tesa_peakoutput( EEG, 'tepName', 'Parietal', 'calcType', 'amplitude', 'winType', 'individual', 'averageWin', [], 'fixedPeak', [], 'tablePlot', 'on' );

%Step 25: Plot the results
pop_tesa_plot( EEG, 'tepType', 'ROI', 'tepName', 'Parietal', 'xlim', [-100 500], 'ylim', [], 'CI','off','plotPeak','on' );

%Save point
filename = [fileprefix,'_cleaned_analysed.set'];
EEG = pop_saveset( EEG, 'filename',filename,'filepath',filepath);
