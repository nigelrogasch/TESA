function tesa_plotcomps( EEG, time, frequency )

    %Calculates time courses
    timeCourse = arrayfun(@(x)mean(reshape(EEG.icaweights(x,:)*EEG.data(:,:),EEG.pnts,EEG.trials),2),1:size(EEG.icawinv,2),'UniformOutput', false);
    
    t1 = time(1,1);
    t2 = time(1,2);
    
    %Calculates FFT of components
    T = 1/EEG.srate;             % Sample time
    L = size(EEG.times,2);            % Length of signal
    NFFT = 2^nextpow2(L);         % Next power of 2 from length of y
    f = EEG.srate/2*linspace(0,1,NFFT/2+1); % Frequencies
    
    f1 = frequency(1,1);
    f2 = frequency(1,2); 
    
    for a=1:size(timeCourse,2);
        y = cell2mat(timeCourse(:,a));
        Y(a,:) = fft(zscore(y),NFFT)/L;
    end;
    
    Yout = (abs(Y).^2); 
    [c index]=min(abs(f-(0.5)));
    
    for x = 1:size(Yout,1);
        freq=f1:0.5:f2;
        for a=1:size(freq,2);
            [c index1]=min(abs(f-((freq(1,a)-0.25))));
            [c index2]=min(abs(f-((freq(1,a)+0.25))));
            Y2(x,a)=mean(Yout(x,index1:index2),2); %creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
        end;
    end;
    
    %Creates matrices for generating sufficient subplots based on 6x6 grid
    
    figNum = 1:1:ceil(size(EEG.icawinv,2)*4/36); %number of figures needed
    
    matNum = reshape(1:36,6,6)';
    
    timeNum = matNum(1:2:end, 1:2:end)'; %positions of time plots
    timeNum = timeNum(:)';
    
    trialNum = matNum(2:2:end, 1:2:end)'; %position for trial plots
    trialNum = trialNum(:)';
    
    fftNum = matNum(2:2:end, 2:2:end)'; %positions of fft plots
    fftNum = fftNum(:)';
    
    topNum = matNum(1:2:end, 2:2:end)'; %positions of topoplots
    topNum = topNum(:)';
    
    compNum = 1:1:size(EEG.icawinv,2); 
    compNum(size(compNum,2)+1:max(figNum)*size(timeNum,2)) = NaN;
    compNum = reshape(compNum,size(timeNum,2),max(figNum))'; %component positions for each figure
    
    [x, tp1] = find(EEG.times == t1);
    [x, tp2] = find(EEG.times == t2);

    %plot the figure
    
    close all;
    
    for a = 1:size(figNum,2)
        figure;
        for b = 1:size(timeNum,2)
            if ~isnan(compNum(a,b))
                temp = cell2mat(timeCourse(1,compNum(a,b)));
                subplot(6,6,timeNum(1,b));
                plot(EEG.times,temp);
                set(gca,'Xlim', [t1 t2], 'Xtick', [t1:t2-t1:t2]);
                title(num2str(compNum(a,b)),'fontweight', 'bold');
                
                temp1 = reshape(EEG.icaweights(compNum(a,b),:)*EEG.data(:,:),EEG.pnts,EEG.trials);
                temp1 = temp1(tp1:tp2,:);
                subplot(6,6,trialNum(1,b));
                imagesc(temp1','XData', [t1 t2]);
                caxis([-max(abs(temp1(:))), max(abs(temp1(:)))]);
                
                subplot(6,6,fftNum(1,b));
                plot(freq,Y2(compNum(a,b),:));
                set(gca,'Xlim', [0 f2], 'Xtick', [0:f2:f2], 'yticklabel', []);
                
                subplot(6,6,topNum(1,b));
                topoplot(EEG.icawinv(:,compNum(a,b)),EEG.chanlocs);
            end
        end
    end
    
    
end