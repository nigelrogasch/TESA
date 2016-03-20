% tesa_pcasuppress()    - suppresses data by top PCA n-dimensions. Note that
%                      PCA is generated from average data, but subtracted
%                      from single trial data.
%                          
%                      Hernandez-Pavon et al (2012) Uncovering neural  
%                      independent components from highly artifactual 
%                      TMS-evoked EEG data. J Neurosci Meth,
%                      209:144-57
% 
% Usage:
%   >>  EEG = tesa_pcasuppress( EEG, timeWin );
%
% Inputs:
%   EEG                 - EEGLAB EEG structure
%   timeWin             - [Vector] - required. Time window for data suppression in ms.
%                       [start, end]
% Outputs:
%   EEG                 - EEGLAB EEG structure
%
% Example:  
%   EEG = tesa_pcasuppress( EEG, [11, 50] );
% 
% See also:
%   tesa_pcacompress, tesa_edm 

% Copyright (C) 2016  Nigel Rogasch & Julio C. Hernandez-Pavon
% Monash University and Aalto University
% nigel.rogasch@monash.edu; julio.hpavon@gmail.com
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

function EEG = tesa_pcasuppress( EEG, timeWin )

%Checks there are enough inputs
if nargin < 2
	error('Not enough input arguments.');
end

%Check that two time points have been specified
if size(timeWin,2) ~= 2
	error('Please provide two time values for the suppression window: [start, end]');
end

%Check that time values given for suppression window are in range of data
if timeWin(1,1) < EEG.times(1,1) || timeWin(1,2) > EEG.times(1,end)
    error('Time values for the suppression window are out of data range. Note that time values are in ms.');
end

%Reshapes data matrix
inMat=reshape(EEG.data,size(EEG.data,1),[],1);

[mint, tInd1]=min(abs(EEG.times-timeWin(1,1))); %Computing the initial time point
[mint, tInd2]=min(abs(EEG.times-timeWin(1,2))); %Computing the terminal time point

%Check that the suppression window does not over lap with removed data
if sum(mean(EEG.data(:,tInd1,:),3)) == 0 || sum(mean(EEG.data(:,tInd2,:),3)) == 0
    error('The window for suppression contains 0s. Please ensure that this window does not overlap with the window of removed data.');
end

% Submatrix according to the time interval where we want to carry out the suppression
 subMat = mean(EEG.data,3); %Mean of the matrix to compute the Submatrix
 Xsub = subMat(:,tInd1:tInd2);%Submatrix
 
%Suppression by PC's
C=0.9999; % Constant 0.9<=C<=1 to preserve the rank of matrix
[U,D,V]=svd(Xsub*Xsub'); % SVD of submatrix to carry out suppression
for N=1:5; %Define the number of PC´s to be removed

    P=eye(size(inMat,1))-C*U(:,1:N)*U(:,1:N)'; % Computing the suppression Matrix
    outMat=P*inMat; % Performing the suppression onto the original matrix

    outData(N).data = reshape(outMat,size(EEG.data,1),size(EEG.data,2),size(EEG.data,3)); % Reshape matrix in to eeglab format
end

titleNames = {'Raw data';'Remove 1 PC';'Remove 2 PCs';'Remove 3 PCs';'Remove 4 PCs';'Remove 5 PCs'};

%Plot the results of removing different componets
figure('Name','Impact of removing PCs on signal','NumberTitle','off'); 
subplot(2,3,1);
plot(EEG.times,mean(EEG.data,3),'b');grid
set(gca,'Xlim',[-100 500]);ylabel('Amplitude (µV)'); xlabel('Time (ms)');
title(titleNames{1,1},'FontWeight','bold');

for a = 1:5;
    subplot(2,3,a+1);
    plot(EEG.times,mean(outData(a).data,3),'b');grid
    set(gca,'Xlim',[-100 500]);ylabel('Amplitude (µV)'); xlabel('Time (ms)');
    title(titleNames{a+1,1},'FontWeight','bold');
end

%Manually select how many components to remove
geometry = {[1 0.5]};
    
uilist = {{'style', 'text', 'string', 'Number of components to remove [0-5]'} ...
          {'style', 'edit', 'string', ''}};

result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Remove PCs');

close;

numPC = str2num(result{1,1});

%Store number of removed PCs
EEG.pcasuppression = numPC;

%Replace data
if isempty(numPC)
    numPC = 0;
elseif numPC > 0 && numPC <= 5
    EEG.data = outData(numPC).data;
elseif numPC > 5
    error('Please enter a number between 0-5');
end

fprintf('Data suppressed by removing top %d dimensions\n', numPC);
        
end
