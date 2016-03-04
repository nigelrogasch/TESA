% tesa_edmpreprocess( )    - This function is needed to run the function 
%                            tesa_edm and to compute the 
%                            preprocessing steps before running  
%                            FastICA as advocated in the following papers
%
%                           Korhonen, Hernandez-Pavon et al (2011) 
%                           Removal of large muscle artifacts 
%                           from transcranial magnetic stimulation-evoked EEG 
%                           by independent component analysis. Med Biol Eng
%                           Compt, 49:397-407.
%
%                           Hernandez-Pavon et al (2012) Uncovering neural  
%                           independent components from highly artifactual 
%                           TMS-evoked EEG data. J Neurosci Meth,
%                           209:144-57
% 
% Usage:
%   >>  EEG = tesa_edmpreprocess( EEG );
%   >>  EEG = tesa_edmpreprocess( EEG , compVal);
%
% Inputs:
%   EEG                 - EEGLAB EEG structure
% 
% Outputs:
%   EEGwhite                 - Whitened data matrix to be used by FastICA
%   inMat                    - Centered matrix   
%   rankMat                  - Rank of the data matrix 
%
% Copyright (C) 2016  Julio Cesar Hernandez Pavon, Aalto University,
% Finland, julio.hpavon@gmail.com
%

function [EEGwhite,inMat,rankMat] = tesa_edmpreprocess(EEG, Nic);
%Reshapes 3D matrix to 2D
EEG1=reshape(EEG.data,size(EEG.data,1),[],1);

% centering
M=size(EEG1,1);
EEG1=EEG1-mean(EEG1,2)*ones(1,size(EEG1,2)); %Centering
T=size(EEG1,2);

% Data to be Whitened
inMat=EEG1;

%Checks that number of ICs is smaller than the rank of the data matrix
rankMat = rank(inMat);
if Nic > rankMat
    error('Number of ICs (%d) is bigger than than the rank of the data matrix (%d). Function terminated.',Nic,rankMat);
end

%Runs singular value decomposition
[U,S,V]=svd(inMat*inMat');
d=diag(sqrt(S));

% Whitening matrix
EEGwhite=sqrt(T)*diag(1./d)*U'*inMat; %Whitening matrix
        
end
