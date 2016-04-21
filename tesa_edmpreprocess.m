% tesa_edmpreprocess( )    - This function is needed to run the function 
%                            tesa_edm and to compute the 
%                            preprocessing steps before running  
%                            ICA as advocated in the following papers
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
%   >>  [EEGwhite, inMat, rankMat] = tesa_edmpreprocess( EEG , Nic);
%
% Inputs:
%   EEG                 - EEGLAB EEG structure
%   Nic                 - [Integer] Number of independent components to look for
% 
% Outputs:
%   EEGwhite                 - Whitened data matrix to be used by FastICA
%   inMat                    - Centered matrix   
%   Mean                     - Mean over channels to be used in tesa_edm   
%
% Copyright (C) 2016  Julio Cesar Hernandez Pavon, Aalto University,
% Finland, julio.hpavon@gmail.com
%

function [EEGwhite,inMat,Meanb] = tesa_edmpreprocess(EEG)
%Reshapes 3D matrix to 2D
EEG1=reshape(EEG.data,size(EEG.data,1),[],1);

% centering
Meanb=mean(EEG1,2)*ones(1,size(EEG1,2));
EEG1=EEG1-Meanb; %Centering
T=size(EEG1,2);

% Data to be Whitened
inMat=EEG1;

%Runs singular value decomposition
[U,S,V]=svd(inMat*inMat');
d=diag(sqrt(S));
m=rank(inMat);
% Whitening matrix
EEGwhite=sqrt(T)*diag(1./d(1:m))*U(:,1:m)'*inMat; %Whitening matrix
        
end
