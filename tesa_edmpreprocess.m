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
%  compVal             - [Integer] Number of dimensions to compress data
%                           Default = 25.
% 
% Outputs:
%   EEGwhite                 - Whitened data matrix to be used by FastICA
%   outMat                   - Compressed data matrix     
%   compVal                  - number of compressed dimensions   
%
% Copyright (C) 2016  Julio Cesar Hernandez Pavon, Aalto University,
% Finland, julio.hpavon@gmail.com
%

function [EEGwhite,outMat,compVal] = tesa_edmpreprocess( EEG, compVal );

%Sets truncate to default value if not specified
if nargin <2;
	compVal = 25;
end
if isempty(compVal)
    compVal = 25;
end

%Reshapes 3D matrix to 2D
EEG1=reshape(EEG.data,size(EEG.data,1),[],1);

%Average reference and centering
M=size(EEG1,1);
EEG1=EEG1-ones(M,1)*mean(EEG1,1);%Reference potential level
EEG1=EEG1-mean(EEG1,2)*ones(1,size(EEG1,2)); %Centering
T=size(EEG1,2);

% Data to be compressed
inMat=EEG1;

%Checks that number of dimensions is larger than compression value
rankMat = rank(inMat);
if rankMat <= compVal
    error('Dimension of data (%d) is lower than the compression value (%d). Function terminated.', rankMat, compVal);
end

%Runs singular value decomposition
[U,S,V]=svd(inMat*inMat');
d=diag(sqrt(S));

%Compresses data by truncating to 'compVal' dimensions 
C=U(:,1:compVal)*U(:,1:compVal)';
outMat=C*inMat;

% Whitening matrix
EEGwhite=sqrt(T)*diag(1./d(1:compVal))*U(:,1:compVal)'*outMat; %Whitening matrix
fprintf('Data compressed to %d dimensions\n', compVal);
        
end
