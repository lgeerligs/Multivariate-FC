%% This script is a demonstration of how to analyze functional connectivity data using distance correlation 
% Linda Geerligs (lindageerligs@gmail.com) - 24-03-2016 
% Many thanks to Rik Henson for providing parts of this script

% script uses functions from SPM12, so this needs to be on your path

%% predefined variables - can be downloaded from http://imaging.mrc-cbu.cam.ac.uk/imaging/Geerligs_DistCor
clear all; 

% TR
TR=1.97;

% cut-off for high-pass filter
HP=0.008; 

% ROI definition file
roifile='Power_ROIs.nii';

% get the data and motion - CSF and WM signals
datafile='EPI_data.nii';
load(['nuisance_signals.mat']);
MO=load(['realignment_parameters.txt']);

%% extract ROI data
hdr=spm_vol(datafile);datacell=[];
for i=1:length(hdr);
    datacell{i}=[datafile ',' num2str(i)];
end

%extract data from each ROI
SR.output_raw=1;
SR.Datafiles{1}=datacell;SR.ROIfiles={roifile};
ROI = roi_extract(SR);

%% prepare for connectivity analyses
%remove voxels that have no signal 
for i=1:length(ROI)
    notempty=find(sum(ROI(i).rawdata)~=0);
    data{i}=ROI(i).rawdata(:,notempty);
end

% get nuisance regressors
Ns=length(MO);
dm = [zeros(1,6); diff(MO,1,1)];
dc = [zeros(1,2); diff([CSF WM],1,1)];
X0r1=[ones(Ns,1) CSF WM MO dm dc  CSF.^2 WM.^2 MO.^2 dm.^2 dc.^2];

% Create a DCT filter and add it to the nuisance regressors
Kall   = spm_dctmtx(Ns,Ns);
nHP = fix(2*(Ns*TR)/(1/HP) + 1);
KHP   = Kall(:,2:nHP);
X0rHP=[X0r1 KHP];

% Create comprehensive model of residual autocorrelation (to counter Eklund et al, 2012, Neuroimage)
T     = (0:(Ns - 1))*TR;                    % time
d     = 2.^(floor(log2(TR/4)):log2(64));    % time constants (seconds)
Q    = {};                                  % dictionary of components
for i = 1:length(d)
    for j = 0:1
        Q{end + 1} = toeplitz((T.^j).*exp(-T/d(i)));
    end
end

%% run nuisance regression and connectivity analysis

for i=1:length(data)

    %remove each voxels' mean
    dat=detrend(data{i},'constant');

    %compute covariance across voxels
    CY = cov(dat');
    [V, ~]   = spm_reml_noprint(CY,X0rHP,Q,1,0,4);   % this is just version of spm_reml with fprintf commented out to speed up
    W_HP       = spm_inv(spm_sqrtm(V));
    
    % get the residuals of the regression+autoregressive model 
    [res_HP]   = spm_ancova1(W_HP*X0rHP,speye(Ns,Ns),W_HP*dat);
    
    M(i,:)=mean(res_HP');
    
    %variance normalize the data
    res_HP=zscore(res_HP);
    
    %get the U-centered distance correlation matrix
    a = pdist2(res_HP, res_HP); 
    D1{i} = Ucenter(a, size(dat,1));
    fprintf('%d.',i)
end
fprintf('\n')

n=size(data{1},1);
Dcor=zeros(length(data), length(data));
for i=1:length(data) 
    for j=1:length(data)
        
        %multivariate U-centered distance correlation
        A=D1{i}; B=D1{j};
        dcovXY = sum(sum(A.*B)) ./ (n*(n-3));
        dvarX = sum(sum(A.*A)) ./ (n*(n-3));
        dvarY = sum(sum(B.*B)) ./ (n*(n-3));
        Dcor(i,j)=dcovXY / sqrt(dvarX * dvarY);
        
        if Dcor(i,j)<0
            Dcor(i,j)=0;
        else
            Dcor(i,j)=sqrt(Dcor(i,j));
        end        
    end
    fprintf('%d.',i)
end
fprintf('\n')

Pcor=corr(M');

%% plot the connectivity estimates
figure; subplot(1,2,1); imagesc(Pcor); caxis([-0.1 1]);colormap(jet); set(gca,'FontSize',20);title('Pcor'); colorbar
set(gca,'dataAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);
subplot(1,2,2); imagesc(Dcor); caxis([0 1]);colormap(jet); set(gca,'FontSize',20);title('Dcor'); colorbar
set(gca,'dataAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);
