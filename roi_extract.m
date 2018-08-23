function ROI = roi_extract(S);

% ROIfiles  = cell array of cell arrays, ie set of files for each ROI mask
% Datafiles = cell array of cell arrays, ie set of files for each subject
% TBA***
% uROIvals = cell array of vector of values that define ROI for each ROI mask

% Extract ROI data from image volume(s), either 4D (e.g., resting state
% timecourses) or 3D (e.g., grey matter signal). Gets data from image files
% themselves; does not require SPM.mat. Will extract data from the image files
% based on the dimensionality of the ROI data. 
%
% Outputs struct ROI with various fields, including 'mean' (ROI mean
% timecourses over voxels) and 'svd' (First temporal singular value; for 4D 
% data only).
% (S.filter (default = 0) = if set to 1, this will filter the voxelwise data prior to averaging or PCA
% (S.HPC (default = 100) = highpass cut-off (in seconds)
% (S.LPC (default = 10) = lowpass cut-off (in seconds)
% S.mask =  binary mask image

% by Rik Henson Feb 2012

%% Parameters:


try ROIfiles = S.ROIfiles;
catch
    error('Must pass cell array of ROI mask files')
end

try Datafiles = S.Datafiles;
catch
    error('Must pass cell array of Data files')
end

try output_raw = S.output_raw;
catch
    output_raw = 0;
end

try no_svd = S.no_svd;
catch
    no_svd = 0;
end

try zero_rel_tol = S.zero_rel_tol;
catch
    zero_rel_tol = 0.5;
end

try zero_abs_tol = S.zero_abs_tol;
catch
    zero_abs_tol = 20;
end

try uROIvals = S.uROIvals;
catch
    uROIvals = {};
end

try filter = S.filter;
catch
    filter = 0;
end

if filter==1
try HPC = S.HPC; catch 
    HPC = 1/0.01;
    warning('Assuming highpass cut-off of %d',HPC); 
end

try LPC = S.LPC; catch 
    LPC = 1/0.2;
    warning('Assuming lowpass cut-off of %d',LPC); 
end

try TR = S.TR; catch 
   error('Must pass TR for filtering')
end

end

try mask = S.mask; mask=spm_vol(mask); usemask=1; catch 
  usemask=0; warning('No whole brain mask image used'); 
end


if isempty(uROIvals);
    for m=1:length(ROIfiles)
        ROIvals{m} = [];
    end
else
    if ~iscell(uROIvals)
        for m=1:length(ROIfiles)
            ROIvals{m} = uROIvals;
        end
    else
        for m=1:length(ROIfiles)
            try
                ROIvals{m} = uROIvals{m};
            catch
                ROIvals{m} = [];
            end
        end
    end
end

do_overwrite = 0;
verbose      = 0;
svd_tol  = 1e-6;



ROI = struct(); 

for s=1:length(Datafiles)
    fprintf('Doing Subject %d/%d \n',s,length(Datafiles));
    
    fl = Datafiles{s}';
    VY = spm_vol(fl);
    VY = [VY{:}];
    
    % Get inverse transform (assumes all data files are aligned)
    Yinv  = inv(VY(1).mat);
        
    nr=0;     
    for m=1:length(ROIfiles)
        fprintf('\tDoing Mask %d/%d \n',m,length(ROIfiles));

        VM = spm_vol(ROIfiles{m});
        [YM,mXYZmm] = spm_read_vols(VM);
        ROIfstem = spm_str_manip(ROIfiles{m},'rt');
        
        if isempty(ROIvals{m})
            ROIvals{m} = setdiff(unique(YM(:)),0); 
        elseif strcmp(ROIvals{m},'>0');  % could add more special options here...
            f = find(YM>0);
            YM(f) = 1;
            ROIvals{m} = [1];
        end
    
        % Transform ROI XYZ in mm to voxel indices in data:
        yXYZind = Yinv(1:3,1:3)*mXYZmm + repmat(Yinv(1:3,4),1,size(mXYZmm,2));
        
        fprintf('\t\tDoing ROI (/%d):',length(ROIvals{m}));  
        for r=1:length(ROIvals{m})
            fprintf('.%d',r);
            nr = nr+1;
                        
            ROI(nr,s).ROIfile  = ROIfstem;
            ROI(nr,s).svd_tol  = svd_tol;           
            ROI(nr,s).zero_rel_tol = zero_rel_tol;           
            
            f = find(YM==ROIvals{m}(r));
                      
            d = spm_get_data(VY,yXYZind(:,f)); % Think this is correct!
            
            if usemask==1
            % use only voxels that are part of the mask image.
            maskdata = spm_get_data(mask,yXYZind(:,f)); 
            maskdata=repmat(maskdata,[length(VY)  1]);
            d=d.*maskdata;
            end
            
            Nvox = size(d,2);
            if verbose
                fprintf('Region %d (%s = %d): %d ',r,ROIfstem,ROIvals{m}(r),length(f));
            end
            ROI(nr,s).ROIval    = ROIvals{m}(r);
            ROI(nr,s).XYZ       = mXYZmm(:,f);
            ROI(nr,s).XYZcentre = mean(mXYZmm(:,f),2);
            ROI(nr,s).numvox    = Nvox;
            
            if filter==1
                Ns=size(d,1);
                K   = spm_dctmtx(Ns,Ns);
                nHP = fix(2*(Ns*TR)/HPC + 1);
                nLP = fix(2*(Ns*TR)/LPC + 1);
                K   = K(:,[2:nHP nLP:Ns]);      % Remove initial constant
                R  = eye(Ns) - K*pinv(K);
                d = R*d;

                ROI(nr,s).LPC=LPC;
                ROI(nr,s).HPC=HPC;
                ROI(nr,s).TR=TR;
            end
            
          
            if output_raw
                ROI(nr,s).rawdata = d;
            end
            
            % Check for zero-variance voxels:
            % * this only makes sense for 4D data! 
            % * remove zeros or NaNs for 3D?
            zero_vox = (var(d)==0);
            zero_count = sum(zero_vox); 
            ROI(nr,s).nonzerovox = Nvox-zero_count;

            
            if (zero_count/Nvox > zero_rel_tol) | ((Nvox - zero_count) < zero_abs_tol)
                %if verbose
                    fprintf('(%d nonzero) voxels -- FAILED (%d percent)!\n',Nvox-zero_count,100*zero_rel_tol);
                %end
                ROI(nr,s).mean     = repmat(NaN,size(d,1),1);
                ROI(nr,s).median   = repmat(NaN,size(d,1),1);
                ROI(nr,s).svd      = repmat(NaN,size(d,1),1);
                ROI(nr,s).svd_vox  = repmat(NaN,size(d,2),1);
                ROI(nr,s).svd_pvar = NaN;
            else
                % Remove zero-variance voxels:
                f = setdiff(f,zero_vox);
                d = spm_get_data(VY,yXYZind(:,f));
                
                if verbose
                    fprintf('(%d nonzero) voxels\n',Nvox-zero_count);
                end
                
                % Remove voxel-wise mean if data from several volumes(?)
%                 if numel(d)>max(size(d))
%                     dd = d-repmat(mean(d,1),size(d,1),1);
%                 end
                
                            
                % MEAN/MEDIAN:
                ROI(nr,s).mean   = mean(d,2);
                ROI(nr,s).median = median(d,2);
                
                % SVD (only for timecourses) (note that singular vectors are L2 normalised by spm_en within spm_svd):
                if ~any(size(d)==1) & ~no_svd
                    [U,S,V] = spm_svd(d,svd_tol);
                    if isempty(S)
                        %if verbose
                            fprintf('..SVD FAILED!\n');
                        %end
                        ROI(nr,s).svd      = repmat(NaN,size(d,1),1);
                        ROI(nr,s).svd_vox  = repmat(NaN,size(d,2),1);
                        ROI(nr,s).svd_pvar = NaN;
                        ROI(nr,s).svd_tol  = svd_tol;
                    else
                        ROI(nr,s).svd      = full(U(:,1));
                        ROI(nr,s).svd_vox  = full(V(:,1));
                        ROI(nr,s).svd_pvar = S(1)/sum(full(S(:)));
                        ROI(nr,s).svd_tol  = svd_tol;
                    end
                    %if isempty(S)
                    %    svd_tol = 1e-9; % could increase tolerance?
                    %    [U,S,V] = spm_svd(dd,svd_tol);
                    %end
                else
                    ROI(nr,s).svd      = NaN;
                    ROI(nr,s).svd_vox  = NaN;
                    ROI(nr,s).svd_pvar = NaN;
                    ROI(nr,s).svd_tol  = NaN;
                end
            end
        end
        fprintf('\n')
    end
end

%         % -- aa: should do output by stream --
%         save(sprintf('ROI_%s.mat',ROIfstem),'ROI');
%         
        %==========================================================
        % PLOT

%         % -- aa: leave this out (?) --
%         
%         % Get data:
%         data = [ROI.mean];
%         failed = isnan(sum(data,1));
%         
%         % Plot ROI centres:
%         xyz = [ROI.XYZcentre];
%         figure(1); clf;
%         subplot(3,1,1);
%         for i=1:size(xyz,2)
%             if any(isnan(ROI(i).mean))
%                 plot3(xyz(1,i),xyz(2,i),xyz(3,i),'rx'); hold on;
%             elseif any(isnan(ROI(i).svd_pvar))
%                 plot3(xyz(1,i),xyz(2,i),xyz(3,i),'go'); hold on;
%             else
%                 plot3(xyz(1,i),xyz(2,i),xyz(3,i),'b.'); hold on;
%             end
%         end
%         axis equal
%         view([90 0])
%         title(sprintf('%s (%d ROIs)',cbuid,500-sum(failed)));
%         axis off
%         
%         % Plot MEANS:
%         data = [ROI.mean];
%         failed = isnan(sum(data,1));
%         data = data(:,~failed);
%         subplot(3,2,3);
%         imagesc(data);
%         axis square
%         title('Means'); ylabel('Time'); xlabel('ROI');
%         
%         % Get covariance:
%         dcov = cov(data);
%         subplot(3,2,5);
%         imagesc(dcov);
%         axis square
%         title('Cov(Means)'); ylabel('ROI'); xlabel('ROI');
%         
%         print(1,'-dpng','Figure_ROI_summary.png');
%         
%         % Plot SVD:
%         data = [ROI.svd];
%         failed = isnan(sum(data,1));
%         data = data(:,~failed);
%         subplot(3,2,4);
%         imagesc(data);
%         axis square
%         title('SVDs'); ylabel('Time'); xlabel('ROI');
%         
%         % Get covariance:
%         dcov = cov(data);
%         subplot(3,2,6);
%         imagesc(dcov);
%         axis square
%         title('Cov(SVDs)'); ylabel('ROI'); xlabel('ROI');
%         
%         print(1,'-dpng','Figure_ROI_summary.png');
%         
%         %==========================================================

%fprintf('\n-+- DONE -+-\n');
