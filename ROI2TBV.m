% function ROI2TBV(structMRI,functMRI,ROI,options)
%
% This function takes an ROI as drawn in MRIcron (based on a structural
% MRI) and coregisters and transforms this ROI to the space and resolution
% of the EPI (functMRI).
%
% Tobias Hauser & Benjamin Chew, 05/2017
function ROI2TBV(structMRI,functMRI,ROI,options)

%% general settings
if nargin < 3
    error('not enough input arguments')
end
if nargin < 4
    options = [];
    disp('using standard settings')
end
if ~isfield(options,'flipX')
    options.flipX = 0;
end
if ~isfield(options,'flipY')
    options.flipY = 1;
end
if ~isfield(options,'flipZ')
    options.flipZ = 0;
end
if ~isfield(options,'shiftX')
    options.shiftX = -1;    % default value: -1 because TBV voxels start from 0
end
if ~isfield(options,'shiftY')
    options.shiftY = -1;    % default value: -1 because TBV voxels start from 0
end
if ~isfield(options,'shiftZ')
    options.shiftZ = -1;    % default value: -1 because TBV voxels start from 0
end
if ~isfield(options,'shiftZ')
    options.shiftZ = -1;    % default value: -1 because TBV voxels start from 0
end
if ~isfield(options,'ROI_threshold')
    options.ROI_threshold = 0.9;    % cutoff for selecting ROI (proportion of max)
end
if ~isfield(options,'out_fname')
    options.out_fname = 'newRoi.roi';
end
if ~isfield(options,'out_dir')
    [p_roi] = fileparts(ROI);
    options.out_dir = [p_roi '\'];
end
if ~isfield(options,'coreg')
    options.coreg = 1;
end
if ~isfield(options,'verbose')
    options.verbose = 1;
end


% addpath('C:\myDocuments\work\Tools\spm\')
% roi_name = 'C:\Users\Tobias Hauser\Desktop\BEN\2017_05_25\T1_RightVentricle_LeftEye.nii';
% epi_name = 'C:\Users\Tobias Hauser\Desktop\BEN\2017_05_25\fMP01859-0006-00001-000001-01.nii';



%% coregister structural & ROI to EPI
if options.coreg
    disp('running coregistration')
    clear matlabbatch
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[functMRI ',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {[structMRI ',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {[ROI ',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    spm_jobman('initcfg');
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch(1));
end


%% load coregistered ROI and EPI
roi = spm_vol(ROI);
[roi_data,roi_XYZ] = spm_read_vols(roi);
disp('roi loaded')

epi = spm_vol(functMRI);
[epi_data,epi_XYZ] = spm_read_vols(epi);
disp('epi loaded')

%% get relevant ROI voxels
roi_vox = roi_XYZ(:,find(roi_data>options.ROI_threshold*max(roi_data(:))));
if options.verbose
    figure(); plot3(roi_vox(1,:),roi_vox(2,:),roi_vox(3,:),'d'); hold on
end

%% find the EPI-voxel index for each ROI-voxel
epi_idx = [];
for v = 1:length(roi_vox)
    [closest,cval,err] = findcEucl(epi_XYZ,roi_vox(:,v));
    
    if err > 1.5    % excludes voxels that are >1.5mm away from origin
        epi_idx(v) = nan;
    else
        epi_idx(v) = closest;
    end
end
epi_idx = unique(epi_idx);  % remove redundant voxels
epi_idx(find(isnan(epi_idx))) = [];
if options.verbose
    plot3(epi_XYZ(1,epi_idx),epi_XYZ(2,epi_idx),epi_XYZ(3,epi_idx),'d','color','r'); legend('original','voxels in EPI')
end

%% transform in coordinates or native voxels
XYZ_out = [];
for v = 1:length(epi_idx)
    [XYZ_out(v,1) XYZ_out(v,2) XYZ_out(v,3)] = ind2sub(size(epi_data),epi_idx(v));
end
if options.verbose
    figure(); plot3(XYZ_out(:,1),XYZ_out(:,2),XYZ_out(:,3),'d'); legend('EPI voxels in native space')
end
disp('coordinates transformed')

%% flip and adjust coordinates for TBV
if options.flipX
    XYZ_out(:,1) = (size(epi_data,1) + 1) - XYZ_out(:,1);
end
if options.flipY
    XYZ_out(:,2) = (size(epi_data,2) + 1) - XYZ_out(:,2);
end
if options.flipZ
    XYZ_out(:,3) = (size(epi_data,3) + 1) - XYZ_out(:,3);
end

XYZ_out(:,1) = XYZ_out(:,1) + options.shiftX;
XYZ_out(:,2) = XYZ_out(:,2) + options.shiftY;
XYZ_out(:,3) = XYZ_out(:,3) + options.shiftZ;

%% write output file for TVB
% write header
txt = ['FileVersion:           5\n\nSaveVoxelsInROIs:      1\n\nSaveSortedVoxelList:   0\n\nSaveROIMaxPSC:         0\n\nNrOfTCPlots:           1\n\nNrOfTCs:               1\nFromSlice:             0\nLeft:                  0\nRight:                 0\nTop:                   0\nBottom:                0\nNrOfVoxels:            '];
txt = [txt num2str(size(XYZ_out,1)) '\n'];
fid=fopen([options.out_dir options.out_fname],'wt');
fprintf(fid, [txt]);
fclose(fid);

% write voxels
dlmwrite([options.out_dir options.out_fname],XYZ_out,'-append',...
    'delimiter','\t','newline','pc')
disp('file written. all done.')

%% helper functions
function [closest,cval,err] = findcEucl(mat,val)
% [closest,cval,err] = findc(mat,val)
%
% returns the index (closest), value (cval) and error (err) of the array closest to val
% in mat
%
% Tobias Hauser, 23/05/2017


eucl_dist = sqrt(sum((mat - repmat(val,1,size(mat,2))).^2));

[err,closest] = min(eucl_dist);
cval = mat(:,closest);
end
end