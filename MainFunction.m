%% INSTRUCTIONS
% For Z/R/MTC_map 
% Draw all the slices at one single run

%% STEP1: Load files 

close all;
clc;
clear all;
addpath('toolbox');

%   Load the mask for the whole braiin
load('mask.mat');
Al_allS = mask;

%   Load the T1map, 64*64*25
load('T1map.mat');

%   Load frequency list
load('crlist.mat');

%   Load frequency list
load('CEST_images.mat');
%%  Motion Correction
%   Convert the original CEST image files to .nii
 
v = squeeze(v);

for n_offset = 1: size(v,4)
        M0_volume_Ses1 = v(:,:,:, n_offset);
        row = size(M0_volume_Ses1,1);
        col = size(M0_volume_Ses1,2);
        N_slice = size(M0_volume_Ses1, 3);
        matrix_size = [row, col, N_slice];
        vox_size = [220/row 220/col 5]; %   mm/pxl, slice thickness is 5 mm
        nii_image = make_nii(M0_volume_Ses1, vox_size, round(matrix_size/2), 16);
        sname = [datapath filesep special_name num2str(n_offset) '.nii'];
        save_nii(nii_image,  sname);
        gz_fname{n_offset} = sname;
end
gzip(gz_fname);

%   Read in the raw data and then run SPM12
for i = 1:69
    sourcepath = [datapath, filesep, special_name, num2str(i), '.nii,1'];
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[datapath, filesep, special_name, num2str(28), '.nii,1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {sourcepath};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'motionCorr_';
    spm_jobman('run', matlabbatch);
end

%   Load motion-corrected .nii files
%   1st frame
nii_file = [datapath filesep 'motionCorr_' special_name num2str(1) '.nii'];
nii_all = niftiread(nii_file);
%   the other frames
for i_dyn = 2: 69
    nii_file = [datapath filesep  'motionCorr_'  special_name num2str(i_dyn) '.nii'];
    nii_img = niftiread(nii_file);
    nii_all = cat(4, nii_all, nii_img);
end

imgs = double(nii_all);
IMGtmp = squeeze(imgs);
cestimgs = IMGtmp;

%% STEP2: Obtain Z map & R map

nslice = size(cestimgs, 3);
offset = 2;
FitParam.CalSNR = 1;
FitParam.ifshowimage = 0;
FitParam.PeakRange = [1, 5];
FitParam.Magfield = 42.58*3; % 3 T
FitParam.satpwr = 0.8;    % saturation power / field strength, in the unit of uT

se = strel('cube',3);
for s = 1: size(cestimgs, 3)   
    fname = [datapath, filesep, special_name, '_ZRMmap_slice', num2str(s)];
    idx = s;
    mymkdir(fname)
    smask = squeeze(Al_allS(:, :, idx));
    bg = smask<.5;
    erode_bg = imerode(bg, se);
    smask = erode_bg<.5;
    scestimgs = squeeze(cestimgs(:, :, idx, :));
    
    FitParam.SNRimg = cestimgs(:, :, idx, offset);
    FitParam.SNRmask = smask;
    FitParam.datapath = fname;
    T1_map = squeeze(T1map(:, :, idx));
    
    [displayimg_M0, FitResult] = amide_process(smask, scestimgs, fullppm, T1_map, FitParam);
   
    %%  Some results you don't want to separate across WM/GM/CSF
    
    % SNR of each slice
    snr_allS(s) = FitResult.SNR; 
    % CESTM0 imag of each slice
    cestOriginal(:, :, idx) = displayimg_M0;
    new_mask(:, :, idx) = FitResult.mask;
    
    %% STEP3 : Plot Z map and R map , then save it and store in an atlas

    % % display Zamide map
    % figure();
    % set(gca,'Position',[0.05 0.05 0.9 0.9]);
    % imagesc(displayimg_M0, [0, 2]); % 0, 4
    % title('Zamide map (%)')
    % colormap(inferno)
    % colormap(gray(255))
    % colorbar('location','Eastoutside','FontSize', 18)
    % set(gca,'dataAspectRatio', [1 1 1]);
    % axis off
    % hold off

    %   amideCEST maps
    Zalist(:, :, idx) = FitResult.ZamideMap;
    %   GuanCEST maps
    Zglist(:, :, idx) = FitResult.ZguanMap;

end
