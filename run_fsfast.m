% This script implements the FreeSurfer/FSFast analysis pipeline for 
% anatomical and functional MRI data. It starts with raw scans produces 
% matlab structures for anatomical cortical surface data and meta-data 
% and processed functional data (up to 1st-level analysis). 
%
% Read the guide "How to analyze fMRI data with FreeSurfer and
% FSFAST" (https://edden-gerber.github.io/analyze_fmri_data) for a full 
% explanation of the script and analysis pipeline. 
%
% Most code sections can only be run on a linux machine with FreeSurfer 
% installed. The open-source MRIcron is used for DICON->NIFTI conversion. 
%
% Folder names and parameters are initialized according to the "Sustained visual
% perception fMRI experiment 1" as an example. This experiment included a
% main slow event-related experimental session to which 1st level GLM
% analysis was applied only to regress out nuisance factors, with the
% residuals imported to matlab for subsequent processing - and a high level 
% visual category localizer session, for which the FsFast's GLM was used to
% apply the contrasts between the stimulus categories (faces, houses,
% etc.), with the resulting voxel-level t-value maps imported to Matlab. 
% You should modify the relevant parts according to your experimental 
% setup. 


%% Set global parameters

% Global parameters
Fs_folder = '/usr/local/freesurfer'; % FreeSurfer folder
Fs_subjects_folder = '/home/hcnl/FreeSurfer/Subjects'; % FreeSurfer subjects folder. 
mricron_folder = '/home/hcnl/Apps/mricron'; % MRIcron program for converting DICOMS to nifti
retinotpy_atlas_folder = '/home/hcnl/shared_folders/duration_gamma_fmri/Analysis/pRF/Atlas';

initial_analysis_folder = '/home/hcnl/shared_folders/duration_gamma_fmri/Analysis/FsFast'; 
initial_scans_folder = '/home/hcnl/shared_folders/duration_gamma_fmri/results'; 
initial_anatomy_folder = '/home/hcnl/shared_folders/duration_gamma_fmri/Analysis/Anatomical'; 

slice_order = 'down' ;  % down, up, odd, even or siemens
smoothing = '4';        % voxel smoothing in mm
hpf_freq = '0.005';     % high pass frequency

% the subject_nonstop_run parameter is set to true if this entire script
% should be run in one go, without stopping for feedback and without
% running the one-time-processing sections. Set it to true before running
% the script if you want to do this, e.g. "subject_nonstop_run=true;
% This would be useful for example if you're calling this script from a 
% loop over all subjects. 
if ~exist('subject_nonstop_run')
    subject_nonstop_run = false;
end


%% Set subject parameters
% Set parameters in code
if ~subject_nonstop_run
    subject_id = [];
end

% use input dialogs (runs if the variable is empty)
if isempty(subject_id)
    subject_id = inputdlg('Enter subject ID: ');
    subject_id = subject_id{1};
end

subject_scans_folder = fullfile(initial_scans_folder,subject_id);
subject_analysis_folder = fullfile(initial_analysis_folder,subject_id);
subject_anatomy_folder = fullfile(initial_anatomy_folder,subject_id);

% The example experiment also had a visual category localizer scan - which
% produced matlab-format metadata saved here: 
subject_category_localizer_results_folder = fullfile(initial_scans_folder,'MAT files','Localizer',subject_id);

subject_log_file = fullfile(subject_analysis_folder,'analysis_log.txt'); % logs all the output to the command window during analysis. Appended if file already exists


%% Create subject folders and subjectname file and convert T1 DICOM to nifti
if ~exist(subject_analysis_folder)
    mkdir(subject_analysis_folder);
end

diary(subject_log_file);
fprintf('\t\tCREATING SUBJECT FOLDERS\n');

% Add FreeSurfer subjectname file
subjectname_filename = fullfile(subject_analysis_folder,'subjectname');
file_handle = fopen(subjectname_filename,'wt');
fprintf(file_handle,subject_id);
fclose(file_handle);

% Create anatomy output folder if it doesn't exist
if ~exist(subject_anatomy_folder,'dir')
    mkdir(subject_anatomy_folder);
end

% Create Matlab output folder if it doesn't exist
if ~exist(fullfile(subject_anatomy_folder,'Matlab'),'dir')
    mkdir(fullfile(subject_anatomy_folder,'Matlab'));
end


% Create images output folder if it doesn't exist
if ~exist(fullfile(subject_anatomy_folder,'Generated_Images'),'dir')
    mkdir(fullfile(subject_anatomy_folder,'Generated_Images'));
end

% Convert anatomy DICOM files to nifti
fprintf('\t\tCONVERTING DICOMS TO NIFTI\n');

t1_scan_dir = dir(fullfile(subject_scans_folder,'*MPRAGE*')); % find anatomy scan folder
t1_scan_dir = fullfile(subject_scans_folder,t1_scan_dir.name);
system([mricron_folder '/dcm2nii -o ' subject_anatomy_folder ' ' t1_scan_dir]);

% from the three output files, keep only the final cropped and re-oriented 
% file and rename it T1.nii.gz:
files = dir(fullfile(subject_anatomy_folder,'*.nii.gz'));
for f=1:length(files)
    fname = fullfile(subject_anatomy_folder,files(f).name);
    if strcmp(files(f).name(1:2),'co')
        movefile(fname,fullfile(subject_anatomy_folder,'T1.nii.gz'));
    else
        delete(fname);
    end
end

diary off


%% Run FreeSurfer on the anatomical scan

% This step creates the segmented volume and surface data for the selected
% brain, and applied surface co-registration to the fsaverage template
% brain. 
% NOTE THAT THIS STEP CAN TAKE UP TO 24 HOURS

diary(subject_log_file);
fprintf('\t\tRUNNING FREESURFER ANATOMICAL ANALYSIS\n');

fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
disp(['Running FreeSurfer for subject ' subject_id '...']);
system([fs_shell_initialize_cmd 'recon-all -all -i ' fullfile(subject_anatomy_folder,'T1.nii.gz') ' -s ' subject_id]);

diary off


%% Inspect the results using freeview
% Uses FreeSurfer's freeview utility. 
% If there are errors in the reconstruction (most likely regions where Fs
% fails to identify white matter, resulting in an erroneous pial surface
% shape), use freeview to correct the errors and re-run the appropriate
% segments of the FreeSurfer analysis. Consult FS's wiki and YouTube videos
% for assistance. 

if ~subject_nonstop_run
    fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
    system([fs_shell_initialize_cmd 'freeview -v "' Fs_subjects_folder '/' subject_id '/mri/T1.mgz" ' ...
        '-f "' Fs_subjects_folder '/' subject_id '/surf/lh.pial" ' '-f "' Fs_subjects_folder '/' subject_id '/surf/rh.pial"']);
end


%% Unpack functional files into FsFast directory structure

% Map existing DICOM files (this takes hours for some reason). 
diary(subject_log_file);
fprintf('\t\tUNPACKING FUNCTIONAL DATA\n');

index_filename = fullfile(subject_analysis_folder,'dcm.index.dat');
system(['dcmunpack -src ' subject_scans_folder ' -index-out ' index_filename]);

% and extract the files - define the folder numbers based on the actual
% scanning order (can change between subjects, for example if the localizer
% is run again after a bathroom break). It is therefore recommended to use
% the "switch/case" selection structure to define this information for each
% subject separately. 

% *** Modify this example section according to your experiment's structure *** 
switch subject_id
    case 's101'
        % Unpack and convert main experiment scans
        system(['dcmunpack -src ' subject_scans_folder ' -index-in ' index_filename ' -targ ' subject_analysis_folder ...
            ' -run 3 main nii f.nii -run 4 main nii f.nii -run 5 main nii f.nii -run 6 main nii f.nii ']);
        % Unpack and convert visual catgegory localizer scans
        system(['dcmunpack -src ' subject_scans_folder ' -index-in ' index_filename ' -targ ' subject_analysis_folder ...
            ' -run 11 catLoc nii f.nii -run 12 catLoc nii f.nii -run 13 catLoc nii f.nii ' ]);
    case 's102'
        % Unpack and convert main experiment scans
        system(['dcmunpack -src ' subject_scans_folder ' -index-in ' index_filename ' -targ ' subject_analysis_folder ...
            ' -run 4 main nii f.nii -run 5 main nii f.nii -run 7 main nii f.nii -run 8 main nii f.nii ']);
        % Unpack and convert visual catgegory localizer scans
        system(['dcmunpack -src ' subject_scans_folder ' -index-in ' index_filename ' -targ ' subject_analysis_folder ...
            ' -run 12 catLoc nii f.nii -run 13 catLoc nii f.nii -run 14 catLoc nii f.nii ' ]);
    otherwise
        error('subject scan folders not defined!.');

end
diary off


%% Import FreeSurfer-generated images and convert to nifti, produce stripped cortex volume images
% This section imports volume images produced by FreeSurfer from the FS
% subjects folder into the anatomy folder, after converting them to Nifti
% format. It also generates new volume images of the left and right
% cortical hemisphere. 
% 
% The following images are created: 
% MR_normalized: This is the original MR scan, after normalization to
%                FreeSurfer's 256x256x256 volume format. 
% MR_T1: this is the same volume after intensity normalization - all voxels
%                of the same tissue type will have the same intensity.
% cortex / cortex_right / cortex_left: These are images of the left, right
%                or both cortical hemispheres after the head, skull and 
%                subcortical structures were stripped away. 

diary(subject_log_file);
fprintf('\t\tCONVERTING FREESURFER-GENERATED VOLUMES TO NIFTI\n');

output_images_folder = fullfile(subject_anatomy_folder,'Generated_Images');

% Create folder in the FreeSurfer subject directory for temporary files
fs_temp_folder = fullfile(Fs_subjects_folder,subject_id,'mri','temp');
if ~exist(fs_temp_folder,'dir')
    mkdir(fs_temp_folder);
end

% Convert the resliced and normalized MR image to nifty format (using a FS 
% function) and save it in the results folder
fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
path_mr_norm_mgz = fullfile(Fs_subjects_folder,subject_id,'/mri/nu.mgz');
path_mr_norm = fullfile(output_images_folder,[subject_id '_MR_normalized.nii.gz']);
system([fs_shell_initialize_cmd 'mri_convert ' path_mr_norm_mgz ' ' path_mr_norm]);

% Do the same for the FreeSurfer-generated T1 image (similar to the 
% normalized image but after intensity "standardization", e.g. all white 
% matter voxels are set to 110, etc.). 
path_mr_t1_mgz = fullfile(Fs_subjects_folder,subject_id,'/mri/T1.mgz');
path_mr_t1 = fullfile(output_images_folder,[subject_id '_MR_T1.nii.gz']);
system([fs_shell_initialize_cmd 'mri_convert ' path_mr_t1_mgz ' ' path_mr_t1]);

% Generate cortex image by masking the T1 image with the cortical ribbon mask
path_ribbon = fullfile(Fs_subjects_folder, subject_id, 'mri/ribbon.mgz');
system([fs_shell_initialize_cmd 'mri_mask "' path_mr_t1_mgz '" "' path_ribbon '" "' fs_temp_folder '/cortex.mgz"']);

% Generate left/right hemisphere cortex volumes (first creating masks for them):
system([fs_shell_initialize_cmd 'mri_binarize --i "' path_ribbon '" ' ...
    '--o "' fs_temp_folder '/ribbon_mask_left.mgz" --min 2 --max 3']);
system([fs_shell_initialize_cmd 'mri_binarize --i "' path_ribbon '" ' ...
    '--o "' fs_temp_folder '/ribbon_mask_right.mgz" --min 41 --max 42']);
system([fs_shell_initialize_cmd 'mri_mask "' path_mr_t1_mgz '" ' ...
    '"' fs_temp_folder '/ribbon_mask_left.mgz" "' fs_temp_folder '/cortex_left.mgz"']);
system([fs_shell_initialize_cmd 'mri_mask "' path_mr_t1_mgz '" ' ...
    '"' fs_temp_folder '/ribbon_mask_right.mgz" "' fs_temp_folder '/cortex_right.mgz"']);

% Convert all scans to Nifti and copy them to the output directory 
system([fs_shell_initialize_cmd 'mri_convert "' fs_temp_folder '/cortex.mgz" "' output_images_folder '/' subject_id '_cortex.nii.gz"']);
system([fs_shell_initialize_cmd 'mri_convert "' fs_temp_folder '/cortex_left.mgz" "' output_images_folder '/' subject_id '_cortex_left.nii.gz"']);
system([fs_shell_initialize_cmd 'mri_convert "' fs_temp_folder '/cortex_right.mgz" "' output_images_folder '/' subject_id '_cortex_right.nii.gz"']);

disp('Finished');
diary off


%% Read pial and inflated cortical surfaces from FreeSurfer format to Matlab
% In this section, the pial and inflated surfaces generated by FreeSurfer
% are read into a Matlab data structure (using FS's read_surf function).
% The surfaces are transformed from the native FS coordinate space to the
% standard RAS anatomical space. The result is a data structure with the
% following fields:
% - pial_left/pial_right: pial (non-inflated) cortical surfaces
% inflated_left/inflated_right: inflated cortical surfaces
%
% Note that since the inflated surfaces are both initially centered around
% the origin, they would overlap if plotted together. Therefore in an
% adidtional step the two inflated hemispheres are pulled apart (this has
% no impact on plotting electrodes on them since for inflated surfaces, the 
% electrodes are mapped to a specific mesh vertex and not to an xyz 
% coordinate. 
% 
% This section is partially based on code by Tal Golan @ Rafael Malach lab.

diary(subject_log_file);
fprintf('\t\tIMPORTING FREESURFER-GENERATED SURFACES INTO MATLAB\n');

% Get transformation matrices: 
disp('Getting transformation matrices...');
fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
[~, cmdout] = system([fs_shell_initialize_cmd 'mri_info --vox2ras ' Fs_subjects_folder '/' subject_id '/mri/orig.mgz']);
transformations.ijk2xyz = affine3d(str2num(cmdout)');

[~, cmdout] = system([fs_shell_initialize_cmd 'mri_info --vox2ras-tkr ' Fs_subjects_folder '/' subject_id '/mri/orig.mgz']);
transformations.ijk2xyz_FsMesh = affine3d(str2num(cmdout)');

brain_data = struct;

% Read pial surfaces to Matlab, transform and save: 
disp('Reading pial surfaces...');
% Right hemisphere
[brain_data.pial_right.vertices,brain_data.pial_right.faces] = read_surf(fullfile(Fs_subjects_folder,subject_id,'surf','rh.pial'));
brain_data.pial_right.faces = brain_data.pial_right.faces + 1;
brain_data.pial_right.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.pial_right.vertices);
brain_data.pial_right.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.pial_right.vertices);
% Left hemisphere
[brain_data.pial_left.vertices,brain_data.pial_left.faces] = read_surf(fullfile(Fs_subjects_folder,subject_id,'surf','lh.pial'));
brain_data.pial_left.faces = brain_data.pial_left.faces + 1;
brain_data.pial_left.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.pial_left.vertices);
brain_data.pial_left.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.pial_left.vertices);

% Read inflated surfaces: 
disp('Reading inflated surfaces...');
% Right hemisphere
[brain_data.inflated_right.vertices,brain_data.inflated_right.faces] = read_surf(fullfile(Fs_subjects_folder,subject_id,'surf','rh.inflated'));
brain_data.inflated_right.faces = brain_data.inflated_right.faces + 1;
brain_data.inflated_right.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.inflated_right.vertices);
brain_data.inflated_right.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.inflated_right.vertices);
% Left hemisphere
[brain_data.inflated_left.vertices,brain_data.inflated_left.faces] = read_surf(fullfile(Fs_subjects_folder,subject_id,'surf','lh.inflated'));
brain_data.inflated_left.faces = brain_data.inflated_left.faces + 1;
brain_data.inflated_left.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.inflated_left.vertices);
brain_data.inflated_left.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.inflated_left.vertices);
% Pull apart the two inflated hemispheres so they do not overlap if plotted together:
[brain_data.inflated_left, brain_data.inflated_right] = pull_apart_inflated_hemispheres(brain_data.inflated_left, brain_data.inflated_right, 10);

% Read curvature
curv_lh = read_curv(fullfile(Fs_subjects_folder,subject_id,'surf','lh.curv'));
curv_rh = read_curv(fullfile(Fs_subjects_folder,subject_id,'surf','rh.curv'));
brain_data.metadata.curvature.lh = curv_lh;
brain_data.metadata.curvature.rh = curv_rh;

% Read thickness
thick_lh = read_curv(fullfile(Fs_subjects_folder,subject_id,'surf','lh.thickness'));
thick_rh = read_curv(fullfile(Fs_subjects_folder,subject_id,'surf','rh.curv'));
brain_data.metadata.thickness.lh = thick_lh;
brain_data.metadata.thickness.rh = thick_rh;

% Read parcellation
brain_data.metadata.parcellation.lh(1:size(brain_data.pial_left.vertices,1)) = nan;
brain_data.metadata.parcellation.rh(1:size(brain_data.pial_right.vertices,1)) = nan;
% Left
[vertices,label,colortable]=read_annotation(fullfile(Fs_subjects_folder,subject_id,'label','lh.aparc.a2009s.annot'));
for p = 1:colortable.numEntries
    v_idx = find(label == colortable.table(p,5));
    brain_data.metadata.parcellation.lh(v_idx) = p;
end
% Right
[vertices,label,colortable]=read_annotation(fullfile(Fs_subjects_folder,subject_id,'label','rh.aparc.a2009s.annot'));
for p = 1:colortable.numEntries
    v_idx = find(label == colortable.table(p,5));
    brain_data.metadata.parcellation.rh(v_idx) = p;
end
brain_data.metadata.parcellation.labels = colortable.struct_names;

% Save
output_path = fullfile(subject_anatomy_folder,'Matlab','brain_data'); % freesurfer mesh data
disp(['Saving brain mesh data to ' output_path ]);
save(output_path,'brain_data');

disp('Done.');
diary off


%% Add retinotopic maps from external V1/V2/V3 atlas
% This section is for working with the V1/V2/V3 atlas described in:
% https://cfn.upenn.edu/aguirre/wiki/public:data_ploscomputbiol_2014_benson
%
% Make sure you cite: NC Benson, OH Butt, R Datta, PD Radoeva, DH Brainard, 
% Benson, N. C., Butt, O. H., Brainard, D. H., & Aguirre, G. K. (2014). 
% Correction of Distortion in Flattened Representations of the Cortical 
% Surface Allows Prediction of V1-V3 Functional Organization from 
% Anatomy. PLoS Computational Biology, 10(3).
% AND/OR:
% GK Aguirre (2012) The retinotopic organization of striate cortex is well 
% predicted by surface topology. Current Biology. 
%
% To work with the atlas, first download the three .mgh files from:
% https://cfn.upenn.edu/aguirre/wiki/public:data_ploscomputbiol_2014_benson
% and put them in a folder as defined in the global definitions of this
% script (retinotpy_atlas_folder). 
%
% Also make sure you have the fsaverage_sym template, which is the surface 
% template the original data is mapped to. If not you can download it from:
% ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/5.1.0/xhemi/fsaverage_sym.tar.gz
% extract it and put it in the FreeSurfer subjects folder. 
%

diary(subject_log_file);
fprintf('\t\tADDING RETINOTOPIC ATLAS METADATA\n');

fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];

% Load subject's brain_data struct:
load(fullfile(subject_anatomy_folder,'Matlab','brain_data'));

system([fs_shell_initialize_cmd ' surfreg --s ' subject_id ' --t  fsaverage_sym  --lh']);
system([fs_shell_initialize_cmd ' surfreg --s ' subject_id ' --t  fsaverage_sym  --lh --xhemi']);


% Register the retinotopy data from the fsaverage_sym template surface to 
% the subject's cortical surface: 
disp('Retinotopy: Angle');
src = fullfile(retinotpy_atlas_folder,'angle-template.sym.mgh');
system([fs_shell_initialize_cmd ' mri_surf2surf --srcsubject fsaverage_sym --trgsubject ' subject_id ...
    ' --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --srcsurfval ' src ' --tval temp.mgh']);
f = MRIread('temp.mgh');
brain_data.metadata.ret_angle.lh = squeeze(f.vol);
system([fs_shell_initialize_cmd ' mri_surf2surf --srcsubject fsaverage_sym --trgsubject ' subject_id '/xhemi'...
    ' --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --srcsurfval ' src ' --tval temp.mgh']);
f = MRIread('temp.mgh');
brain_data.metadata.ret_angle.rh = squeeze(f.vol);

disp('Retinotopy: Eccentricity');
src = fullfile(retinotpy_atlas_folder,'eccen-template.sym.mgh');
system([fs_shell_initialize_cmd ' mri_surf2surf --srcsubject fsaverage_sym --trgsubject ' subject_id ...
    ' --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --srcsurfval ' src ' --tval temp.mgh']);
f = MRIread('temp.mgh');
brain_data.metadata.ret_eccentricity.lh = squeeze(f.vol);
system([fs_shell_initialize_cmd ' mri_surf2surf --srcsubject fsaverage_sym --trgsubject ' subject_id '/xhemi'...
    ' --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --srcsurfval ' src ' --tval temp.mgh']);
f = MRIread('temp.mgh');
brain_data.metadata.ret_eccentricity.rh = squeeze(f.vol);

disp('Retinotopy: Area');
src = fullfile(retinotpy_atlas_folder,'areas-template.sym.mgh');
system([fs_shell_initialize_cmd ' mri_surf2surf --srcsubject fsaverage_sym --trgsubject ' subject_id ...
    ' --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --srcsurfval ' src ' --tval temp.mgh']);
f = MRIread('temp.mgh');
brain_data.metadata.ret_area.lh = squeeze(f.vol);
system([fs_shell_initialize_cmd ' mri_surf2surf --srcsubject fsaverage_sym --trgsubject ' subject_id '/xhemi'...
    ' --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --srcsurfval ' src ' --tval temp.mgh']);
f = MRIread('temp.mgh');
brain_data.metadata.ret_area.rh = squeeze(f.vol);

% Delete temp file
delete('temp.mgh');

% Save
output_path = fullfile(subject_anatomy_folder,'Matlab','brain_data'); % freesurfer mesh data
disp(['Saving brain mesh data to ' output_path ]);
save(output_path,'brain_data');

disp('Done.');
diary off


%% Run pre-processing pipeline

diary(subject_log_file);
% If you have different parts to your experiment, define a preprocessing command for each one: 
% main exp
fprintf('\t\tPREPROCESSING: MAIN EXP\n');
system(['preproc-sess -d ' initial_analysis_folder ' -s ' subject_id ' -fsd main -surface ' subject_id ' lhrh -mni305 -fwhm ' smoothing ' -per-run -sliceorder ' slice_order ]);

% category localizer
fprintf('\t\tPREPROCESSING: CATEGORY LOCALIZER\n');
system(['preproc-sess -d ' initial_analysis_folder ' -s ' subject_id ' -fsd catLoc -surface ' subject_id ' lhrh -mni305 -fwhm ' smoothing ' -per-run -sliceorder ' slice_order ]);

diary off


%% Create paradigm files
% This step is needed when applying FSFast's 1st-level GLM analysis to your
% experimental paradigm. In this example experiment, this is relevant only
% for the visual category localizer session. 
% Due to the nature of this step, the code is very experiment-specific and
% included here only as an example. 

diary(subject_log_file);
fprintf('\t\tCREATING CATEGORY LOCALIZER PARADIGM FILES\n');

catLoc_out_folder = fullfile(subject_analysis_folder,'catLoc');

in_files = dir(fullfile(subject_category_localizer_results_folder,'*.mat'));
out_dir = dir(fullfile(catLoc_out_folder,'0*'));
for f=1:3
    Load = load(fullfile(subject_category_localizer_results_folder,in_files(f).name));
    
    cond = Load.theSubject.trials.cond';
    t = (0:0.5:299.5)';
    
    f_id = fopen(fullfile(catLoc_out_folder,out_dir(f).name,'paradigm_file.par'),'w+');
    for l = 1:600
            fprintf(f_id, '%.1f\t%d\t%.1f\t%d\n',t(l),cond(l),0.5,1);
    end
    fclose(f_id);
end

diary off


%% Create analysis configuration
% According to the YouTube lecture on FsFast, this step is intended to be
% run once for the entire experiment (that is, not for every subject).
% However, this step also defines the cortical surface to which the
% functional data should be mapped; I assume that the original intent was
% that all subjects will be mapped to the common template cortical surface,
% which makes it easier to run 2nd-level statistics. However in the example
% experiment shown here each brain is meant to be analyzed separately
% (based on subject-level ROI definition), and therefore the analysis
% should be applied on individual subjects' cortical surfaces. Because of
% this, this step is re-run for each subject in this analysis (overwriting
% the exisiting files using -force). 

diary(subject_log_file);
fprintf('\t\tCREATING ANALYSIS CONFIGURATIONS: MAIN EXP\n');
% Main exp lh
system(['mkanalysis-sess -force -analysis main.lh -fsd main -surface ' subject_id ' lh -fwhm ' smoothing ' -stc ' slice_order ' -notask -spmhrf 0 -refeventdur 12 -polyfit 2 -hpf ' hpf_freq ' -mcextreg -TR 1.5  -per-run']);

% Main exp rh
system(['mkanalysis-sess -force -analysis main.rh -fsd main -surface ' subject_id ' rh -fwhm ' smoothing ' -stc ' slice_order ' -notask -spmhrf 0 -refeventdur 12 -polyfit 2 -hpf ' hpf_freq ' -mcextreg -TR 1.5  -per-run']);

% Main exp mni
system(['mkanalysis-sess -force -analysis main.mni -fsd main -mni305 -fwhm ' smoothing ' -stc ' slice_order ' -notask -spmhrf 0 -refeventdur 12 -polyfit 2 -hpf ' hpf_freq ' -mcextreg -TR 1.5  -per-run']);


fprintf('\t\tCREATING ANALYSIS CONFIGURATIONS: CATEGORY LOCALIZER\n');
% category localizer lh
system(['mkanalysis-sess -force -analysis catLoc.lh -fsd catLoc -surface ' subject_id ' lh -fwhm ' smoothing ' -stc ' slice_order ' -paradigm paradigm_file.par -event-related -spmhrf 0 -refeventdur 0.5 -polyfit 2 -hpf ' hpf_freq ' -delay 9 -mcextreg -TR 1.5 -nconditions 10 -per-run']);

% category localizer rh
system(['mkanalysis-sess -force -analysis catLoc.rh -fsd catLoc -surface ' subject_id ' rh -fwhm ' smoothing ' -stc ' slice_order ' -paradigm paradigm_file.par -event-related -spmhrf 0 -refeventdur 0.5 -polyfit 2 -hpf ' hpf_freq ' -delay 9 -mcextreg -TR 1.5 -nconditions 10 -per-run']);

% category localizer mni
system(['mkanalysis-sess -force -analysis catLoc.mni -fsd catLoc -mni305 -fwhm ' smoothing ' -stc ' slice_order ' -paradigm paradigm_file.par -event-related -spmhrf 0 -refeventdur 0.5 -polyfit 2 -hpf ' hpf_freq ' -delay 9 -mcextreg -TR 1.5 -nconditions 10 -per-run']);

diary off


%% Create contrasts for 1st level analysis
% This needs to run just once. Again, this is very paradigm-specific and
% included only as an example. 

diary(subject_log_file);

if ~subject_nonstop_run
    fprintf('\t\tCREATING 1ST LEVEL CONTRASTS FOR CATEGORY LOCALIZER ANALYSIS\n');
    
    analysis_names = {'lh','rh','mni'};
    % category localizer lh
    for i = 1:3
        str = analysis_names{i};
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast letters_bodyparts -a 1 -a 2 -c 3 -c 4']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast letters_faces -a 1 -a 2 -c 5 -c 6']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast letters_places -a 1 -a 2 -c 7 -c 8']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast letters_objects -a 1 -a 2 -c 9 -c 10']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast bodyparts_faces -a 3 -a 4 -c 5 -c 6']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast bodyparts_places -a 3 -a 4 -c 7 -c 8']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast bodyparts_objects -a 3 -a 4 -c 9 -c 10']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast faces_places -a 5 -a 6 -c 7 -c 8']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast faces_objects -a 5 -a 6 -c 9 -c 10']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast places_objects -a 7 -a 8 -c 9 -c 10']);

        % One-vs-all-others contrasts:
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast letters -a 1 -a 2 -c 3 -c 4 -c 5 -c 6 -c 7 -c 8 -c 9 -c 10']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast bodyparts -a 3 -a 4 -c 1 -c 2 -c 5 -c 6 -c 7 -c 8 -c 9 -c 10']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast faces -a 5 -a 6 -c 1 -c 2 -c 3 -c 4 -c 7 -c 8 -c 9 -c 10']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast places -a 7 -a 8 -c 1 -c 2 -c 3 -c 4 -c 5 -c 6 -c 9 -c 10']);
        system(['mkcontrast-sess -analysis catLoc.' str ' -contrast objects -a 9 -a 10 -c 1 -c 2 -c 3 -c 4 -c 5 -c 6 -c 7 -c 8']);
    end
end

diary off


%% Run first level analysis
% Runs the GLM analysis (and saves the residuals if requested). 

diary(subject_log_file);
fprintf('\t\tRUNNING 1ST LEVEL ANALYSIS: MAIN EXP\n');
% Main lh
system(['selxavg3-sess -d ' initial_analysis_folder ' -s ' subject_id ' -analysis main.lh -no-con-ok -svres' ]);

% Main rh
system(['selxavg3-sess -d ' initial_analysis_folder ' -s ' subject_id ' -analysis main.rh -no-con-ok -svres' ]);

% Main mni
system(['selxavg3-sess -d ' initial_analysis_folder ' -s ' subject_id ' -analysis main.mni -no-con-ok -svres' ]);


fprintf('\t\tRUNNING 1ST LEVEL ANALYSIS: CATEGORY LOCALIZER\n');
% catLoc lh
system(['selxavg3-sess -d ' initial_analysis_folder ' -s ' subject_id ' -analysis catLoc.lh' ]); % don't keep residuals

% catLoc rh
system(['selxavg3-sess -d ' initial_analysis_folder ' -s ' subject_id ' -analysis catLoc.rh' ]); % don't keep residuals

% catLoc MNI
system(['selxavg3-sess -d ' initial_analysis_folder ' -s ' subject_id ' -analysis catLoc.mni' ]); % don't keep residuals

diary off


%% Extract zipped residual data files
% Runs gzip to convert the residual data files from .nii.gz to uncompressed 
% .nii format You can modify this to include any other compressed nifty 
% file you want to read later into Matlab (e.g. pre-processed functional 
% data files before GLM, t-values for GLM contrasts, etc.) 
diary(subject_log_file);
fprintf('\t\tEXTACTING NII.GZ RESIDUAL DATA FILES\n');

system(['gunzip ' fullfile(subject_analysis_folder,'main','main.lh','res','*.gz')]);
system(['gunzip ' fullfile(subject_analysis_folder,'main','main.rh','res','*.gz')]);

diary off
