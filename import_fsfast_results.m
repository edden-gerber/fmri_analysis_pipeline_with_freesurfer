%% Set subject parameters
subject_id = 's110';
fs_analysis_folder = fullfile('L:\Experiments\Duration Gamma\fMRI\fMRI 1 - Long Stimuli\Analysis\FreeSurfer',subject_id,'main');

%% Load and visualize anatomical data
% The ECoG/fMRI Visualization Toolbox 
% (https://github.com/edden-gerber/ecog_fmri_visualization_matlab)
% is used to visualize the 3D cortical mesh and metadata, so add it to 
% your path: 
addpath('...\ecog_fmri_visualization_matlab');

% Load anatomical brain data
brain_data_path = fullfile('...\Analysis\Anatomical',subject_id,'Matlab','brain_data.mat');
load(brain_data_path);


hem = {'right','left'};
hem_label = {'rh','lh'};
for h = 1:2
    subplot(1,2,h);
    % plot 3D cortical mesh
    plot_mesh_brain(brain_data.(['pial_' hem{h}]));
    
    % just for fun, also paint V1-V3 retinal angle data
    ang = brain_data.metadata.ret_angle.(hem_label{h})';
    ang(ang==0) = nan; % nan values are not painted on the surface
    paint_mesh(ang,0.5,1); % plot angle as a new layer (->different colormap) with 0.5 transparency
    colormap jet;
    caxis([0 180]);
    colorbar;
end


%% Read functional data - after preprocessing and before 1st level GLM
% Each preprocessing step in FSFast produces a separate nifty-format data 
% file, named as a concatenated string of labels indicating the applied
% processing steps. So for example fmcpr.down.sm4 indicates data after
% motion correction, slice timing correction (with "down" slice order),
% and 4 mm smoothing.
processed_file_prefix = 'fmcpr.down.sm4'; 

% The read_gzip_nifti function uses MRIread.m from the FreeSurfer
% Matlab function library, so add it to your path:
addpath('...\FreeSurfer Matlab Functions');

% Preprocessing is applied to each block separately, so the data needs to
% be read from each block folder and concatenated. Block folders are named
% by a 3-digit block number (001 etc). 
% This should return all the block folders (assuming there are less than
% 100 so they all begin with 0):
block_folders = dir([fs_analysis_folder '\0*']);
num_blocks = length(block_folders);
hem_label = {'rh','lh'};
for h = 1:2
        hem = hem_label{h};
        data.(hem) = [];
        for i=1:num_blocks
            fprintf('Reading block %g...\n',i);
            file_path = fullfile(block_folders(i).folder,block_folders(i).name);
            file_name = [procroc_file_prefix '.' subject_id '.' hem '.nii.gz'];
            % This will read the nifti file. If it is still in
            % gzip-compression, it will create an unzipped file first (may
            % take some time).
            block_data = read_gzip_nifti(fullfile(file_path, file_name));
            % convert to single to decrease data size:
            block_data = single(block_data);
            % Concatenate blocks (transposing so that the time/TR dimension is first, voxel is second):
            data.(hem) = [data.(hem) ; block_data'];            
        end
end

%% Read functional data - residuals after 1st level GLM
% Use this code to read the residual data resulting from GLM analysis
% performed with FSFast (e.g. to regress out nuisance variables). 

% The read_gzip_nifti function uses MRIread.m from the FreeSurfer
% Matlab function library, so add it to your path:
addpath('...\FreeSurfer Matlab Functions');

firstlevel_analysis_folder.rh = fullfile(fs_analysis_folder,'main.rh');
firstlevel_analysis_folder.lh = fullfile(fs_analysis_folder,'main.lh');

hem_label = {'rh','lh'};
for h = 1:2
        hem = hem_label{h};
        data.(hem) = [];
        residuals_block_folders = dir([firstlevel_analysis_folder.(hem) '\res\res*']);
        num_blocks = length(residuals_block_folders);

        for i=1:num_blocks
            fprintf('Reading block %g...\n',i);
            file_path = fullfile(residuals_block_folders(i).folder,residuals_block_folders(i).name);
            % This will read the nifti file. If it is still in
            % gzip-compression, it will create an unzipped file first (may
            % take some time).
            block_data = read_gzip_nifti(fullfile(file_path));
            % convert to single to decrease data size:
            block_data = single(block_data);
            % Concatenate blocks (transposing so that the time/TR dimension is first, voxel is second):
            data.(hem) = [data.(hem) ; block_data'];            
        end
        
        % The GLM residuals fluctuate around zero (since the intercept was
        % removed as a factor), and so they cannot be used to produce a "%
        % change" signal. To compute this value, we can read the "voxel
        % offset" file, which contains the mean intensity of each voxel.
        % The percent change signal can be computed by adding the fixed 
        % offset to each voxel's time series, and dividing by the mean
        % value during baseline epochs. 
        voxel_offset.(hem) = read_gzip_nifti(fullfile(firstlevel_analysis_folder.(hem_label{h}),'h-offset.nii.gz'));
end


%% Map functional data to another cortical surface
% All cortical surfaces created with FreeSurfer are registered to the
% common FreeSurfer cortical template (fsaverage). This means that 
% functional data from any individual subject can be mapped to the common
% template and thus multiple subjects' data can be directly compared on the
% voxel-level through the common template (any other brain can also serve 
% as the common template, but other surfaces will be registered to it not
% directly by through fsaverage). 
% This code section demonstrates how to convert functional data from one
% cortical surface to another using FreeSurfer's mri_surf2surf command.
% You'll need to run it on a Linux machine with FS installed. 

% FreeSurfer initialization parameters:
Fs_folder = '/usr/local/freesurfer'; % FreeSurfer folder
Fs_subjects_folder = '/home/my_user/FreeSurfer/Subjects'; % FreeSurfer subjects folder. If using a virtual machine, note 
fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];

% Step 1: mri_surf2surf takes a nifty file as input. If your functional 
% data is already saved as a nifty file, skip to step 2. Otherwise you'll
% need to save it as a nifty file first. 
% Easiest way to do it is to load an existing nifti file and modify it.
% Chose any nifti file from your dataset: 
template_nifti = MRIread('my_nifti_file.nii');
% Replace the volume data and save (do this for left and right hemisphere data):
template_nifti.vol = my_data_lh;
MRIwrite(template_nifti,'original_surface_nifti_file_lh.nii');
template_nifti.vol = my_data_rh;
MRIwrite(template_nifti,'original_surface_nifti_file_rh.nii');

% Step 2: Convert surface-mapped data to new surface
origin_surface = 'my_subject'; % This should be the subject name as it is saved in FreeSurfer's subjects folder.
target_surface = 'fsaverage'; % use fsaverage if you want to map all subjects to a single cortical template
origin_data_file_lh = 'original_surface_nifti_file_lh.nii'; % this should be the nifti file with the original functional data mapped to origin_surface (left hemisphere surface)
target_data_file_lh = 'new_surface_nifti_file_lh.nii'; % this will hold the functional data mapped to target_surface (left hemisphere surface)
origin_data_file_rh = 'original_surface_nifti_file_rh.nii'; % Right hemisphere
target_data_file_rh = 'new_surface_nifti_file_rh.nii'; % Right hemisphere
% Run FreeSurfer's mri_surf2surf command for each hemisphere:
system([fs_shell_initialize_cmd ' mri_surf2surf --srcsubject ' origin_surface ' --srcsurfval ' origin_data_file_lh ...
    ' --trgsubject ' target_surface ' --trgsurfval ' target_data_file_lh ' --hemi lh']);
system([fs_shell_initialize_cmd ' mri_surf2surf --srcsubject ' origin_surface ' --srcsurfval ' origin_data_file_rh ...
    ' --trgsubject ' target_surface ' --trgsurfval ' target_data_file_rh ' --hemi rh']);

% Step 3: load the new functional data into Matlab
new_data_lh = MRIread(target_data_file_lh); 
new_data_lh = new_data_lh.vol;
new_data_rh = MRIread(target_data_file_rh); 
new_data_rh = new_data_rh.vol;

%% Read 1st level GLM results 
% The results of the 1st level analysis are stored as nifty files in each
% hemisphere's analysis folder. Read them using MRIread (or read_gzip_nifti
% if you haven't unzipped them. 

% This is just an example of reading beta values. You can do the same with
% planned contrast t-value etc. 
betas = read_gzip_nifti(fullfile(fs_analysis_folder,'main.rh','beta.nii.gz'));


%% Read head movement data
% Head movement data is derived for each block during preprocessing
fprintf('Loading head movement data\n');
% This should return all the block folders (assuming there are less than
% 100 so they all begin with 0):
block_folders = dir([fs_analysis_folder '\0*']);
num_blocks = length(block_folders);

head_mov_data = [];
for b = 1:num_blocks
    % fmcpr.mcdat is the cryptically named file holding the head movement
    % data
    head_mov_file_name = fullfile(fs_analysis_folder,block_folders(b).name,'fmcpr.mcdat');
    % format is tab-delimited text
    d = dlmread(head_mov_file_name);
    % the first column is just a TR counter
    d = d(:,2:10);
    % concatenate blocks
    head_mov_data = [head_mov_data ; d]; 
end
% arrange the results in a table:
% (head movement data format taken from:
% https://www.mail-archive.com/freesurfer@nmr.mgh.harvard.edu/msg24594.html)
head_mov_data = table(head_mov_data(:,1),head_mov_data(:,2),head_mov_data(:,3),head_mov_data(:,4),head_mov_data(:,5),head_mov_data(:,6),head_mov_data(:,7),head_mov_data(:,8),head_mov_data(:,9),...
    'VariableNames',{'roll_deg_ccw','pitch_deg_ccw','yaw_deg_ccw','disp_superior_mm','disp_left_mm','disp_posterior_mm','rms_input_to_reference','rms_output_to_reference','total_translation_mm'});


%% Import 1st level GLM regressor structure
% If you ran 1st level analysis on your data, you can load the regressor
% matrix as it's stored as a matlab variable ("X.mat").
% The matrix will be size TxN, where T is the total number of TRs in the
% experiment, and N is the number of regressors. N=b*r where b is the
% number of blocks and r is the number of regression factors. The first
% factors should be nuisance variables as you defined them with the
% mkanalysis-sess command - polynomial regressors, head movement regressors
% etc. 
% It is not completely clear how the 3 head movement regressors are related
% to the 9 head movement variables extracted in the previous code section.
% The documentation of mkanalysis-sess mentions that the movement data is
% orthogonalized and the three top components are used, but I could not
% replicate how exactly this was done. 

firstlevel_analysis_folder.rh = fullfile(fs_analysis_folder,'main.rh');
firstlevel_analysis_folder.lh = fullfile(fs_analysis_folder,'main.lh');

% should be the same for rh and lh: 
Load = load(fullfile(firstlevel_analysis_folder.rh,'X.mat'));
regression_matrix = Load.X;
