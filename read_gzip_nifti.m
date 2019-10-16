function [ data, header ] = read_gzip_nifti( file_path )
% Unzipps a gzipped nifti data file (*.nii.gz) and reads the data. Creates
% an unzipped file if it doesn't already exists. If input if an unzipped
% file it just reads it. The first output is the volume data extracted from
% the nifti structure (.vol), the second is the entire structure. 
% NOTE: unzipping may perhaps be faster using an external application. 

if strcmp(file_path(end-2:end),'.gz')
    gz = true;
    if strcmp(file_path(end-6:end-3),'.nii')
        ni = true;
    else
        ni = false;
    end
else
    if strcmp(file_path(end-3:end),'.nii')
        gz = false;
        ni = true;
    end
end

if ~ni
    error('input path should be to a .nii or .nii.gz file');
end

if gz
    unzipped_file_path = file_path(1:end-3);
    if ~exist(unzipped_file_path,'file')
        % unzip file
        fprintf('unzipping gz file... (this may be slow, consider unzipping outside Matlab)');
        gunzip(file_path);
        fprintf(repmat('\b',1,20));
    end
else
    unzipped_file_path = file_path;
end

% Read nifti
nifti = MRIread(unzipped_file_path);

% Extract data
data = squeeze(nifti.vol);
header = nifti;

fprintf('\n');
end