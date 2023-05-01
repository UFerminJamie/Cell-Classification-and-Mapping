%% extract_intensity.m
% Extracts intensity features of all patches and stores them in a giant
% feature array.

% Same as extract_intensity_v2.m, except there are additional graphing
% things in here.

% In feature array, each row is one cell and every two columns correspond
% to one channel. Columns are organized as:
% [std_ch1 mean_ch1 std_ch2 mean_ch2 ... std_ch17 mean_ch17]
% create a cell array of strings

%% Load in all patches
% Load file names into array
clc, clear, close all;
% Start timer
tic

fileinfo = dir('/blue/pinaki.sarder/j.fermin/cell-mapping-data/patch_folder_reg_data');
fnames = {fileinfo.name};
fnames = {fnames{3:end}}';% Clearing '.' and '..' from structure (first 2 elements)
fname = string(fnames)';

%fnames_ = cell2table(fnames);
%writetable(fnames_, 'filenames3004_1.csv');
% adding data path
addpath('/blue/pinaki.sarder/j.fermin/cell-mapping-data/patch_folder_reg_data')

%% Clear out corrupt files
% MATLAB has issues with some of the files (not sure why) Clearing all of
% these ones out for now Could combine this with the loop later to extract
% features, would save a few minutes of run-time (not worth it right now)

fnames_t = {};
errors = [];
bad_files = {};
num_bad = 0;

for i = 1:length(fnames) % Iterate through patch files
    i
    
    try % Checking if MATLAB can successfully read the tiff file
        imread(fnames{i});
    catch error % If not, do not consider this file for analysis
        errors = [errors error]; % Storing errors
        num_bad = num_bad + 1;
        bad_files{num_bad} = fnames{i}; % Storing bad file names
        continue;
    end  
    
    % Program reaches here if the file is valid
    fnames_t{i - num_bad} = fnames{i};
end

% Update file names array
fnames = fnames_t; clear fnames_t;

% fnames_ = cell2table(fnames);
% writetable(fnames_, 'filenames_reg_data.csv')


%% Make giant array of intensity features ---- mean and SD
num_f = 2; % Features per channel -- 2 for now
num_ch = 16; % Number of channels
features_mean_sd = zeros(length(fnames), num_f*num_ch);

for j = 1:length(fnames) % Loop through patches
    fname = fnames{j};
    j
    
    ft_idx = 1;
    dapi = imread(fname,1);
    for i = 3:18 % Loop through channels
 
        try
            patch = double(imread(fname,i));
            dapim_region = double(patch) .* double(dapi > 0);
        
            %patch = patch/max(patch,[],'all'); % Scale 0 to 1
            %figure(10), imagesc(patch);
            %ch_hist = imhist(patch); % Generate intensity histogram of channel
            %features(j,ft_idx) = std(patch,[],'all'); ft_idx = ft_idx + 1;
            features_mean_sd(j,ft_idx) = std(dapim_region,[],'all'); ft_idx = ft_idx + 1;
            features_mean_sd(j,ft_idx) = mean(dapim_region,'all'); ft_idx = ft_idx + 1;
            
        catch
            warning('Error occurred on iteration %d. Skipping.', i);
            ft_idx = ft_idx + 2;
            break;
        end

    end
end

%% Make giant array of intensity features ---- mean only
num_f = 1; % Features per channel -- 2 for now
num_ch = 16; % Number of channels
features = zeros(length(fnames), num_f*num_ch);

for j = 1:length(fnames) % Loop through patches
    fname = fnames{j};
    j
    
    ft_idx = 1;
    dapi = imread(fname,1);
    for i = 3:18 % Loop through channels
        
        try
            patch = double(imread(fname,i));
            dapim_region = double(patch) .* double(dapi > 0);
        
            %patch = patch/max(patch,[],'all'); % Scale 0 to 1
            %figure(10), imagesc(patch);
            %ch_hist = imhist(patch); % Generate intensity histogram of channel
            %features(j,ft_idx) = std(patch,[],'all'); ft_idx = ft_idx + 1;
            features(j,ft_idx) = mean(dapim_region,'all'); 
            ft_idx = ft_idx + 1;
        catch
            warning('Error occurred on iteration %d. Skipping.', i);
            ft_idx = ft_idx + 1;

            break;
        end
    end
end

%% Save feature matrix and file names
% initialize a new table
T3 = table();
T2 = transpose(fnames);
% loop over the rows of the original table
for i = 1:size(T2,1)
    % split the string using "_" separator
    str_cells = split(T2(i), '_');
    
    % extract the first cell
    first_cell = str_cells(1);
    
    % add the first cell as a new row to the new table
    T3 = [T3; table(str2num(first_cell{1}))];
end

% mean features
mu = mean(features);
sigma = std(features);
data_norm = (features - mu) ./ sigma;
T_mean = array2table(data_norm);
% normalize each column using mean and SD
T_mean_norm = normalize(T_mean, 'range');
T_ = horzcat(T3, T_mean_norm);
T_ = sortrows(T_, 'Var1', 'ascend');
% Transpose



mu = mean(features_mean_sd);
sigma = std(features_mean_sd);
data_norm = (features_mean_sd - mu) ./ sigma;
T_meansd = array2table(data_norm);
% normalize each column using mean and SD
T_meansd_norm =normalize(T_meansd, 'range');
T_2 = horzcat(T3, T_meansd_norm);
T_2= sortrows(T_2, 'Var1', 'ascend');
% Transpose


% convert the table to an array
A = table2array(T_);
A = transpose(A);
% save the array to a CSV file
filename = 'features_mean_0105_new.csv';
writematrix(A, filename);

% convert the table to an array
B = table2array(T_2);
B = transpose(B);
% save the array to a CSV file
filename = 'features_mean_SD_0105_new.csv';
writematrix(B, filename);

% End timer and display elapsed time
elapsed_time = toc;
disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);
