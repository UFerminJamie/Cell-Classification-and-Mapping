%% extract_intensity.m
% Extracts intensity features of all patches and stores them in a giant
% feature array.

% Same as extract_intensity_v2.m, except there are additional graphing
% things in here.

% In feature array, each row is one cell and every two columns correspond
% to one channel. There are currently 37 features per cell (3 morphological
% derived from the DAPI mask and 2 intensity features derived from each of
% the 17 channels). Columns are organized as: 
% [mask_aspect_ratio mask_per mask_area std_ch1 mean_ch1 std_ch2 mean_ch2 ... std_ch17 mean_ch17]

%% Load in all patches
% Load file names into array
fileinfo = dir('patch_folder_10px_with_mask');
fnames = {fileinfo.name};
fnames = {fnames{3:end}}; % Clearing '.' and '..' from structure (first 2 elements)

%% Clear out corrupt files (all files seem to work now... can skip)
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

%% Animation to display all 18 channels of one patch
idx = 10001;
fname = fnames{idx};
patch = zeros([size(imread(fname)) 18]);

figure(5),
for i = 1:18
    patch(:,:,i) = imread(fname,i);
    patch(:,:,i) = patch(:,:,i)/max(patch(:,:,i),[],'all'); % Normalize
    imagesc(patch(:,:,i)), axis image, title("Channel " + i), colorbar
    pause(0.5);
end

%% Make giant array of intensity features
num_ch = 18; % Number of channels
num_f_ch = 2; % Features per channel
num_f = num_ch*(num_f_ch - 1) + 3; % Total features per cell -- 2 intensity (per channel), 3 morphological per cell
features = zeros(length(fnames), num_f);

for j = 1:length(fnames) % Loop through patches
    fname = fnames{j};
    j
        
    % GET MORPHOLOGICAL FEATURES
    % Clear out any other blobs included in each patch (this is done by
    % assuming the cell-of-interest is the largest blob)
    dapi_mask = double(imread(fname,1));
    dapi_mask = logical(dapi_mask(:,:,1));
    dapi_mask = imfill(dapi_mask,'holes'); % fill all blobs
    dapi_mask = bwareafilt(dapi_mask,1); % gets largest blob
    morph = regionprops('table', dapi_mask, 'MajorAxisLength', 'MinorAxisLength', 'Area', 'Perimeter');

    ft_idx = 1;
    features(j,ft_idx) = morph.MajorAxisLength / morph.MinorAxisLength; ft_idx = ft_idx + 1; % Aspect Ratio
    features(j,ft_idx) = morph.Perimeter; ft_idx = ft_idx + 1; % Perimeter
    features(j,ft_idx) = morph.Area; ft_idx = ft_idx + 1; % Area

    
    % GET INTENSITY FEATURES
    for i = 2:num_ch % Loop through channels
        patch = double(imread(fname,i));
        patch = patch/max(patch,[],'all'); % Scale 0 to 1
        %figure(10), imagesc(patch);
        %ch_hist = imhist(patch); % Generate intensity histogram of channel
        features(j,ft_idx) = std(patch,[],'all'); ft_idx = ft_idx + 1;
        features(j,ft_idx) = mean(patch,'all'); ft_idx = ft_idx + 1;
    end
end


%% Preliminary Clustering
k_clusters = kmeans(features(:,:),10,'MaxIter',1000);
% Will give a warning about empty entries -- some have std and mean of NaN
% for empty channels (only a few entries). Just ignore the warning.

%%
% Plot clusters
figure(200), clf
hold on
for i=1:17
    subplot(3,6,i)
    P = gscatter(features(:,2*(i-1)+1),features(:,2*i),k_clusters(:)); 
    set(P,'MarkerSize',1);
    b = gca; legend(b,'off');
end
hold off

%% Save feature matrix and file names
save('feature_mat.mat', 'features'); % Matrix of features
%save('fnames_bad.mat', 'bad_files'); % Matrix of bad file names
%save('fnames_good.mat', 'fnames'); % Matrix of working files
