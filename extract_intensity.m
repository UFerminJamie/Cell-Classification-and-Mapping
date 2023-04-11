%% extract_intensity.m
% Extracts intensity features of all patches and stores them in a giant
% feature array.

% Same as extract_intensity_v2.m, except there are additional graphing
% things in here.

% In feature array, each row is one cell and every two columns correspond
% to one channel. Columns are organized as:
% [std_ch1 mean_ch1 std_ch2 mean_ch2 ... std_ch17 mean_ch17]

%% Load in all patches
% Load file names into array
fileinfo = dir('patch_folder_10px');
fnames = {fileinfo.name};
fnames = {fnames{3:end}}; % Clearing '.' and '..' from structure (first 2 elements)

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

%% Animation to display all 17 channels of one patch
idx = 10000;
fname = fnames{idx};
patch = zeros([size(imread(fname)) 17]);

figure(5),
for i = 1:17
    patch(:,:,i) = imread(fname,i);
    patch(:,:,i) = patch(:,:,i)/max(patch(:,:,i),[],'all'); % Normalize
    imagesc(patch(:,:,i)), axis image, title("Channel " + i), colorbar
    pause(0.5);
end

%% Make giant array of intensity features
num_f = 2; % Features per channel -- 2 for now
num_ch = 17; % Number of channels
features = zeros(length(fnames), num_f*num_ch);

for j = 1:length(fnames) % Loop through patches
    fname = fnames{j};
    j
    
    ft_idx = 1;
    for i = 1:17 % Loop through channels
        patch = double(imread(fname,i));
        patch = patch/max(patch,[],'all'); % Scale 0 to 1
        %figure(10), imagesc(patch);
        %ch_hist = imhist(patch); % Generate intensity histogram of channel
        features(j,ft_idx) = std(patch,[],'all'); ft_idx = ft_idx + 1;
        features(j,ft_idx) = mean(patch,'all'); ft_idx = ft_idx + 1;
    end
end

%% Plot standard deviations versus means
i = 1:2:34; % odd indices = standard deviation
j = 2:2:34; % even indices = mean
std_all = features(1:500,i); 
std_all = std_all(:);
means_all = features(1:500,j); 
means_all = means_all(:);

% Plotting features for visualization (color = channel)
figure(100), clf
hold on
for ii=1:17 % each color is one channel
    scatter(std_all(ii:17:17*5), means_all(ii:17:17*5), 25, 'filled')
end
hold off

% Plotting features for visualization (color = cell)
figure(150), clf
hold on
for ii=1:5 % each color is one cell (plotting 50)
    scatter(std_all((ii-1)*17+1:ii*17), means_all((ii-1)*17+1:ii*17), 25, 'filled')
end
hold off

%% Preliminary Clustering
k_clusters = kmeans(features(:,:),10,'MaxIter',1000);
% will give a warning about empty entries -- some have std and mean of NaN
% for empty channels (only a few entries)

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
save('fnames_bad.mat', 'bad_files'); % Matrix of bad file names
save('fnames_good.mat', 'fnames'); % Matrix of working files
