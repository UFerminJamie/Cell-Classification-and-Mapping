umap_euclidean = csvread('UMAP_clustering_correlation_nn15_mindist0.003.csv');
umap_correlation = csvread('UMAP_clustering_correlation_nn15_mindist0.003.csv');
figure, scatter(umap_correlation(:,1),umap_correlation(:,2),10,'filled');
figure, scatter(umap_euclidean(:,1),umap_euclidean(:,2),10,'filled');

%% KMeans Euclidean

figure,
for j = 1:20
    k_clusters = kmeans(umap_euclidean(:,:),j,'MaxIter',1000);

    % Save clustering results
    fname = "k_means_euc_" + j + "_clusters";
    save(fname, 'k_clusters');
    
    P = gscatter(umap_euclidean(:,1),umap_euclidean(:,2),k_clusters(:)); 
    set(P,'MarkerSize',5);
    b = gca; legend(b,'off');
    title("Clusters = " + j);
    pause(1)
end

%% KMeans Correlation

figure,
for j = 1:20
    k_clusters = kmeans(umap_correlation(:,:),j,'MaxIter',1000);

    % Save clustering results
    fname = "k_means_corr_" + j + "_clusters";
    save(fname, 'k_clusters');
    
    P = gscatter(umap_correlation(:,1),umap_correlation(:,2),k_clusters(:)); 
    set(P,'MarkerSize',5);
    b = gca; legend(b,'off');
    title("Clusters = " + j);
    pause(1)
end

%% DBSCAN (not doing well ): )
idx = dbscan(umap_euclidean,0.1,50); 
figure(15),
P = gscatter(umap_euclidean(:,1),umap_euclidean(:,2),idx(:)); 
set(P,'MarkerSize',5); b = gca; legend(b,'off');
title('Clustering with DBSCAN (UMAP Euclidean)');

