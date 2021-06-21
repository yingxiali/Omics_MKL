
 %% load data
clc
clear all
cd C:\Users\yingxiali\Desktop\mRMR_MKL
load('C:\Users\yingxiali\Desktop\mRMR_MKL\data\all_data.mat');
 %%
tic
Start =20;
step = 10;
End = 500;
auc = [];
best_auc = [];
best_N = [];
positive_count = sum(clinical + 1) / 2;
positive_indc = crossvalind('Kfold', positive_count, 10);
negative_indc = crossvalind('Kfold', length(clinical) - positive_count,10);
[~, sort_indc] = sort(clinical);
indcs = [negative_indc; positive_indc];
indcs(sort_indc) = indcs;
[~, pos] = sort(indcs);
results = zeros(size(clinical));
rowDist = cell(10, 1);
%%
%poly & gaussian
C = 300;
kernelt={'gaussian' 'gaussian' 'gaussian' 'gaussian' 'gaussian' 'poly' 'poly' 'poly' 'poly' 'poly'};
kerneloptionvect = { [0.001 0.002 0.005 0.01 0.05 0.1 0.25 0.5 1 2 5 7 10 12 15 17 20] ...
           [0.001 0.002 0.005 0.01 0.05 0.1 0.25 0.5 1 2 5 7 10 12 15 17 20] [0.001 0.002 0.005 0.01 0.05 0.1 0.25 0.5 1 2 5 7 10 12 15 17 20] [0.001 0.002 0.005 0.01 0.05 0.1 0.25 0.5 1 2 5 7 10 12 15 17 20 ] [0.001 0.002 0.005 0.01 0.05 0.1 0.25 0.5 1 2 5 7 10 12 15 17 20] [1 2 3] [1 2 3] [1 2 3] [1 2 3] [1 2 3]};
variablevec = {'random' 'random' 'random' 'random' 'random' 'random' 'random' 'random' 'random' 'random'  };

%%
C = 300;
for ratio=Start:step:End
    ratio
    [data, feature_name,ranks,rank,feature_indc] = SelectFeature_rank( all_data, clinical, ratio, names );
    ranks = [ranks ranks];
% mkl cv
   parfor i=1:10
    [rowDist{i}, ~,variableveccell] = ...
        mklclass_random(data(indcs~=i,:), clinical(indcs~=i), ...
                    data(indcs==i,:), clinical(indcs==i), ...
                    kernelt,kerneloptionvect,variablevec, ranks, C); %#ok
end
result = cell2mat(rowDist);
result(pos) = result;

%experiment.tcga = fastAUC(class == 1, result, 1, 'tcga', 0);
[auc(ratio),fpr,tpr] = fastAUC((clinical+1)/2==1, result, 0, strcat(num2str(ratio), '_mkl_big_h'));
end
toc
%%
AUC = sort(auc);
auc(2,:) = 1:1:End;
idx = any(auc==0);
Auc = auc(:, ~idx);
x = Auc(2,:);
y = Auc(1,:) ;
sort(y)
 plot(x,y);
 %% find best N = 70 
 ratio = 70
 [data_N, feature_name_N,ranks_N,rank_N,feature_indc_N] = SelectFeature_rank( all_data, clinical, ratio, names );

