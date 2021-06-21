function [data, feature_name,ranks,rank,feature_indc] = SelectFeature_rank( data, clinical, ratio, names )
% Feature selection

% this function selects feature from the matrix
% all_data is a value matrix
% matched_clinical is the original matched_clinical dataset
% n is the feature number to be selected
%
% it's output is a matrix with selected features



%% discretize by using mean+/-alpha*std. Here let alpha=0.5 for simplify
m = mean(data);
s = std(data);
disc_all_data = zeros(size(data));
for i = 1: size(data,2)
  disc_all_data(data(:, i)>=m(i)+s(i)/2, i) = 1;
  disc_all_data(data(:, i)<=m(i)-s(i)/2, i) = -1;
end

%% mrmr
if ratio < 1
    seln = floor(size(all_data,2)*n);
else
    seln = ratio;
end
feature_indc =  mrmr_miq_d(disc_all_data, clinical, seln);

%% rank

random_num= sort(feature_indc);
a = [];b = [];c = [];d = [];e = [];
k1=1; k2=1; k3=1; k4=1; k5=1;
for i = 1:ratio
    com = random_num(i);
    if com <= 216
        a(k1) = i;  % µ±Ç°µÄÁÐÊý
        k1 = k1+1;
        
    elseif com <= 811
        b(k2) = i;
        k2 = k2+1;
    elseif com <= 16562
        c(k3) = i;
        k3 = k3+1;
    elseif com <= 30182
        d(k4) = i;
        k4 = k4+1;
    elseif com <= 56170
        e(k5) = i;
        k5 = k5+1;
    end
    
end
if isempty(a)
    a = [0 0];
end
if isempty(b)
    b = [0 0];
end
if isempty(c)
    c = [0 0];
end
if isempty(d)
    d = [0 0];
end
if isempty(e)
    e = [0 0];
end

ranks = { [a(1):a(end)], [b(1):b(end)], [c(1):c(end)], [d(1):d(end)],  [e(1):e(end)]};
rank = [a(1) a(end); b(1) b(end); c(1) c(end); d(1) d(end); e(1) e(end)  ];

%%
%data = all_data(:, feature_indc);
data = all_data(:, random_num);
feature_name = names(random_num);



end

