function [Xnorm, D, U, V] = multi_class_lda_ncomps_full_dataset_py(smallXfilename, fullXfilename, classesfilename, ncomps, outputfilename)

smallX = load(smallXfilename);
fullX = load(fullXfilename);
smallX = smallX';
fullX = fullX';
classes = load(classesfilename);

% Perform the steps to obtain the LDA components using the "smallmatrix" - nonduplicating conditions
X = smallX;
a = sum(abs(X), 1);
b = sum(abs(X), 2);

X = X(b>0, a> 0);
X(X==0) = normrnd(nanmean(X(:)), nanstd(X(:)), length(find(X==0)), 1);

% size(X)
% size(a)
% size(b)
% size(classes)
classes = classes(b>0);

class_labels = unique(classes);
dataMean = nanmean(X, 1);
Sw = zeros(size(X, 2));
Sb = zeros(size(X, 2));

for i = 1:length(class_labels)
    ind = find(classes == class_labels(i));
    if numel(ind) == 0
        continue;
    end
    
    classMean = nanmean(X(ind, :));
    Sw = Sw + cov(X(ind, :), 1);
    Sb = Sb + numel(ind)*(classMean - dataMean)'*(classMean - dataMean);
end

eig_mat = pinv(Sw)*Sb;


[U, D, V] = svd(eig_mat);
a = diag(D)/max(diag(D));
%  stopind = max(find(a>=perc));
stopind = ncomps;

N = V(:, 1:stopind);

% Here is where the full matrix is reintroduced - right before the LDA components are removed!
Xnorm = fullX;
a = sum(abs(Xnorm), 1);
b = sum(abs(Xnorm), 2);

% Note that components will not be removed from "badly behaving" rows/columns,
% aka rows and columns that have no interactions whatsoever
Xnorm(b>0, a>0) = Xnorm(b>0, a>0) - Xnorm(b>0, a>0)*N*N';
Xnorm = Xnorm';
save(outputfilename, 'Xnorm', '-ASCII');
