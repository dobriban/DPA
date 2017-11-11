function [num_sel,runtime,threshold] = PA_BE(X,p_thresh,deflate)
% Estimate number of factors by Buja& Eyuboglu's Parallel Analysis

use_p_thresh = 0;
if ~exist('p_thresh','var')
    p_thresh = 95;
end
ti = tic;
[n,p] = size(X);

if ~exist('deflate','var')
    deflate = 0;
    epsi = 0;
end

if deflate==1
    if ~exist('epsi','var')
        epsi = epsi_thresh(n,p);
    end
end


s = svd(X);
num_perm=19;%20;%19;
X_perm= zeros(n,p);
%get eigenvalues of permutations
s_perm = zeros(min(n,p),num_perm);
rng(2);
for i=1:num_perm
    toc(ti);
    for j=1:p
        pe = randperm(n);
        X_perm(:,j) = X(pe,j);
    end
    s_perm(:,i) = svd(X_perm);
end

factor_retained=1;
num_factor=0;

while (factor_retained==1)&&(num_factor<min(n,p))
    num_factor= num_factor+1;
    %threshold is empirical percentile of permutation distribution
    if use_p_thresh==1
        thresh = prctile(s_perm(num_factor,:),p_thresh);
    else
        thresh = max(s_perm(num_factor,:)); 
    end
    %if data singular value is less than the threshold, don't keep it
    test = (s(num_factor)^2<=thresh^2*(1+epsi));
    threshold = thresh*(1+epsi);
    
    if(test==1)
        factor_retained=0;
        num_factor= num_factor-1; %decrease by 1 to go back to the pre-error state
    end
end
num_sel  = num_factor;
runtime = toc(ti);

