function [num_sel, runtime,thresh] = ddpa(X,deflate,epsi)
%Deterministic Parallel Analysis with deflation

%Inputs
%X - data matrix, each row is an independent sample
%deflate - (optional) binary variable, indicating if deflation should be performed
%                default: deflate=1, do deflation
%epsi - (optional) positive real variable, indicating the increase in threshold at each step
%               of DDPA.

%Examples:
% ddpa(X): DDPA+
% ddpa(X,0): DPA
% ddpa(X,0,epsi): DPA with increased threshold
% ddpa(X,1): DDPA+
% ddpa(X,1,0): DDPA

if ~exist('deflate','var')
    deflate = 1;
end

ti = tic;

[n,p] = size(X);
gamma = p/n;
it = 1;
num_sel = 0;



%DPA: no deflate
if deflate==0
    s = svd(X);
    t = sum(X.^2, 1)';
    u = upper_edge(t, gamma);
    num_sel = sum(s.^2>u);
    thresh = u^(1/2);
    %DDPA: deflate
else
    [U,s,V] = svd(X);
    s = diag(s);
    while (it==1)&&(~isempty(s))
        %Option 1 (Default): if epsi is not provided, then select above PGD threshold
        if ~exist('epsi','var')
            [above_threshold] = perry_gavish_donoho(s, gamma);
            
            %Option 2: if epsi is provided, then select above u*(1+epsi)
        else
            t = sum(X.^2, 1)';
            u = upper_edge(t, gamma);
            thresh = u^(1/2)*(1+epsi);
            above_threshold = (s(1)>thresh);
            
        end
        
        if(above_threshold==1)
            num_sel = num_sel+1;
            X = X - U(:,1)*s(1)*V(:,1)'; %deflate
            s = s(2:length(s));
            U = U(:,2:size(U,2));
            V = V(:,2:size(V,2));
        else %stop iteration
            it=0;
            
            %determine threshold
            if ~exist('epsi','var')
                thresh = s(1); %below the thresold
                grid = linspace(thresh,2*thresh,50);
                for i=1:50
                    ab_thr = perry_gavish_donoho([grid(i); s(2:length(s))], gamma);
                    if ab_thr==1
                        thresh = grid(i);
                        break
                    end
                end
            end
        end
    end
end

runtime = toc(ti);
if thresh==Inf
    thresh=0;
end

