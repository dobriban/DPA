function [above_thresh] = perry_gavish_donoho(s, gamma,t)
%checks if top spike of s is above the perry_gavish_donoho thresold

addpath('C:\Git\EigenEdge\Code');

lambda = s(1)^2;
evals = s(2:length(s)).^2;

[ell,cos_right,cos_left] = general_spiked_estimate(lambda,evals,gamma);

%The formula is that
%sigma_1(X)<2*theta*c_1*c_2
%and ell  = theta^2 [?]
%this seems correct, double-checked in my paper with Will and Amit
test_stat = lambda - 4*ell*cos_left*cos_right;

if test_stat<0
    above_thresh = 1;
else
    above_thresh = 0;
end


%thresh = (4*ell*cos_left*cos_right)^(1/2);


%self-consistency: not bad
if exist('t','var')
    w = ones(length(t),1)/length(t);
    [lambda_forw,cos_right_forw,cos_left_forw] = general_spiked_forward(ell, t, w, gamma);
end
