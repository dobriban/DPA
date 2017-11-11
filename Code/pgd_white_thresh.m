function [thresh_lambda,thresh_ell] = pgd_white_thresh(gamma)
%compute the perry_gavish_donoho threshold for white noise

%from eq. 11. of Gavish, Donoho - The Optimal Hard Threshold for Singular Values is  4 over root 3 - 2014 - IEEE IT
addpath('C:\Git\EigenEdge\Code');

a = gamma;
b = a+1+sqrt(a^2+14*a+1);
thresh = 2*(a+1)+8*a/b;
thresh_lambda = thresh^(1/2);
thresh_ell= standard_spiked_inverse(thresh_lambda^2,gamma);
thresh_ell = thresh_ell^(1/2);