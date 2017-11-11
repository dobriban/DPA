function epsi = epsi_thresh(n,p)

sigma = (n^(1/2)+p^(1/2))*(n^(-1/2)+p^(-1/2))^(1/3);
mu = (n^(1/2)+p^(1/2))^2;

%test = 'old';
test = 'new';

switch test
    case 'old'
        tau = (sigma/n)^(1/2);
    case 'new'
        tau = sigma/mu;
end

epsi = 3*tau;
