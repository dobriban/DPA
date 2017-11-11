%Experiment to evaluate accuracy of DDPA with increased
%threshold
cd('C:\Dropbox\Projects\Deflation\Experiments\Experiment 10 - Depa Increase Thesh')
addpath('C:\Dropbox\Projects\Deflation\Code');
beep off
%% Expt 1: one-factor model
%% Set parameters
n_mc = 1e1;
rng(2);
n = 500;
gamma  = 0.6;
p = floor(gamma*n);
m = 1;

%% Effect of signal strength
%Plot mean and +/-sd of number of factors selected as a function of signal
%strength

rng(2);
l = 10;
%theta_arr = linspace(0.1,1,l);
sig=linspace(0.2,6,l);
theta_arr =  gamma^(1/2)*sig;
mean_num_selected_ddpa =  zeros(l,1);
var_num_selected_ddpa =  zeros(l,1);
mean_runtime =  zeros(l,1);
mean_runtime_depa =  zeros(l,1);

for k=1:l
    theta = theta_arr(k); %factor strength
    for i=1:n_mc
        Lambda = randn(p,m);
        Lambda = normc(Lambda);
        ep = randn(n,p);
        eta  = randn(n,m);
        eta  = normc(eta);
        X =theta*eta*Lambda'+n^(-1/2)*ep;
        s = svd(X);
        
        %DDPA
        deflate = 1;
        epsi = 0;
        num_selected_ddpa(i) =ddpa(X,deflate,epsi);
        
        %DDPA+ (PG threshold)
        deflate = 1;
        num_selected_ddpa_plus(i) = ddpa(X,deflate);
    end
    mean_num_selected_ddpa(k) = mean(num_selected_ddpa);
    var_num_selected_ddpa(k) = var(num_selected_ddpa);
    mean_num_selected_ddpa_plus(k) = mean(num_selected_ddpa_plus);
    var_num_selected_ddpa_plus(k) = var(num_selected_ddpa_plus);
    
end

%%
rng(2);
savefigs =1; a = {'-','--','-.',':'};
sd_num_selected_ddpa = var_num_selected_ddpa.^(1/2);
sd_num_selected_ddpa_plus = var_num_selected_ddpa_plus.^(1/2);
figure, hold on
h1 = errorbar(sig,mean_num_selected_ddpa,sd_num_selected_ddpa,'linewidth',3,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = errorbar(sig,mean_num_selected_ddpa_plus,sd_num_selected_ddpa_plus,'linewidth',3,'color',rand(1,3));
set(h2,'LineStyle',a{2});
xlabel('Signal Strength')
ylabel('# Factors Selected')
set(gca,'fontsize',20)
xlim([min(sig), max(sig)]);
%figure, hold on

legend([h1,h2],{'DDPA','DDPA+'},'location','Best')

if savefigs==1
    filename = sprintf( './DEPA-Increase-Thesh-n=%d-p=%d-n-iter=%d.png',n,p,n_mc);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

%% Expt 2: three-factor model, one strong
%% Set parameters
n_mc = 1e1;
rng(2);
n = 500;
p = 300;
m = 3;
num_selected = zeros(n_mc,1);
num_selected_inc = zeros(n_mc,1);

%% Effect of signal strength
%Plot mean and +/-sd of number of factors selected as a function of signal
%strength

rng(2);
l = 20;
gamma  = p/n;
sig1 = 6*ones(l,1);
sig2 = 10*ones(l,1);
sig3 = linspace(10,70,l);
theta_arr1 =  gamma^(1/2)*sig1;
theta_arr2 =  gamma^(1/2)*sig2;
theta_arr3 =  gamma^(1/2)*sig3;
mean_num_selected =  zeros(l,1);
var_num_selected =  zeros(l,1);
mean_num_selected_inc =  zeros(l,1);
var_num_selected_inc =  zeros(l,1);

theta = zeros(m);
for k=1:l
    theta(1,1) = theta_arr1(k); %factor strength
    theta(2,2) = theta_arr2(k);
    theta(3,3) = theta_arr3(k);
    for i=1:n_mc
        Lambda = randn(p,m);
        Lambda = normc(Lambda);
        ep = randn(n,p);
        eta  = randn(n,m);
        eta  = normc(eta);
        X =eta*theta*Lambda'+n^(-1/2)*ep;
        
        %DDPA
        deflate = 1;
        epsi = 0;
        num_selected(i) =ddpa(X,deflate,epsi);
        
        %DDPA+
        deflate = 1;
        num_selected_inc(i) = ddpa(X,deflate);
    end
    mean_num_selected(k) = mean(num_selected);
    var_num_selected(k) = var(num_selected);
    mean_num_selected_inc(k) = mean(num_selected_inc);
    var_num_selected_inc(k) = var(num_selected_inc);
end


%%
rng(2);
savefigs =1; a = {'-','--','-.',':'};
sd_num_selected = var_num_selected.^(1/2);
sd_num_selected_inc = var_num_selected_inc.^(1/2);
figure, hold on
h1 = errorbar(sig3,mean_num_selected,sd_num_selected,'linewidth',3,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = errorbar(sig3,mean_num_selected_inc,sd_num_selected_inc,'linewidth',3,'color',rand(1,3));
set(h2,'LineStyle',a{2});
xlabel('Signal strength')
ylabel('# Factors Selected')
set(gca,'fontsize',20)
xlim([min(sig3), max(sig3)]);
%figure, hold on

legend([h1,h2],{'DDPA','DDPA+'},'location','Best')

if savefigs==1
    filename = sprintf( './DEPA-Increase-Thesh-n=%d-p=%d-n-iter=%d-n-fac=%d.png',n,p,n_mc,m);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end