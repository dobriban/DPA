%Experiment to evaluate accuracy of Derandomized PA vs PA
cd('C:\Dropbox\Projects\Deflation\Experiments\Experiment 5 - PA vs Depa')
addpath('C:\Dropbox\Projects\Deflation\Code');
beep off
%% Set parameters
n_mc = 1e1;
rng(2);
n = 500;
gamma  = 0.6;
p = floor(gamma*n);
m = 1;
num_selected = zeros(n_mc,1);

%% Effect of signal strength
%Plot mean and +/-sd of number of factors selected as a function of signal
%strength

rng(2);
l = 10;
%theta_arr = linspace(0.1,1,l);
sig=linspace(0.2,6,l);
theta_arr =  gamma^(1/2)*sig;
mean_num_selected =  zeros(l,1);
var_num_selected =  zeros(l,1);
mean_num_selected_depa =  zeros(l,1);
var_num_selected_depa =  zeros(l,1);
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
        
        %PA
        [num_selected(i),runtime(i)] = PA_BE(X);
        
        %DPA
        ti = tic;
        deflate = 0;
        num_selected_depa(i) = ddpa(X,deflate);
        runtime_depa(i) = toc(ti);
        
    end
    mean_num_selected(k) = mean(num_selected);
    var_num_selected(k) = var(num_selected);
    mean_num_selected_depa(k) = mean(num_selected_depa);
    var_num_selected_depa(k) = var(num_selected_depa);
    
    mean_runtime(k) =  mean(runtime);
    var_runtime(k) =  var(runtime);
    mean_runtime_depa(k) =  mean(runtime_depa);
    var_runtime_depa(k) =  var(runtime_depa);
    
end

%%
rng(2);
savefigs =1; a = {'-','--','-.',':'};
sd_num_selected = var_num_selected.^(1/2);
sd_num_selected_depa = var_num_selected_depa.^(1/2);
figure, hold on
h1 = errorbar(sig,mean_num_selected,sd_num_selected,'linewidth',3,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = errorbar(sig,mean_num_selected_depa,sd_num_selected_depa,'linewidth',3,'color',rand(1,3));
set(h2,'LineStyle',a{2});
xlabel('Signal Strength')
ylabel('# Factors Selected')
set(gca,'fontsize',20)
xlim([min(sig), max(sig)]);
%figure, hold on

legend([h1,h2],{'PA','DPA'},'location','Best')

if savefigs==1
    filename = sprintf( './PA-vs-DEPA-n=%d-p=%d-n-iter=%d.png',n,p,n_mc);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

