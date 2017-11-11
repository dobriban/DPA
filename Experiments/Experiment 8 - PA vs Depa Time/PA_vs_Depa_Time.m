%Experiment to evaluate accuracy and running time of Derandomized PA vs PA
cd('C:\Dropbox\Projects\Deflation\Experiments\Experiment 8 - PA vs Depa Time')
addpath('C:\Dropbox\Projects\Deflation\Code');
beep off
%% Set parameters
clear, clc
n_mc = 1;
m = 1;
num_selected = zeros(n_mc,1);
runtime = zeros(n_mc,1);
num_selected_depa = zeros(n_mc,1);
runtime_depa = zeros(n_mc,1);
gamma  = 0.6;
theta =  gamma^(1/2)*6;
l = 7;
n_arr =  500*linspace(1,7,7);
n_perm = 20;
%% Effect of signal strength
%Plot mean and +/-sd of number of factors selected as a function of signal
%strength

rng(2);
mean_num_selected =  zeros(l,1);
var_num_selected =  zeros(l,1);
mean_num_selected_depa =  zeros(l,1);
var_num_selected_depa =  zeros(l,1);
mean_runtime =  zeros(l,1);
mean_runtime_depa =  zeros(l,1);

for k=1:l
    n = n_arr(k);
    p = floor(gamma*n);
    for i=1:n_mc
        ti0 = tic;
        Lambda = randn(p,m);
        Lambda = normc(Lambda);
        ep = randn(n,p);
        eta  = randn(n,m);
        eta  = normc(eta);
        X =theta*eta*Lambda'+n^(-1/2)*ep;
        s = svd(X);
        ti1 = toc(ti0);
        
        %PA
        ti = tic;
        s_perm = zeros(n_perm,1);
        for q=1:n_perm
            X_perm= zeros(n,p);
            %get eigenvalues of permutations
            for j=1:p
                pe = randperm(n);
                X_perm(:,j) = X(pe,j);
            end
            svals = svd(X_perm);
            s_perm(q) = svals(1);
        end
        s_perm = sort(s_perm,'descend');
        alpha = 0.05;
        s_thresh = s_perm(max(floor(alpha*n_perm),1));
        num_selected(i) = sum(s>s_thresh);
        runtime(i) = toc(ti)+ti1;
        
        %DPA
        ti = tic;
        deflate = 0;
        num_selected_depa(i) = ddpa(X,deflate);
        runtime_depa(i) = toc(ti)+ti1;
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

%% number of factors selected
rng(2);
savefigs =1; a = {'-','--','-.',':'};
sd_num_selected = var_num_selected.^(1/2);
sd_num_selected_depa = var_num_selected_depa.^(1/2);
figure, hold on
h1 = errorbar(n_arr,mean_num_selected,sd_num_selected,'linewidth',3,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = errorbar(n_arr,mean_num_selected_depa,sd_num_selected_depa,'linewidth',3,'color',rand(1,3));
set(h2,'LineStyle',a{2});
xlabel('n')
ylabel('# Factors Selected')
set(gca,'fontsize',20)
xlim([min(n_arr), max(n_arr)]);
%figure, hold on

legend([h1,h2],{'PA','DPA'},'location','Best')

if savefigs==1
    filename = sprintf( './PA-vs-DEPA-n=%d-p=%d-n-iter=%d.png',n,p,n_mc);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

%% runtime log-log
rng(2);
savefigs =1; a = {'-','--','-.',':'};
figure, hold on
h1 = plot(n_arr,log(mean_runtime),'linewidth',3,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(n_arr,log(mean_runtime_depa),'linewidth',3,'color',rand(1,3));
set(h2,'LineStyle',a{2});
xlabel('n')
ylabel('log(time)')
set(gca,'fontsize',20)
xlim([min(n_arr), max(n_arr)]);
%figure, hold on

legend([h1,h2],{'PA','DPA'},'location','Best')

if savefigs==1
    filename = sprintf( './PA-vs-DEPA-log-time-n=%d-p=%d-n-iter=%d.png',n,p,n_mc);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end
%% evaluate running time
speed_boost = mean_runtime./mean_runtime_depa;
