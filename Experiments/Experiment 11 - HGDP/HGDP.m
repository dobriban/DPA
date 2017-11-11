%% HGDP PA
cd('C:\Dropbox\Projects\Deflation\Experiments\Experiment 11 - HGDP')
addpath('C:\Dropbox\Projects\Deflation\Experiments\Experiment 11 - HGDP\hgdp\chr22')
addpath('C:\Dropbox\Projects\Deflation\Code')
%% load the raw data processed into Matlab format
clear; 
normalized = 1;
switch normalized
    case 1
        %normalized
        load('chr22_geno_norm.mat');
        X = c22_norm; clear('c22_norm');
    case 0
        %centered
        load('chr22_geno_centered.mat');
        X = c22_centered; clear('c22_centered');
end

%set to 0 missing entries;
X(isnan(X))=0;
[n,p] = size(X);
X = n^(-1/2)*X;
%% All 4 methods
[num_sel_pa,~,thresh_pa] = PA_BE(X); %PA
[num_sel_dpa,~,thresh_dpa] = ddpa(X,0); %DPA
[num_sel_ddpa,~,thresh_ddpa] = ddpa(X,1,0); %DDPA
[num_sel_ddpa_plus,~,thresh_ddpa_plus] = ddpa(X); %DDPA+
%
s = svd(X);
%%
rng(2);
savefigs =1; a = {'-','--','-.',':'};
figure, hold on
hist(s,2*floor(sqrt(p)));
xlim([min(s(s>0)), max(s)]);
%pa=s(num_sel_pa); %your point goes here
pa=thresh_pa; %your point goes here
x=[pa,pa];
y=get(gca,'Ylim');
h1 = plot(x,y,'linewidth',2,'color',rand(1,3));
set(h1,'LineStyle',a{1});
%pa=s(num_sel_dpa); %your point goes here
pa=thresh_dpa; %your point goes here
x=[pa,pa];
y=get(gca,'Ylim');
h2 = plot(x,y,'linewidth',2,'color',rand(1,3));
set(h2,'LineStyle',a{2});
%pa=s(num_sel_ddpa); %your point goes here
pa=thresh_ddpa; %your point goes here
x=[pa,pa];
y=get(gca,'Ylim');
h3 = plot(x,y,'linewidth',2,'color',rand(1,3));
set(h3,'LineStyle',a{3});
%pa=s(num_sel_ddpa_plus); %your point goes here
pa=thresh_ddpa_plus; %your point goes here
x=[pa,pa];
y=get(gca,'Ylim');
h4 = plot(x,y,'linewidth',2,'color',rand(1,3));
set(h4,'LineStyle',a{4});
legend([h1,h2,h3,h4],{'PA','DPA','DDPA','DDPA+'},'location','Best')
xlabel('Singular value')
if savefigs==1
    filename = sprintf( './hgdp-hist-all.png');
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end



%% scatterplots
tic;
[U,s,V] = svd(X);
toc;
s =diag(s);
%% these show the clustering of the samples
k = 12;
savefigs =1;
for i=1:k
    figure, scatter(U(:,2*i-1),U(:,2*i),'filled');
    xlabel(2*i-1);
    ylabel(2*i);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    if savefigs==1
        filename = sprintf( './hgdp-scatter-i=%d.png',i);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        close(gcf)
    end
    set(gca,'fontsize',20)
end

%% self-consistency: sample splitting
X = n^(1/2)*X; %un-normalize
%%
rng(2);
num_subsample = 50;
tic;
for i=1:num_subsample
    toc;
    pe = randperm(n);
    ind = zeros(n,1); ind(1:floor(n/2)) = ones(floor(n/2),1);
    Xs=X(pe(ind==1),:);
    num_sel_de_sub(i) = ddpa(Xs,0);
    num_sel_def_sub(i) = ddpa(Xs);
end
m_de_half =  mean(num_sel_de_sub);
m_def_half =  mean(num_sel_def_sub);
sd_de_three_quarter =  var(num_sel_de_sub)^(1/2);
sd_def_three_quarter =  var(num_sel_def_sub)^(1/2);

%%
alpha = 3/4;
rng(2);
num_subsample = 50;
tic;
for i=1:num_subsample
    toc;
    pe = randperm(n);
    ind = zeros(n,1); ind(1:floor(alpha*n)) = ones(floor(alpha*n),1);
    Xs=X(pe(ind==1),:);
    num_sel_de_sub(i) = ddpa(Xs,0);
    num_sel_def_sub(i) = ddpa(Xs);
end
m_de_three_quarter =  mean(num_sel_de_sub);
m_def_three_quarter =  mean(num_sel_def_sub);
sd_de_three_quarter =  var(num_sel_de_sub)^(1/2);
sd_def_three_quarter =  var(num_sel_def_sub)^(1/2);


%% plot two histograms
rng(2);
X_perm= zeros(n,p);
for j=1:p
    pe = randperm(n);
    X_perm(:,j) = X(pe,j);
end
s_perm = svd(X_perm);
s = svd(X);
savefigs =1;
figure, hold on
histogram(s_perm,floor(sqrt(p)));
histogram(s,floor(sqrt(p)),'FaceColor','red');
legend({'X_\pi','X'},'location','Best')
xlabel('Singular value')

if savefigs==1
    filename = sprintf( './hgdp-perm.png');
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end
