%% HGDP - process data
cd('C:\Dropbox\Projects\Depa\Experiments\Experiment 11 - HGDP\hgdp\chr22')
%% process cell arrays into numeric arrays
%in this section I assume we have processed and saved the HGDP data
%into matlab cell arrays
load('chr22_geno.mat','chr22')
chr22 = chr22'; %transpose to get nxp
[n,p] = size(chr22);
chr22_num = zeros(n,p);
print_iter=1;
tic
for i=1:p
    %pr
    if (mod(i,floor(p/50))==0)
        if print_iter==1
            str = sprintf('Exp %d out of %d.\n',i,p);
            fprintf(str);
            toc
        end
    end
    chr22_num(:,i) = transform_snp_col_char_to_num(chr22(:,i));
end
%%
save('chr22_geno_num.mat','chr22_num')
%% 
%in this section I assume we have processed and saved the HGDP data
%into matlab numeric arrays
clear
load('chr22_geno_num.mat','chr22_num')
[n,p] = size(chr22_num);
%% plot allele frequencies
allele_freq = nanmean(chr22_num)/2;
MAF = min(allele_freq,1-allele_freq);
hist(MAF,floor(sqrt(p)));
xlabel('MAF'); ylabel('freq'); sf('MAF_HGDP_chr22');

%% de-mean
%take the mean of the non-missing data
%center the non-missing data
num_non_missing = sum(~isnan(chr22_num));
means = nanmean(chr22_num);
c22_centered = chr22_num - ones(n,1)*means;
m = nanmean(c22_centered);
%%
save('chr22_geno_centered.mat','c22_centered')
clear
%% standardize
%standardize non-missing data
nanvars = nanvar(c22_centered);
c22_norm = c22_centered./(ones(n,1)*sqrt(nanvars));
no = nanvar(c22_norm);
%%
save('chr22_geno_norm.mat','c22_norm')
clear


