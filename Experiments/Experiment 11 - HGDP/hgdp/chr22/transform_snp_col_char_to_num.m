function chr_num = transform_snp_col_char_to_num(chr)
%transform a SNP column from character array to numerical array
%used to process HGDP data

missing = strcmp(chr,'--');
snps = sort(unique(chr(~missing)));
L = length(chr);
chr_num = zeros(L,1);
for i=1:L
    if strcmp(chr(i),'--')
        chr_num(i) = NaN;
    else
        chr_num(i) = find(strcmp(snps,chr(i))==1)-1;
    end
end