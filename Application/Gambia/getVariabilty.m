% check for the number of positions where there are variations

folders{1} = 'data/maffilter/chr3L/'; 
folders{2} = 'data/maffilter/chr3R/';
folders{3} = 'data/maffilter/chrX/';

for l = 1 : length(folders)
    loci{l} = dir([folders{l} '*.fa']); 
end
    
for l = 1 : length(folders)
    tmp = strsplit(folders{l}, '/');
    f = fopen([tmp{3} '.tsv'], 'w');
    for i = 1 : length(loci{l})
        if mod(i,100)==0
            disp(i)
        end
        loci{l}(i).name = [folders{l} loci{l}(i).name];
        fasta = fastaread(loci{l}(i).name);

        Seqs = char(fasta.Sequence);
        % get the number of gaps        
        nrgaps = sum(Seqs(:)'=='-') + sum(Seqs(:)'=='N');
        maxgaps = max(sum(Seqs'=='-') + sum(Seqs'=='N'));
        % check the number of variable sites
        diffs=0;
        for j = 1 : length(Seqs)
            if length(unique(Seqs(:,j)))>1
                diffs=diffs+1;
            end
        end
        fprintf(f, '%s\t%d\t%d\t%d\t%d\n',loci{l}(i).name, nrgaps, maxgaps, diffs, length(Seqs));
    end
    fclose(f);
end

