%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script genereates Sequences for the simulated trees using BPP. The
% sequences are generated using an HKY model and seq-gen v1.3.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% add software folder to path
addpath('../../../Software')

% delete Sequences folder and create it newly
system('rm -r Sequences');
system('mkdir Sequences/')

bpp_files = dir('BPP/*.trees'); % get the trees

% print the relative rates to a file
f_rates = fopen('Sequences/relative_rates.txt', 'w');
    fprintf(f_rates, '%s\n', strtrim(sprintf('gene%d\t',[1:50])));


for j = 1 : length(bpp_files)
    % read in the tree file
    
    %sample the r=elative rates of each loci
    relative_rates = exprnd(1,50,1);
    relative_rates = relative_rates./sum(relative_rates)*50;
    relative_rates = relative_rates;
    
    tmp = strsplit(bpp_files(j).name, '_');
    tmp2 = strsplit(tmp{end}, '.');
    
    fprintf(f_rates, '%s\t%s\n', tmp2{1}, strtrim(sprintf('%.12f\t',relative_rates)));
    
    f = fopen(sprintf('BPP/%s',bpp_files(j).name),'r');
    geneNumber=1;
    c = 1;
    tree_nr = 1;
    while ~feof(f)
        line = fgets(f);
        g = fopen('tree.tree','w');
        fprintf(g, line);fclose(g);
        command = sprintf('%s -mHKY -t3.0 %s -l 1000 -s %.12f < tree.tree > %s',...
            '../../Software/seq-gen1.3.3',...
                '-f0.3,0.2,0.2,0.3',relative_rates(geneNumber),'seq.fasta');
        system(command);
        tree_nr = tree_nr + 1;

        % read the sequence file
        s = fopen('seq.fasta', 'r');
        % read file
        seq = textscan(s, '%s'); fclose(s);
        if length(seq{1}) > 0
            for k = 3 : 2 : length(seq{1})
                tmp = strsplit(seq{1}{k},'^');
                Data(c).Header = sprintf('gene%d.%s_%s',geneNumber,tmp{2},tmp{1});
                Data(c).Sequence = seq{1}{k+1};
                c = c + 1;
            end       
        else
            s = fopen('no_variation.fasta', 'r');
            % read file
            seq = textscan(s, '%s'); fclose(s);

            for k = 3 : 2 : length(seq{1})
                tmp = strsplit(seq{1}{k},'^');
                Data(c).Header = sprintf('gene%d.%s',geneNumber,tmp{1});
                Data(c).Sequence = seq{1}{k+1};
                c = c + 1;
            end       
        end
        geneNumber = geneNumber + 1;
    end   
    fclose(f);

    % write the fasta file to the target folder
    fastawrite(sprintf('Sequences/%s',...
        strrep(bpp_files(j).name,'trees','fasta')),Data);
    system('rm tree.tree');system('rm seq.fasta');
end
% fclose(f_rates);
