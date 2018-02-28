%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script first simulates trees with a given number of leaves and a
% given effective population size. These tree then constitute the species
% tree from which a gene tree is simulated using BPP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% delete BPP folder and create it newly
system('rm -r BPP');
system('mkdir BPP');
pause(1)


migration_rate = 0.05; % define the migration rates


geneNe = 1; % define the Ne for the gene tree



distance = [0 5 10 20 0 0 0;
            0 0 10 20 0 0 0;
            0 0 0 20 0 0 5;
            0 0 0 0 0 10 15;
            0 0 0 0 0 0 0;
            0 0 0 0 0 0 0;
            0 0 0 0 0 0 0;];
        
distance = distance*1;

for i = 1 : size(distance)
    for j = (i+1) : size(distance)
        abs_mig(i,j) = migration_rate/(distance(i,j));
        abs_mig(j,i) = migration_rate/(distance(i,j));
    end
end
% set the random number generator
rng(1);

% sample the coalescent times and the individuals 
for i = 1:100
    % build the input file for BPP to get the expected node order
    f = fopen(sprintf('BPP/BPP_S%d.txt',i),'w');
    relativeNe =  lognrnd(-0.1250,0.2,1,7)*8;
    fprintf(f, '# Ne = %s\n', sprintf('%f ', relativeNe));
    
    for a = 1 : size(distance,1)
        for b = 1 : size(distance,2)
            mig_rate(a,b) = abs_mig(a,b);
            if isnan(mig_rate(a,b))
                rel_rate(b) = 0;
                mig_rate(a,b) = 0;
            else
               rel_rate(b) = exprnd(1);
               mig_rate(a,b) = mig_rate(a,b) * rel_rate(b) ;
            end
        end
        fprintf(f, '# %s\n', sprintf('%f ', rel_rate));
    end   
    
    
    fprintf(f, 'seed = %d\n\n', randi(100000));
    fprintf(f, 'treefile = BPP_S%d.trees\n',i);
    fprintf(f, 'species&tree = 4 A B C D\n');
    fprintf(f, '\t\t 20 10 5 15\n\n');
    fprintf(f, '(((A #%f,B #%f):5 #%f,C #%f):10 #%f, D #%f):20 #%f;\n',...
        relativeNe(1), relativeNe(2), relativeNe(7), relativeNe(3), relativeNe(6), relativeNe(4), relativeNe(5));
    fprintf(f, 'migration = 7\n');
    fprintf(f, '\n  ');    
    nodes = {'A', 'B', 'C', 'D', 'ABCD', 'ABC', 'AB'};
    for s = 1 : length(nodes)
        fprintf(f, '%s  ',nodes{s});
    end
    for s1 = 1 : length(nodes)
        fprintf(f, '\n%s  ',nodes{s1});
        for s2 = 1 : length(nodes)
            fprintf(f, '%.10f  ', mig_rate(s1,s2)*(relativeNe(s2)/4));
        end
    end
    fprintf(f, '\nloci&length = 500 1\n');
    fclose(f);
end
fclose('all');
delete('speciesTree_old.ctl');system('rm *.trees');

cd('BPP/')
% run all files in the BPP folder
ctl_files = dir('*.txt');
for i = 1 : length(ctl_files)
    system(sprintf('../../../Software/MCcoal %s',ctl_files(i).name))
end
cd('..');
