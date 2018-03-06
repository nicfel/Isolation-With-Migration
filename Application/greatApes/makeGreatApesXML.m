%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make IM xml to estimate the population size and migration rates for a
% fixed gene and species tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% add software folder to path
addpath('../../Software')

% get all the sequences
sequences = dir('Sequences/*.fasta'); 

system('rm -r IM');
system('mkdir IM');

species = cell(0,0);
% convert the fasta data
for j = 1 : length(sequences)
    tmp = strsplit(sequences(j).name,'.');
    gene_name = tmp{1};
    fasta = fastaread(['Sequences/' sequences(j).name]);
    for k = 1 : length(fasta)
        name = strrep(fasta(k).Header,'-','_');
        genes.(gene_name).name{k,1} = name;
        tmp = strsplit(name,'_');
        genes.(gene_name).species{k,1} = tmp{1};
        species{end+1,1} = tmp{1};
        genes.(gene_name).sequence{k,1} = fasta(k).Sequence;        
    end
end       

species = unique(species);
% make xml files for all possible combination of species
for i = 4 : 2 : length(species)
    if i == 4
        system('mkdir IM/chimps');
        combinations = {'Ppa'    'Pts'    'Ptt'    'Ptv'};
        folder = [pwd '/IM/chimps/'];
    else
        system('mkdir IM/greatapes');
        folder = [pwd '/IM/greatapes/'];
        combinations = combnk(species,i);
    end
    
    for j = 1 : size(combinations,1)
        clear use_species
        use_species = combinations(j,:);
        for rep = 0 : 9
            name = 'Chimp';
            for k = 1 : length(use_species)
                name = [name '_' use_species{k}];
            end
            filename = [name '_rep' num2str(rep)];

            makeAIMxml(folder,filename, 200000000, 1000000,...
                genes, use_species, [], [],...
                true,true,length(sequences), true);
        end
    end
end


    
