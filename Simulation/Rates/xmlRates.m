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

bpp_files = dir('BPP/*.txt'); % get the BPP files

system('rm -r IM');
system('mkdir IM');

use_genes = [500];

% loop over all files
for i = 1 : length(bpp_files)
    disp(bpp_files(i).name)
    % filename of the tree
    tree_fname = strrep(bpp_files(i).name, '.txt','.trees');
    
    % open the tree and control file
    t = fopen(['BPP/' tree_fname], 'r');
    c = fopen(['BPP/' bpp_files(i).name], 'r');
    
    % find the species tree in the control file
    while ~feof(c)
        line = fgets(c);
        % # is only used in the tree definition
        if ~isempty(strfind(line,'#'))
            species_tree = line;
        end
    end
    
    % read in all the simulated gene trees
    gene_tree = cell(0,0);
    while ~feof(t)
        tmp = regexprep(fgets(t), '\[TH = (\d*).(\d*)\]','');
        gene_tree{end+1,1} =  strrep(tmp,' ','');clear tmp
    end
    
        % make unique gene names
    for j = 1 : length(gene_tree)
        gene_name = ['gene' num2str(j)];clear tmp
        tmp = phytreeread(gene_tree{j});
        leafname_ori = get(tmp, 'leafnames');
        for k = 1 : length(leafname_ori)
            split_name = strsplit(leafname_ori{k}, '^');
            genes.(gene_name).name{k,1} = [split_name{2} '_' split_name{1}];
            genes.(gene_name).species{k,1} = split_name{2};
            genes.(gene_name).sequence{k,1} = '???';
        end
        % replace the leafnames in the gene tree
        tree = gene_tree{j,1};
        for k = 1 : length(leafname_ori)
            tree = strrep(tree, leafname_ori{k},  genes.(gene_name).name{k,1});
        end
        gTree{j,1} = regexprep(tree,'\n','');
    end
    
    % convert species tree to newick format
    st = regexprep(species_tree, '#(\d*).(\d*)','');
    st = regexprep(st, ' ','');
    tmp = phytreeread(st);
    
    % CHECK IF THIS WORKS FOR TREES WITH MORE THAN 3 LEAVES
    speciesTree = getnewickstr(tmp);
    species = get(tmp, 'leafnames');
    
    for j = 1 : length(use_genes)
        % define target folder for xmls
        folder = [pwd '/IM'];
        filename = strrep(bpp_files(i).name, '.txt', '');
        filename = strrep(filename, 'BPP', ['IM_' num2str(use_genes(j)) 'Genes']);
        
        % get the migration and Ne to make the prior
        tmp = strsplit(bpp_files(i).name, '_');      
        % create xml
        makeAIMxml(folder,filename, 2000000, 2000,...
            genes, species, gTree, speciesTree,...
            false, false, use_genes(j), false);

    end
    
    fclose('all');
end