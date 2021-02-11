% combines mcc tree files into one file (for densitree)
nr_loci = {'s200'};

for l = 1 : length(nr_loci)
    % get all species mcc trees (to know which runs are there)
    species_trees = dir(['loci/mcc/' nr_loci{l} '/*species.trees']);
    % build the repo for the combined mcc trees
    system(['rm -r loci/ancestors/'  nr_loci{l}]);
    system(['mkdir loci/ancestors/'  nr_loci{l}]);


    for s = 1 : length(species_trees)
        % get the name base of the different runs
        name_base = strrep(species_trees(s).name, '_species.trees', '');

        % get the corresponding gene trees    
        loci_trees = dir(['loci/mcc/' nr_loci{l} '/' name_base '*.trees']);

        f = fopen(['loci/ancestors/' nr_loci{l} '/' name_base '.csv'],'w');
        
        fprintf(f, 'species1,species2,distance,tree\n');
        
        numval = cell(0,0);
        species = cell(0,0);
        
        % keeps track of the number mapping
        map = cell(0,0);
        prune_leave = false(0,0);
        for i = 1 : length(loci_trees)
            g = fopen(['loci/mcc/' nr_loci{l} '/' loci_trees(i).name]);
            while ~feof(g)
                line = fgets(g);
                if contains(line, 'tree TREE')
                    tmp = strsplit(loci_trees(i).name, '_');
                    tree_string =  strsplit(strtrim(line));
                    % read in the tree
                    ptree = phytreeread(regexprep(tree_string{end}, '\[(.*?)\]', ''));
                    % check which leaves to prune
                    leaf_names = get(ptree, 'leafnames');
                    
                    for a = 1 : length(numval)-1
                        for b = a+1:length(numval)
                            % get the subtree
                            inda = find(ismember(leaf_names, numval{a}));
                            indb = find(ismember(leaf_names, numval{b}));
                            keep_leave = true(size(leaf_names));
                            keep_leave([inda, indb])  = false;
                            tree = prune(ptree, keep_leave);
                            distances = pdist(tree);
                            if contains(tmp{end}, 'species')
                                fprintf(f, '%s,%s,%f,species\n', species{a}, species{b}, distances);
                            else
                                fprintf(f, '%s,%s,%f,gene\n', species{a}, species{b}, distances);
                            end
                        end
                    end                    
                    break;
                else
                    if i == 1
                        split_line = strsplit(strtrim(line));
                        if length(split_line)==2
                            val = str2double(split_line{1});
                            if ~isnan(val)
                                numval{end+1,1} = split_line{1};
                                species{end+1,1} = strrep(split_line{2}, ',','');
                            end
                        end   
                    end
                end
            end  
            fclose(g);   

        end
        fclose(f);
    end
end
