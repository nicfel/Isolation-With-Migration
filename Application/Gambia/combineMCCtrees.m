% combines mcc tree files into one file (for densitree)
nr_loci = {'s200'};

% specify if outgroupes should be pruned
pruneSpecies={'AchrA1_AchrA1', 'AepiE1_AepiE1'};
pruneSpecies={};

for l = 1 : length(nr_loci)
    % get all species mcc trees (to know which runs are there)
    species_trees = dir(['loci/mcc/' nr_loci{l} '/*species.trees']);
    % build the repo for the combined mcc trees
    system(['rm -r loci/combinedmcc/'  nr_loci{l}]);
    system(['mkdir loci/combinedmcc/'  nr_loci{l}]);


    for s = 1 : length(species_trees)
        % get the name base of the different runs
        name_base = strrep(species_trees(s).name, '_species.trees', '');

        % get the corresponding gene trees    
        loci_trees = dir(['loci/mcc/' nr_loci{l} '/' name_base '*.trees']);

        f = fopen(['loci/combinedmcc/' nr_loci{l} '/' name_base '.trees'],'w');

        % keeps track of the number mapping
        map = cell(0,0);
        prune_leave = false(0,0);
        first = true;
        for i = 1 : length(loci_trees)
            g = fopen(['loci/mcc/' nr_loci{l} '/' loci_trees(i).name]);
            while ~feof(g)
                line = fgets(g);
                if contains(line, 'tree TREE')
                    tmp = strsplit(loci_trees(i).name, '_');
%                     if isempty(pruneSpecies)
%                         fprintf(f, '%s',strrep(line, 'tree TREE1', ['tree ' strrep(tmp{end}, '.trees', '')]));
%                     else
                        tree_string =  strsplit(strtrim(line));
                        % read in the tree
                        ptree = phytreeread(regexprep(tree_string{end}, '\[(.*?)\]', ''));
                        % check which leaves to prune
                        leaf_names = get(ptree, 'leafnames');
                        keep_leave = true(size(leaf_names));
                        for j = 1:length(leaf_names)
                            if prune_leave(str2double(leaf_names{j}))
                                keep_leave(j) = false;
                            end
                        end
                        pruned_tree = subtree(ptree, keep_leave);
                        pruned_tree_string = getnewickstr(pruned_tree);   
                        if contains(tmp{end}, 'species')
                            pruned_tree_string = strrep(pruned_tree_string, ':', '[&location=1]:');
                            pruned_tree_string = strrep(pruned_tree_string, ';', '[&location=2];');
                        else
                            pruned_tree_string = strrep(pruned_tree_string, ':', '[&location=0]:');
                            pruned_tree_string = strrep(pruned_tree_string, ';', '[&location=2];');
                        end
                        fprintf(f, 'tree %s = %s\n', strrep(tmp{end}, '.trees', ''), pruned_tree_string);
%                     end
                    break;            
                else
                    if i == 1       
                        if first
                            if length(pruneSpecies)==2
                                fprintf(f, '#NEXUS\n');
                                fprintf(f, 'Begin taxa;\n');
                                fprintf(f, '	Dimensions ntax=6;\n');
                                fprintf(f, '		Taxlabels\n');
                                fprintf(f, '			1 \n');
                                fprintf(f, '			3 \n');
                                fprintf(f, '			5 \n');
                                fprintf(f, '			6 \n');
                                fprintf(f, '			7 \n');
                                fprintf(f, '			8\n');
                                fprintf(f, '			;\n');
                                fprintf(f, 'End;\n');
                                fprintf(f, 'Begin trees;\n');
                                fprintf(f, '\n');
                                first = false;
                            else
                                fprintf(f, '#NEXUS\n');
                                fprintf(f, 'Begin taxa;\n');
                                fprintf(f, '	Dimensions ntax=8;\n');
                                fprintf(f, '		Taxlabels\n');
                                fprintf(f, '			1 \n');
                                fprintf(f, '			2 \n');
                                fprintf(f, '			3 \n');
                                fprintf(f, '			4 \n');
                                fprintf(f, '			5 \n');
                                fprintf(f, '			6 \n');
                                fprintf(f, '			7 \n');
                                fprintf(f, '			8\n');
                                fprintf(f, '			;\n');
                                fprintf(f, 'End;\n');
                                fprintf(f, 'Begin trees;\n');
                                fprintf(f, '\n');
                                first = false;
                            end
                        end
                        
                        
                        split_line = strsplit(strtrim(line));
                        if length(split_line)==2
                            numval = str2double(split_line{1});
                            if ~isnan(numval)
                                map{numval} = strrep(split_line{2}, ',', '');
                                if sum(ismember(pruneSpecies, map{numval}))==1
                                    prune_leave(numval) = true;
                                else
                                    prune_leave(numval) = false;
                                end
                            end
                        end
                        print=true;
                        for j = 1 : length(pruneSpecies)
                            if contains(line, pruneSpecies{j})
                                print=false;
                            end
                        end
                        if print
                            if contains(line, 'ntax')
                                tmp = regexp(line, '(\d*)', 'match');
%                                 fprintf(f, strrep(line, tmp{1}, num2str(str2double(tmp{1})-length(pruneSpecies))));
                            else
%                                 fprintf(f, line);
                            end
                        end

   
                    end
                end
            end  
            fclose(g);   

        end
        fprintf(f, 'END;');
        fclose(f);
    end
end
