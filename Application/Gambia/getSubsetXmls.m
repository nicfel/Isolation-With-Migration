function [] = getSubsetXmls(nr_samples, nrSubset,...
    nrRep, temp, chainLength, autodelay, chromosomeWeights,...
    maxNr,chrName, subs, mig_prior, fixed)
% get at random subset of the anophelese sequences including the x
% chromosome


loci_names = {'chr2R', 'chr3L', 'chr3R', 'chrX'};


for subset = 1 : length(nr_samples)
    for r = 1 : subs
        clear loci

        system(sprintf('mkdir xmls/s%d',nr_samples(subset)));


        % read in the loci names
        for i = 2 : length(loci_names)
            f = fopen([loci_names{i} '.tsv']); t = textscan(f, '%s\t%d\t%d\t%d\t%d'); fclose(f);
            gaps = double(t{:,3})./double(t{:,5});
            diff = double(t{:,4})./double(t{:,5});
            loci{i} = t{1}(find(diff>0 & diff<0.5),1);hold on
            maxNr(i) = length(loci{i});
            disp(length(loci{i}))
        end
            
%         % make the names complete with directory
%         for i = 1 : length(loci{1})
%             loci{1}(i).name = ['data/maffilter/chr2R/' loci{1}(i).name];
%         end
%         for i = 1 : length(loci{2})
%             loci{2}(i).name = ['data/maffilter/chr3L/' loci{2}(i).name];
%         end
%         for i = 1 : length(loci{3})
%             loci{3}(i).name = ['data/maffilter/chr3R/' loci{3}(i).name];
%         end
%         for i = 1 : length(loci{4})
%             loci{4}(i).name = ['data/maffilter/chrX/' loci{4}(i).name];
%         end
      

        % sample how many loci of which chromosome to take
        chromosome = randsample(4, nr_samples(subset), true, chromosomeWeights);
        
        % for each chromsome, sample random loci
        for i = 1 : length(loci)  
            indices{i} = randsample(maxNr(i),sum(chromosome==i));
        end
        system('rm -r sequencestmp');
        system('mkdir sequencestmp');

        for i = 1 : length(indices)
            for j = 1 : length(indices{i})
                fasta = fastaread(loci{i}{indices{i}(j)});
                c = 1;
                for k = 1 : length(fasta)
                    if ~contains(fasta(k).Header, 'AgamP4')
                        
%                     && ...
%                         isempty(strfind(fasta(k).Header, 'AchrA1'))
                        newFasta(c) = fasta(k);
                        newFasta(c).Header = [fasta(k).Header '_' fasta(k).Header];
                        c = c+1;
                    end
                end
                nogap = true(1,length(newFasta(1).Sequence));
                for k = 1 : length(newFasta(1).Sequence)
                    for l = 1 : length(newFasta)
                        if strcmp(newFasta(l).Sequence(k),{'-'})
                            nogap(k) = false;
                            break;
                        end 
                    end
                end
                
                for l = 1 : length(newFasta)
                    newFasta(l).Sequence = newFasta(l).Sequence(nogap);
                end

                
                
                
                new_name = strsplit(loci{i}{indices{i}(j)}, '/');
                fastawrite(['sequencestmp/' new_name{end}], newFasta);
            end
        end

        % read in all files from the sequences folder.
        loci = dir('sequencestmp/*.fa');


        loci_length = zeros(length(loci),1);


        for mig_prior_it = 1 : length(mig_prior)
            name_base = sprintf('anopheles_%s_s%d_sub%d_mig%d', chrName, nr_samples(subset),r,mig_prior_it);
    
            if fixed
                f = fopen('template_fixed.xml');
            else
                f = fopen('template.xml');
            end
            g = fopen(sprintf('xmls/%s.xml', name_base), 'w');
            while ~feof(f)
                line = fgets(f);
                if ~isempty(strfind(line, 'insert_sequences'))
                    %% print the loci
                    for l = 1 : length(loci)
                        fprintf(g, '\t<data id="%s" name="alignment">\n', strrep(loci(l).name, '.fa',''));
                        fas = fastaread(['sequencestmp/' loci(l).name]);
                        for k = 1 : length(fas)
                            fprintf(g, '\t\t<sequence id="seq_%s_%s" taxon="%s" totalcount="4" value="%s"/>\n', strrep(loci(l).name, '.fa',''), fas(k).Header, fas(k).Header, fas(k).Sequence);
                        end
                        fprintf(g, '\t</data>\n');
                        loci_length(l) = length(fas(1).Sequence);
                    end
                elseif ~isempty(strfind(line, 'overallChainLength'))
                    tmp_line = strrep(line, 'temperatureScalerInput', num2str(temp(subset))); 
                    tmp_line3 = strrep(tmp_line, 'resampleEveryInput', num2str(10000));
                    tmp_line3 = strrep(tmp_line3, 'logHeatedChains="true"', 'logHeatedChains="true"');
                    fprintf(g, '%s', strrep(tmp_line3, 'overallChainLength', num2str(chainLength(subset))));
                elseif ~isempty(strfind(line, 'insert_init_tree'))
                    %% insert starting trees
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t\t<tree id="Tree.t:%s" name="stateNode">\n', strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t<taxonset id="TaxonSet.%s" spec="TaxonSet">\n', strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t\t<alignment idref="%s"/>\n', strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t</taxonset>\n');
                        fprintf(g, '\t\t\t\t</tree>\n');

                    end

                elseif ~isempty(strfind(line, 'insert_genetreeconstrait'))
                     %% insert gene tree constraints
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t\t<distribution id="genetreeconstraint.%s" monophyletic="true"  spec="beast.math.distributions.MRCAPrior" tree="@Tree.t:%s">\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t<taxonset id="ph1.%s" spec="TaxonSet">\n', strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t<taxon idref="AaraD1_AaraD1" spec="Taxon"/>\n');
                        fprintf(g, '\t\t\t\t\t<taxon idref="AmelC2_AmelC2" spec="Taxon"/>\n');
                        fprintf(g, '\t\t\t\t\t<taxon idref="AgamS1_AgamS1" spec="Taxon"/>\n');
                        fprintf(g, '\t\t\t\t\t<taxon idref="AmerM2_AmerM2" spec="Taxon"/>\n');
                        fprintf(g, '\t\t\t\t\t<taxon idref="AcolM1_AcolM1" spec="Taxon"/>\n');
                        fprintf(g, '\t\t\t\t\t<taxon idref="AquaS1_AquaS1" spec="Taxon"/>\n');
                        fprintf(g, '\t\t\t\t</taxonset>\n');
    %                     fprintf(g, '\t\t\t\t<Uniform id="Uniform.%s" name="distr" upper="0.05"/>\n', strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t</distribution>\n');
                    end               
                elseif  ~isempty(strfind(line, 'insert_mig_prior'))
                    fprintf(g, strrep(line, 'insert_mig_prior', num2str(mig_prior(mig_prior_it))));
                elseif  ~isempty(strfind(line, 'insert_kappas'))
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t\t<parameter id="kappa.s:%s" lower="0.0" name="stateNode">2.0</parameter>\n', strrep(loci(l).name, '.fa',''));
                    end        
                elseif  ~isempty(strfind(line, 'insert_gammas'))
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t\t<parameter id="gamma.s:%s" lower="0.0" name="stateNode">1.0</parameter>\n', strrep(loci(l).name, '.fa',''));
                    end        

                elseif  ~isempty(strfind(line, 'insert_mutation_rates'))
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t\t<parameter id="mutationRate.s:%s" lower="0.0" name="stateNode">1.0</parameter>\n', strrep(loci(l).name, '.fa',''));
                    end        
                elseif  ~isempty(strfind(line, 'insert_gene_tree'))
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t\t<geneTree idref="Tree.t:%s"/>\n', strrep(loci(l).name, '.fa',''));
                    end 
                 elseif  ~isempty(strfind(line, 'insert_gene_distribution'))
                    for l = 1 : length(loci)
                        ploidy = 2;
                        if ~isempty(strfind(loci(l).name, 'chrX'))
                            ploidy = 1.5;
                        end
                        fprintf(g, '\t\t\t\t<distribution id="geneTree.t:%s" spec="starbeast2.GeneTreeWithMigration" ploidy="%.1f" populationModel="@popModelAIM.Species" tree="@Tree.t:%s"/>\n', strrep(loci(l).name, '.fa',''),ploidy, strrep(loci(l).name, '.fa',''));
                    end
                 elseif  ~isempty(strfind(line, 'insert_kappa_prior'))
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t\t<prior id="KappaPrior.s:%s" name="distribution" x="@kappa.s:%s">\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.kappa.%s" name="distr">\n', strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t\t<parameter id="RealParameter.M.%s" estimate="false" name="M">1.0</parameter>\n', strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t\t<parameter id="RealParameter.S.%s" estimate="false" name="S">1.25</parameter>\n', strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t</LogNormal>\n');
                        fprintf(g, '\t\t\t\t</prior>\n');
                    end 
                 elseif  ~isempty(strfind(line, 'insert_gamma_prior'))
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t\t<prior id="GammaPrior.s:%s" name="distribution" x="@gamma.s:%s">\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t<Exponential id="Exponential.gamma.%s" name="distr"/>\n', strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t</prior>\n');
                    end    
                 elseif  ~isempty(strfind(line, 'insert_tree_likelihood'))
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t\t<distribution id="treeLikelihood.%s" spec="TreeLikelihood" data="@%s" tree="@Tree.t:%s" useAmbiguities="true">\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t<siteModel id="SiteModel.s:%s" spec="SiteModel" gammaCategoryCount="4" shape="@gamma.s:%s"  mutationRate="@mutationRate.s:%s">\n', strrep(loci(l).name, '.fa',''),strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t\t<parameter id="proportionInvariant.s:%s" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>\n', strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t\t<substModel id="hky.s:%s" spec="HKY" kappa="@kappa.s:%s">\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t\t\t<frequencies id="empiricalFreqs.s:%s" spec="Frequencies" data="@%s"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t\t</substModel>\n');
                        fprintf(g, '\t\t\t\t\t</siteModel>\n');
                        fprintf(g, '\t\t\t\t\t<branchRateModel id="StrictClock.c:%s" spec="beast.evolution.branchratemodel.StrictClockModel">\n', strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t\t<parameter id="strictClockRate.c:%s" estimate="false" lower="0.0" name="clock.rate">1.0</parameter>\n', strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t</branchRateModel>\n');
                        fprintf(g, '\t\t\t\t</distribution>\n');
                    end           
                  elseif  ~isempty(strfind(line, 'insert_down'))           
                     for l = 1 : length(loci)
                        fprintf(g, '\t\t\t\t\t<down idref="Tree.t:%s"/>\n', strrep(loci(l).name, '.fa',''));
                     end
                  elseif  ~isempty(strfind(line, 'insert_operators'))           
                     for l = 1 : length(loci)
                        fprintf(g, '\t\t\t\t\t<operator id="TreeScaler.t:%s" spec="ScaleOperator" scaleFactor="0.95" optimise="true" tree="@Tree.t:%s" weight="3.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t<operator id="TreeRootScaler.t:%s" spec="ScaleOperator" rootOnly="true" scaleFactor="0.7"  optimise="true"  tree="@Tree.t:%s" weight="3.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t<operator id="UniformOperator.t:%s" spec="Uniform" tree="@Tree.t:%s" weight="15.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t<operator id="SubtreeSlide.t:%s" spec="SubtreeSlide" size="0.001"  optimise="true"  tree="@Tree.t:%s" weight="15.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t<operator id="Narrow.t:%s" spec="Exchange" tree="@Tree.t:%s" weight="15.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t<operator id="Wide.t:%s" spec="Exchange" isNarrow="false" tree="@Tree.t:%s" weight="15.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t\t\t\t<operator id="WilsonBalding.t:%s" spec="WilsonBalding" tree="@Tree.t:%s" weight="15.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                     end
                  elseif  ~isempty(strfind(line, 'insert_mut_scaler')) 
                        fprintf(g, '\t\t\t\t<operator id="FixMeanMutationRatesOperator" spec="DeltaExchangeOperator"  autoOptimize="true" delta="0.9" weight="%d.0">\n', length(loci));
                        for l = 1 : length(loci)
                            fprintf(g, '\t\t\t\t\t<parameter idref="mutationRate.s:%s"/>\n', strrep(loci(l).name, '.fa',''));
                        end
                        fprintf(g, '\t\t\t\t\t<weightvector id="weightparameter" spec="parameter.IntegerParameter" dimension="%d" estimate="false" lower="0" upper="0">', length(loci));
                        for l = 1 : length(loci)
                            fprintf(g, '\t%d', loci_length(l));
                        end               
                        fprintf(g, '</weightvector>\n');
                        fprintf(g, '\t\t\t\t</operator>\n');

                elseif  ~isempty(strfind(line, 'insert_kappa_scaler'))           
                  for l = 1 : length(loci)
                      fprintf(g, '\t\t\t\t<operator id="KappaScaler.s:%s" spec="ScaleOperator" parameter="@kappa.s:%s" scaleFactor="0.75" weight="1.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                  end  
                elseif  ~isempty(strfind(line, 'insert_gamma_scaler'))           
                  for l = 1 : length(loci)
                      fprintf(g, '\t\t\t\t<operator id="GammaScaler.s:%s" spec="ScaleOperator" parameter="@gamma.s:%s" scaleFactor="0.75" weight="1.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                  end                  
                elseif  ~isempty(strfind(line, 'log_tree_likelihood'))           
                    for l = 1 : length(loci)
                         fprintf(g, '\t\t\t<log idref="treeLikelihood.%s"/>\n', strrep(loci(l).name, '.fa',''));
                    end       
                elseif  ~isempty(strfind(line, 'log_tree_height'))           
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t<log id="TreeHeight.t:%s" spec="beast.evolution.tree.TreeHeightLogger" tree="@Tree.t:%s"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                    end 
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t<log id="TreeLength.t:%s" spec="starbeast2.TreeLengthLogger" tree="@Tree.t:%s"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                    end 
                elseif  ~isempty(strfind(line, 'log_kappa'))           
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t<log idref="kappa.s:%s"/>\n', strrep(loci(l).name, '.fa',''));
                    end 
                elseif  ~isempty(strfind(line, 'log_gamma'))           
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t<log idref="gamma.s:%s"/>\n', strrep(loci(l).name, '.fa',''));
                    end 

                elseif  ~isempty(strfind(line, 'log_mut'))           
                    for l = 1 : length(loci)
                        fprintf(g, '\t\t\t<log idref="mutationRate.s:%s"/>\n', strrep(loci(l).name, '.fa',''));
                    end 
                elseif  ~isempty(strfind(line, 'log_gene_tree'))           
                     for l = 1 : length(loci)
                        fprintf(g, '\t\t<logger id="treelog.t:%s" fileName="$(tree).trees" logEvery="%d" mode="tree">\n', strrep(loci(l).name, '.fa',''), round(chainLength(subset)/100));
                        fprintf(g, '\t\t\t<log id="TreeWithMetaDataLogger.t:%s" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:%s"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                        fprintf(g, '\t\t</logger>\n');
                     end
                elseif  ~isempty(strfind(line, 'autoDelayInput'))
                    fprintf(f, strrep(line, 'autoDelayInput', num2str(autodelay(subset))));
                elseif  ~isempty(strfind(line, 'logEvery="1000000"'))
                    fprintf(g, '%s', strrep(line, 'logEvery="1000000"', ['logEvery="'  num2str(round(chainLength(subset)/100)) '"']));
                else
                    fprintf(g, '%s', line);
                end
            end  
            fclose('all');

            %% make replictes of the beast xmls
            for i = 0 : nrRep-1
                f = fopen(sprintf('xmls/%s.xml',name_base));
                g = fopen(sprintf('xmls/s%d/%s_rep%d.xml',nr_samples(subset),name_base,i), 'w');
                while ~feof(f)
                    line = fgets(f);
                    if ~isempty(strfind(line, 'fileName="'))
                        fprintf(g, strrep(line, 'fileName="', ['fileName="' name_base '_rep' num2str(i), '_']));
                    elseif ~isempty(strfind(line, 'stateFileName="'))
                        fprintf(g, strrep(line, 'stateFileName="', ['stateFileName="' name_base '_rep' num2str(i) '_']));

                    else
                        fprintf(g, line);
                    end
                end
                fclose(f);
                fclose(g);
            end
            system(sprintf('rm xmls/%s.xml', name_base));
        end
    end
end
system('rm -r sequencestmp');
end
