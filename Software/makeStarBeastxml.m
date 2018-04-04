function [  ] = makeStarBeastxml(folder,filename, MCMClength, printEvery,...
    genes, species, geneTree, speciesTree,...
    estimateGeneTree, estimateSpeciesTree, estimateMigration,...
    estimatePopulation, useNrGenes, estimateAllNe)
%Uses input arguments such as genes and species to make an IM input xml
%   makes an IM input xml with the name of filename

f = fopen([folder '/' filename '.xml'], 'w');

%% print header
fprintf(f, '<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate=''StarBeast2'' beautistatus=''noAutoSetClockRate'' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" required="starbeast2 v0.13.5" version="2.4">\n\n');

%% print the gene names and sequences
geneNames = fieldnames(genes);
for i = 1 : useNrGenes
    fprintf(f, '\n');
    fprintf(f, '\t<data id="%s.sequences" name="alignment">\n',geneNames{i});
    for j = 1 : length(genes.(geneNames{i}).name)
        if ~isempty(find(ismember(species, genes.(geneNames{i}).species{j})))
            fprintf(f, '\t\t<sequence id="%s.%s.%d" taxon="%s" totalcount="4" value="%s"/>\n',...
                geneNames{i},genes.(geneNames{i}).name{j},j, genes.(geneNames{i}).name{j}, genes.(geneNames{i}).sequence{j});
        end
    end
    fprintf(f, '\t</data>\n');
end

%% print the species names


fprintf(f, '\n\t<map name="Uniform" >beast.math.distributions.Uniform</map>\n');
fprintf(f, '\t<map name="Exponential" >beast.math.distributions.Exponential</map>\n');
fprintf(f, '\t<map name="LogNormal" >beast.math.distributions.LogNormalDistributionModel</map>\n');
fprintf(f, '\t<map name="Normal" >beast.math.distributions.Normal</map>\n');
fprintf(f, '\t<map name="Beta" >beast.math.distributions.Beta</map>\n');
fprintf(f, '\t<map name="Gamma" >beast.math.distributions.Gamma</map>\n');
fprintf(f, '\t<map name="LaplaceDistribution" >beast.math.distributions.LaplaceDistribution</map>\n');
fprintf(f, '\t<map name="prior" >beast.math.distributions.Prior</map>\n');
fprintf(f, '\t<map name="speciesPrior" >beast.isolationWithMigration.distribution.SpeciesNePrior</map>\n');
fprintf(f, '\t<map name="InverseGamma" >beast.math.distributions.InverseGamma</map>\n');
fprintf(f, '\t<map name="OneOnX" >beast.math.distributions.OneOnX</map>\n\n');

% %% print gene to species mapping
% for i = 1 : useNrGenes
%     fprintf(f, '\n');
%     fprintf(f, '\t<traitSet id="species.map.%s" spec="beast.evolution.tree.TraitSet" traitname="species">\n',geneNames{i});
%     for j = 1 : (length(genes.(geneNames{i}).name)-1)
%         fprintf(f, '%s=%s,',genes.(geneNames{i}).name{j},genes.(geneNames{i}).species{j});
%     end
%     fprintf(f, '%s=%s\n',genes.(geneNames{i}).name{end},genes.(geneNames{i}).species{end});    
%     fprintf(f, '\t\t<taxa spec="TaxonSet" alignment=''@%s.sequences''/>\n',geneNames{i});
%     fprintf(f, '\t</traitSet>\n');
% end

% print MCMC specification
fprintf(f, '\t<run id="mcmc" spec="MCMC" chainLength="%d" numInitializationAttempts="1000">\n', MCMClength);


%% print state nodes
fprintf(f, '\t\t<state id="state" storeEvery="%d">\n', printEvery);
fprintf(f, '\t\t\t<stateNode id="Tree.t:Species" spec="starbeast2.SpeciesTree">\n');
fprintf(f, '\t\t\t\t<taxonset id="taxonsuperset" spec="TaxonSet">\n');

% get unique taxa
taxa = cell(0,0);
for i = 1 : length(geneNames)
    for j = 1 : length(genes.(geneNames{i}).name)
        if ~isempty(find(ismember(species, genes.(geneNames{i}).species{j})))
            taxa{end+1} = genes.(geneNames{i}).name{j};
        end
    end
end

taxa = unique(taxa);
taxa = sort(taxa);
taxonset = cell(0);
for i = 1 : length(taxa)
    tmp = strsplit(taxa{i},'_');
    ind = find(ismember(species,tmp{1}));  
    taxonset{end+1,ind} = taxa{i};
end


for i = 1 : size(taxonset,2)
    fprintf(f, '\t\t\t\t\t<taxon id="%s" spec="TaxonSet">\n',species{i});
    for j = 1 : size(taxonset,1)
        if ~isempty(taxonset{j,i})
            fprintf(f, '\t\t\t\t\t\t<taxon id="%s" spec="Taxon"/>\n',taxonset{j,i});
        end
    end
    fprintf(f, '\t\t\t\t\t</taxon>\n');
end


fprintf(f, '\t\t\t\t</taxonset>\n');
fprintf(f, '\t\t\t</stateNode>\n');
fprintf(f, '\t\t\t<parameter id="speciationRate.t:Species" lower="0.0" name="stateNode">1.0</parameter>\n');
for i = 1 : useNrGenes
    fprintf(f, '\t\t\t<tree id="Tree.%s" name="stateNode">\n', geneNames{i});
    fprintf(f, '\t\t\t\t<taxonset id="TaxonSet.%s" spec="TaxonSet">\n', geneNames{i});
    fprintf(f, '\t\t\t\t\t<alignment idref="%s.sequences"/>\n', geneNames{i});
    fprintf(f, '\t\t\t\t</taxonset>\n');
    fprintf(f, '\t\t\t</tree>\n');
end
fprintf(f, '\t\t\t<parameter id="constPopSizes.Species" lower="0.0" name="stateNode">1.0</parameter>\n');
fprintf(f, '\t\t\t<parameter id="constPopMean.Species" lower="0.0" name="stateNode">1.0</parameter>\n');
fprintf(f, '\t\t\t<parameter id="kappa.s:26" lower="0.0" name="stateNode">2.0</parameter>\n');
fprintf(f, '\t\t\t<parameter id="freqParameter.s:26" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n');
fprintf(f, '\t\t\t<parameter id="clockRate" name="stateNode">1.0</parameter>\n');
for i = 1 : useNrGenes
   fprintf(f, '\t\t\t<parameter id="kappa.%s" lower="0.0" name="stateNode">6.0</parameter>\n', geneNames{i});
   fprintf(f, '\t\t\t<parameter id="mutationRate.%s" name="stateNode">1</parameter>\n', geneNames{i});
   fprintf(f, '\t\t\t<parameter id="freqParameter.%s" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n', geneNames{i});
end
fprintf(f, '\t\t</state>\n');

%% initialize gene trees
fprintf(f, '\t\t<init id="SBI" spec="starbeast2.StarBeastInitializer" birthRate="@speciationRate.t:Species" estimate="false" speciesTree="@Tree.t:Species">\n');

for i = 1 : useNrGenes
fprintf(f, '\t\t\t<geneTree idref="Tree.%s"/>\n', geneNames{i});
end
fprintf(f, '\t\t\t<populationModel id="popModelBridge.Species" spec="starbeast2.PassthroughModel">\n');
fprintf(f, '\t\t\t\t<childModel id="constPopModel.Species" spec="starbeast2.ConstantPopulations" populationSizes="@constPopSizes.Species" speciesTree="@Tree.t:Species"/>\n');
fprintf(f, '\t\t\t</populationModel>\n');
fprintf(f, '\t\t</init>\n');

fprintf(f, '\t\t<distribution id="posterior" spec="util.CompoundDistribution">\n');

fprintf(f, '\t\t\t<distribution id="speciescoalescent" spec="starbeast2.MultispeciesCoalescent">\n');
for i = 1 : useNrGenes
    fprintf(f, '\t\t\t\t<distribution id="geneTree.t:%s" spec="starbeast2.GeneTree" populationModel="@popModelBridge.Species" speciesTree="@Tree.t:Species" tree="@Tree.%s"/>\n', geneNames{i}, geneNames{i});
end
fprintf(f, '\t\t\t</distribution>\n');
fprintf(f, '\t\t\t<distribution id="prior" spec="util.CompoundDistribution">\n');
fprintf(f, '\t\t\t\t<distribution id="YuleModel.t:Species" spec="beast.evolution.speciation.YuleModel" birthDiffRate="@speciationRate.t:Species" tree="@Tree.t:Species"/>\n');
fprintf(f, '\t\t\t\t<prior id="constPopMeanPrior.Species" name="distribution" x="@constPopMean.Species">\n');
fprintf(f, '\t\t\t\t\t<OneOnX id="OneOnX" name="distr"/>\n');
fprintf(f, '\t\t\t\t</prior>\n');
fprintf(f, '\t\t\t\t<prior id="constPopSizesPrior.Species" name="distribution" x="@constPopSizes.Species">\n');
fprintf(f, '\t\t\t\t\t<Gamma id="Gamma" beta="@constPopMean.Species" mode="ShapeMean" name="distr">\n');
fprintf(f, '\t\t\t\t\t\t<parameter id="constPopShape.Species" estimate="false" lower="0.0" name="alpha">2.0</parameter>\n');
fprintf(f, '\t\t\t\t\t</Gamma>\n');
fprintf(f, '\t\t\t\t</prior>\n');
fprintf(f, '\t\t\t\t<prior id="speciationRatePrior.t:Species" name="distribution" x="@speciationRate.t:Species">\n');
fprintf(f, '\t\t\t\t\t<Exponential id="Exponential.3" name="distr">\n');
fprintf(f, '\t\t\t\t\t\t<parameter id="RealParameter.2" estimate="false" name="mean">2.0</parameter>\n');
fprintf(f, '\t\t\t\t\t</Exponential>\n');
fprintf(f, '\t\t\t\t</prior>\n');
for i = 1 : useNrGenes
    fprintf(f, '\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@kappa.%s">\n', geneNames{i});
    fprintf(f, '\t\t\t\t\t\t<distr spec="beast.math.distributions.LogNormalDistributionModel" M="1" S="1.25"/>\n');
    fprintf(f, '\t\t\t\t</distribution>\n');
end
fprintf(f, '\t\t\t</distribution>\n');
fprintf(f, '\t\t\t<distribution id="likelihood" spec="util.CompoundDistribution">\n');

for i = 1 : useNrGenes
    fprintf(f, '\t\t\t\t<distribution id="treeLikelihood.%s" spec="ThreadedTreeLikelihood" data="@%s.sequences" tree="@Tree.%s">\n', geneNames{i}, geneNames{i}, geneNames{i});
    fprintf(f, '\t\t\t\t\t<siteModel id="SiteModel:%s" mutationRate="@mutationRate.%s" spec="SiteModel">\n', geneNames{i}, geneNames{i});
    fprintf(f, '\t\t\t\t\t\t<substModel id="hky.s:%s" spec="HKY" kappa="@kappa.%s">\n', geneNames{i}, geneNames{i});
    fprintf(f, '\t\t\t\t\t\t\t<frequencies id="estimatedFreqs.%s" spec="Frequencies" frequencies="@freqParameter.%s"/>\n', geneNames{i}, geneNames{i});
    fprintf(f, '\t\t\t\t\t\t</substModel>\n');
    fprintf(f, '\t\t\t\t\t</siteModel>\n');
    fprintf(f, '\t\t\t\t\t<branchRateModel id="StrictClock.c:%s" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate"/>\n', geneNames{i});
    fprintf(f, '\t\t\t\t</distribution>\n');
end
fprintf(f, '\t\t\t</distribution>\n');
fprintf(f, '\t\t</distribution>\n');

%% print Operators
for i = 1 : useNrGenes
    fprintf(f, '\t\t<operator id="KappaScaler.%s" spec="ScaleOperator" optimise="true" parameter="@kappa.%s" scaleFactor="0.5" weight="10"/>\n', geneNames{i}, geneNames{i});
    fprintf(f, '\t\t<operator id="FrequenciesExchanger.%s" spec="DeltaExchangeOperator" delta="0.2" weight="10">\n', geneNames{i});
    fprintf(f, '\t\t\t<parameter idref="freqParameter.%s"/>\n', geneNames{i});
    fprintf(f, '\t\t</operator>\n');
end

fprintf(f, '\t\t<operator id="FixMeanMutationRatesOperator" spec="DeltaExchangeOperator" delta="0.75" weight="500.0">\n');
for i = 1 : useNrGenes
    fprintf(f, '\t\t\t<parameter idref="mutationRate.%s"/>\n', geneNames{i});
end
fprintf(f, '\t\t\t<weightvector id="weightparameter" spec="parameter.IntegerParameter" dimension="%d" estimate="false" lower="0" upper="0">1</weightvector>\n',useNrGenes);
fprintf(f, '\t\t</operator>\n');


for i = 1 : useNrGenes
    fprintf(f, '\n');
    fprintf(f, '\t\t<operator id="TreeScaler.t:%s" spec="ScaleOperator" scaleFactor="0.95" tree="@Tree.%s" weight="3.0"/>\n', geneNames{i}, geneNames{i});
    fprintf(f, '\t\t<operator id="TreeRootScaler.t:%s" spec="ScaleOperator" rootOnly="true" scaleFactor="0.7" tree="@Tree.%s" weight="3.0"/>\n', geneNames{i}, geneNames{i});
    fprintf(f, '\t\t<operator id="UniformOperator.t:%s" spec="Uniform" tree="@Tree.%s" weight="15.0"/>\n', geneNames{i}, geneNames{i});
    fprintf(f, '\t\t<operator id="SubtreeSlide.t:%s" spec="SubtreeSlide" size="0.002" tree="@Tree.%s" weight="15.0"/>\n', geneNames{i}, geneNames{i});
    fprintf(f, '\t\t<operator id="Narrow.t:%s" spec="Exchange" tree="@Tree.%s" weight="15.0"/>\n', geneNames{i}, geneNames{i});
    fprintf(f, '\t\t<operator id="Wide.t:%s" spec="Exchange" isNarrow="false" tree="@Tree.%s" weight="15.0"/>\n', geneNames{i}, geneNames{i});
    fprintf(f, '\t\t<operator id="WilsonBalding.t:%s" spec="WilsonBalding" tree="@Tree.%s" weight="15.0"/>\n', geneNames{i}, geneNames{i});
end


if estimateSpeciesTree
fprintf(f, '\n');

fprintf(f, '\t\t<operator id="Reheight.t:Species" spec="starbeast2.NodeReheight2" taxonset="@taxonsuperset" tree="@Tree.t:Species" weight="75.0">\n');
for i = 1 : useNrGenes
    fprintf(f, '\t\t\t<geneTree idref="geneTree.t:%s"/>\n', geneNames{i});
end
fprintf(f, '\t\t</operator>\n');

fprintf(f, '\t\t<operator id="coordinatedUniform.t:Species" spec="starbeast2.CoordinatedUniform" speciesTree="@Tree.t:Species" weight="15.0">\n');
for i = 1 : useNrGenes
    fprintf(f, '\t\t\t<geneTree idref="Tree.%s"/>\n', geneNames{i});
end
fprintf(f, '\t\t</operator>\n');

fprintf(f, '\t\t<operator id="coordinatedExponential.t:Species" spec="starbeast2.CoordinatedExponential" speciesTree="@Tree.t:Species" weight="15.0">\n');
for i = 1 : useNrGenes
    fprintf(f, '\t\t\t<geneTree idref="Tree.%s"/>\n', geneNames{i});
end
fprintf(f, '\t\t</operator>\n');

fprintf(f, '\t\t<operator id="TreeScaler.t:Species" spec="ScaleOperator" scaleFactor="0.95" tree="@Tree.t:Species" weight="3.0"/>\n');
fprintf(f, '\t\t<operator id="TreeRootScaler.t:Species" spec="ScaleOperator" rootOnly="true" scaleFactor="0.7" tree="@Tree.t:Species" weight="3.0"/>\n');
fprintf(f, '\t\t<operator id="UniformOperator.t:Species" spec="Uniform" tree="@Tree.t:Species" weight="15.0"/>\n');
fprintf(f, '\t\t<operator id="SubtreeSlide.t:Species" spec="SubtreeSlide" size="0.002" tree="@Tree.t:Species" weight="15.0"/>\n');
fprintf(f, '\t\t<operator id="Narrow.t:Species" spec="Exchange" tree="@Tree.t:Species" weight="15.0"/>\n');
fprintf(f, '\t\t<operator id="Wide.t:Species" spec="Exchange" isNarrow="false" tree="@Tree.t:Species" weight="15.0"/>\n');
fprintf(f, '\t\t<operator id="WilsonBalding.t:Species" spec="WilsonBalding" tree="@Tree.t:Species" weight="15.0"/>\n');

fprintf(f, '\t\t<operator id="updownAll:Species" spec="UpDownOperator" scaleFactor="0.75" weight="6.0">\n');
fprintf(f, '\t\t\t<up idref="speciationRate.t:Species"/>\n');
fprintf(f, '\t\t\t<down idref="Tree.t:Species"/>\n');
for i = 1 : useNrGenes
    fprintf(f, '\t\t\t<down idref="Tree.%s"/>\n', geneNames{i});
end
fprintf(f, '\t\t\t<down idref="constPopSizes.Species"/>\n');
fprintf(f, '\t\t\t<down idref="constPopMean.Species"/>\n');
fprintf(f, '\t\t</operator>\n');
fprintf(f, '\t\t<operator id="speciationRateScale.t:Species" spec="ScaleOperator" parameter="@speciationRate.t:Species" scaleFactor="0.5" weight="1.0"/>\n');
%% print loggers
for i = 1 : useNrGenes
    fprintf(f, '\n');
    fprintf(f, '\t\t<logger id="%s.Logger" fileName="%s.%s.trees" logEvery="%d" mode="tree">\n', geneNames{i}, filename, geneNames{i}, printEvery);
    fprintf(f, '\t\t\t <log id="%s.Tree.Logger" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.%s"/>\n', geneNames{i}, geneNames{i});
    fprintf(f, '\t\t</logger>\n');
end
fprintf(f, '\t\t<logger id="speciesTreeLogger" fileName="%s.species.trees" logEvery="%d" mode="tree">\n', filename, printEvery);
fprintf(f, '\t\t\t<log id="SpeciesTreeLoggerX" spec="starbeast2.SpeciesTreeLogger" populationmodel="@constPopModel.Species" speciesTree="@Tree.t:Species"/>\n');
fprintf(f, '\t\t</logger>\n');



fprintf(f, '\t\t<logger id="tracelog" fileName="%s.log" logEvery="%d" model="@posterior" sanitiseHeaders="true" sort="smart">\n', filename, printEvery);
fprintf(f, '\t\t\t<log idref="posterior"/>\n');
fprintf(f, '\t\t\t<log idref="prior"/>\n');
for i = 1 : useNrGenes
    fprintf(f, '\t\t\t<log idref="treeLikelihood.%s"/>\n', geneNames{i});
end

for i = 1 : useNrGenes
    fprintf(f, '\t\t\t<log idref="kappa.%s"/>\n', geneNames{i});
    fprintf(f, '\t\t\t<log idref="freqParameter.%s"/>\n', geneNames{i});
    fprintf(f, '\t\t\t<log idref="mutationRate.%s"/>\n', geneNames{i});
end
for i = 1 : useNrGenes
    fprintf(f, '\t\t\t<log id="Height.%s" spec="beast.evolution.tree.TreeHeightLogger" tree="@Tree.%s"/>\n', geneNames{i}, geneNames{i});
end

fprintf(f, '\t\t\t<log idref="speciescoalescent"/>\n');
fprintf(f, '\t\t\t<log idref="speciationRate.t:Species"/>\n');
fprintf(f, '\t\t\t<log idref="YuleModel.t:Species"/>\n');
fprintf(f, '\t\t\t<log idref="constPopMean.Species"/>\n');

fprintf(f, '\t\t</logger>\n');
fprintf(f, '\t\t<logger id="screenlog" logEvery="%d">\n', printEvery);
fprintf(f, '\t\t\t<log idref="posterior"/>\n');
fprintf(f, '\t\t</logger>\n');

fprintf(f, '\t</run>\n');
fprintf(f, '</beast>\n');
fclose(f);
end

