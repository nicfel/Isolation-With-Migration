function [  ] = makeAIMxml(folder,filename, MCMClength, printEvery,...
    genes, species, geneTree, speciesTree,...
    estimateGeneTree, estimateSpeciesTree,...
    useNrGenes, isBSSVS)
%Uses input arguments such as genes and species to make an IM input xml
%   makes an IM input xml with the name of filename

f = fopen([folder '/' filename '.xml'], 'w');

%% print header
fprintf(f, '<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate=''StarBeast2'' beautistatus=''noAutoSetClockRate'' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" required="starbeast2 v0.13.0" version="2.4">\n\n');

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
fprintf(f, '\n');
fprintf(f, '\t<data id="species" name="alignment">\n');
for i = 1 : length(species)
    fprintf(f, '\t\t<sequence id="species.%s" taxon="%s" totalcount="4" value="??"/>\n',...
        species{i},species{i});
end
fprintf(f, '\t</data>\n');


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

% print MCMC specification
fprintf(f, '\t<run id="mcmc" spec="MCMC" chainLength="%d" numInitializationAttempts="1000">\n', MCMClength);

%% print state nodes
fprintf(f, '\t\t<state id="state" storeEvery="%d">\n', printEvery);
fprintf(f, '\t\t\t<stateNode id="species.tree" spec="starbeast2.SpeciesTree">\n');
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
for i = 1 : useNrGenes
    fprintf(f, '\t\t\t<tree id="Tree.%s" name="stateNode">\n',geneNames{i});
    fprintf(f, '\t\t\t\t<taxonset spec="TaxonSet" alignment=''@%s.sequences''/>\n',geneNames{i});
    fprintf(f, '\t\t\t</tree>\n');
end
if ~estimateGeneTree
    fprintf(f, '\t\t\t<parameter id="NeMean" name="stateNode" dimension="1" lower="0">1</parameter>\n');
else
    fprintf(f, '\t\t\t<parameter id="NeMean" name="stateNode" dimension="1" lower="0">0.0001</parameter>\n');
end
em = lognrnd(-2.9957,1,1);
while em > 0.195
    em = lognrnd(-2.9957,1,1);
end
if isBSSVS
    fprintf(f, '\t\t\t<parameter id="eM" name="stateNode" dimension="1" upper="2">1</parameter>\n');
else
    fprintf(f, '\t\t\t<parameter id="eM" name="stateNode" dimension="1" upper="1.0">%.3f</parameter>\n', em);
end
fprintf(f, '\t\t\t<parameter id="m" name="stateNode" dimension="0" lower="0">1</parameter>\n');
if ~estimateGeneTree
    fprintf(f, '\t\t\t<parameter id="Ne" name="stateNode" dimension="0" lower="0">1</parameter>\n');
else
    fprintf(f, '\t\t\t<parameter id="Ne" name="stateNode" dimension="0" lower="0">0.0001</parameter>\n');
end
fprintf(f, '\t\t\t<parameter id="bddiff" name="stateNode" dimension="0" lower="0">%.3f</parameter>\n', normrnd(100,10,1));

if isBSSVS
    fprintf(f, '\t\t\t<stateNode id="migration.indicator" spec="parameter.BooleanParameter" dimension="0">1</stateNode>');
end

if estimateGeneTree
    if isBSSVS
        fprintf(f, '\t\t\t<parameter id="kappa" lower="0.0" name="stateNode">6.0</parameter>\n');
        for i = 1 : useNrGenes
           fprintf(f, '\t\t\t<parameter id="clockRate.%s" name="stateNode">1</parameter>\n', geneNames{i});
        end
    else
        for i = 1 : useNrGenes
           fprintf(f, '\t\t\t<parameter id="kappa.%s" lower="0.0" name="stateNode">6.0</parameter>\n', geneNames{i});
           fprintf(f, '\t\t\t<parameter id="clockRate.%s" name="stateNode">1</parameter>\n', geneNames{i});
           fprintf(f, '\t\t\t<parameter id="freqParameter.%s" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n', geneNames{i});
        end
    end
end
fprintf(f, '\t\t</state>\n');

%% initialize gene trees
if isempty(geneTree)
    for i = 1 : useNrGenes
        fprintf(f, '\t\t<init id="RandomTree.t:%s" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.%s" taxa="@%s.sequences">\n',geneNames{i},geneNames{i},geneNames{i});
        fprintf(f, '\t\t\t<populationModel id="ConstantPopulation.%s" spec="ConstantPopulation">\n',geneNames{i});
        fprintf(f, '\t\t\t\t<parameter id="randomPopSize.%s" name="popSize">0.001</parameter>\n',geneNames{i});
        fprintf(f, '\t\t\t</populationModel>\n');
        fprintf(f, '\t\t</init>\n');
    end
else
    for i = 1 : useNrGenes
        fprintf(f, '\t\t<init spec="beast.util.TreeParser" id="NewickTree.t:%s" adjustTipHeights="true"\n', geneNames{i});
        fprintf(f, '\t\tinitial="@Tree.%s" taxa="@%s.sequences" IsLabelledNewick="true"\n', geneNames{i}, geneNames{i});
        fprintf(f, '\t\tnewick="%s"/>\n', geneTree{i});
    end
end

%% initialize species trees
if estimateSpeciesTree
    fprintf(f, '<!-- true species tree: %s -->\n', speciesTree);
    fprintf(f, '\t\t<init id="RandomTree.t:Species" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@species.tree" taxa="@species">\n');
    fprintf(f, '\t\t\t<populationModel id="ConstantPopulation.Species" spec="ConstantPopulation">\n');
    fprintf(f, '\t\t\t\t<parameter id="randomPopSize.Species" name="popSize">0.001</parameter>\n');
    fprintf(f, '\t\t\t</populationModel>\n');
    fprintf(f, '\t\t</init>\n');
else
    fprintf(f, '\t\t<init spec="beast.util.TreeParser" id="NewickTree.t:Species" adjustTipHeights="true"\n');
    fprintf(f, '\t\tinitial="@species.tree" taxa="@species" IsLabelledNewick="true"\n');
    fprintf(f, '\t\tnewick="%s"/>\n', speciesTree);
end

% print distributions
fprintf(f, '\t\t<distribution id="posterior" spec="util.CompoundDistribution">\n');
fprintf(f, '\t\t\t<distribution id="prior" spec="util.CompoundDistribution">\n');
fprintf(f, '\t\t\t\t<distribution id="YuleModel.speciesTree" spec="beast.evolution.speciation.YuleModel" birthDiffRate="@bddiff" tree="@species.tree"/>\n');
fprintf(f, '\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@Ne">\n');
fprintf(f, '\t\t\t\t\t<distr spec="beast.math.distributions.LogNormalDistributionModel" meanInRealSpace="true" S="0.25">\n');
fprintf(f, '\t\t\t\t\t\t<M idref="NeMean"/>\n');
fprintf(f, '\t\t\t\t\t</distr>\n');
fprintf(f, '\t\t\t\t</distribution>\n');
fprintf(f, '\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@NeMean">\n');
fprintf(f, '\t\t\t\t\t<distr spec="beast.math.distributions.OneOnX"/>\n');
fprintf(f, '\t\t\t\t</distribution>\n');
if isBSSVS
    fprintf(f, '\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@m">\n');
    fprintf(f, '\t\t\t\t\t<distr spec="beast.math.distributions.Exponential"/>\n');
    fprintf(f, '\t\t\t\t</distribution>\n');
    fprintf(f, '\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@eM">\n');
    fprintf(f, '\t\t\t\t\t<distr spec="beast.math.distributions.Exponential" mean="0.05"/>\n');
    fprintf(f, '\t\t\t\t</distribution>\n');
else
    fprintf(f, '\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@m">\n');
    fprintf(f, '\t\t\t\t\t<distr spec="beast.math.distributions.Exponential"/>\n');
    fprintf(f, '\t\t\t\t</distribution>\n');
    fprintf(f, '\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@eM">\n');
    fprintf(f, '\t\t\t\t\t<distr spec="beast.math.distributions.Exponential" mean="0.1"/>\n');
    fprintf(f, '\t\t\t\t</distribution>\n');
end

if isBSSVS
    fprintf(f, '\t\t\t\t<prior id="nonZeroRatePrior.Migration" name="distribution">\n');
    fprintf(f, '\t\t\t\t\t<x id="nonZeroMigrationRates" spec="util.Sum">\n');
    fprintf(f, '\t\t\t\t\t\t<arg idref="migration.indicator"/>\n');
    fprintf(f, '\t\t\t\t\t</x>\n');
%     fprintf(f, '\t\t\t\t\t<distr spec="beast.math.distributions.Poisson" lambda="2"/>\n');

    fprintf(f, '\t\t\t\t\t<distr spec="beast.math.distributions.Uniform" lower="0.0" upper="%d"/>\n', (length(species)-1)*(length(species)-1)*2);
    fprintf(f, '\t\t\t\t</prior>\n');
end

fprintf(f, '\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@bddiff">\n');
fprintf(f, '\t\t\t\t\t\t<distr spec="beast.math.distributions.OneOnX"/>\n');
fprintf(f, '\t\t\t\t</distribution>\n');


if estimateGeneTree
    if isBSSVS
        fprintf(f, '\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@migration.indicator">\n');
        fprintf(f, '\t\t\t\t\t\t<distr spec="beast.math.distributions.Exponential" mean="0.339623"/>\n');
        fprintf(f, '\t\t\t\t</distribution>\n');
        fprintf(f, '\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@kappa">\n');
        fprintf(f, '\t\t\t\t\t\t<distr spec="beast.math.distributions.LogNormalDistributionModel" M="1" S="1.25"/>\n');
        fprintf(f, '\t\t\t\t</distribution>\n');
    else
        for i = 1 : useNrGenes
            fprintf(f, '\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@kappa.%s">\n', geneNames{i});
            fprintf(f, '\t\t\t\t\t\t<distr spec="beast.math.distributions.LogNormalDistributionModel" M="1" S="1.25"/>\n');
            fprintf(f, '\t\t\t\t</distribution>\n');
        end
    end
end
fprintf(f, '\t\t\t\t<distribution id="speciescoalescent" spec="starbeast2.MultispeciesCoalescent">\n');
fprintf(f, '\t\t\t\t\t<distribution id="geneTree.t:%s" spec="starbeast2.GeneTreeWithMigration" tree="@Tree.%s">\n', geneNames{1}, geneNames{1});
fprintf(f, '\t\t\t\t\t\t<populationModel id="popModelBridge.Species" spec="starbeast2.ConstantWithGeneFlow">\n');
fprintf(f, '\t\t\t\t\t\t\t<Ne idref="Ne"/>\n');
fprintf(f, '\t\t\t\t\t\t\t<m idref="m"/>\n');
if isBSSVS
    fprintf(f, '\t\t\t\t\t\t\t<indicator idref="migration.indicator"/>\n');
end
% if isBSSVS
%     fprintf(f, '\t\t\t\t\t\t\t<migrationModel id="migModel" spec="starbeast2.AllEqual" effectiveMigrants="@eM" speciesTree="@species.tree"/>\n');
% else
    fprintf(f, '\t\t\t\t\t\t\t<migrationModel id="migModel" spec="starbeast2.MinimalBranchLength"  minimalBranchLength="0.00005" effectiveMigrants="@eM" speciesTree="@species.tree"/>\n');
% end
fprintf(f, '\t\t\t\t\t\t</populationModel>\n');
fprintf(f, '\t\t\t\t\t</distribution>\n');


for i = 2 : useNrGenes
    fprintf(f, '\t\t\t\t\t<distribution id="geneTree.t:%s" spec="starbeast2.GeneTreeWithMigration" populationModel="@popModelBridge.Species" tree="@Tree.%s"/>\n', geneNames{i}, geneNames{i});
end
fprintf(f, '\t\t\t\t</distribution>\n');


fprintf(f, '\t\t\t</distribution>\n');
fprintf(f, '\t\t\t<distribution id="likelihood" spec="util.CompoundDistribution">\n');
if estimateGeneTree
    fprintf(f, '\t\t\t\t<distribution id="likelihood.GeneTrees" spec="util.CompoundDistribution">\n');

    for i = 1 : useNrGenes
        fprintf(f, '\t\t\t\t\t<distribution id="treeLikelihood.%s" spec="ThreadedTreeLikelihood" data="@%s.sequences" tree="@Tree.%s" useAmbiguities="true">\n', geneNames{i}, geneNames{i}, geneNames{i});
        fprintf(f, '\t\t\t\t\t\t<siteModel id="SiteModel:%s" spec="SiteModel">\n', geneNames{i});
        if isBSSVS
            fprintf(f, '\t\t\t\t\t\t\t<substModel id="hky.s:%s" spec="HKY" kappa="@kappa">\n', geneNames{i});
            fprintf(f, '\t\t\t\t\t\t\t\t<frequencies id="empiricalFreqs.%s" spec="Frequencies" data="@%s.sequences"/>\n', geneNames{i}, geneNames{i});
        else
            fprintf(f, '\t\t\t\t\t\t\t<substModel id="hky.s:%s" spec="HKY" kappa="@kappa.%s">\n', geneNames{i}, geneNames{i});
            fprintf(f, '\t\t\t\t\t\t\t\t<frequencies id="estimatedFreqs.%s" spec="Frequencies" frequencies="@freqParameter.%s"/>\n', geneNames{i}, geneNames{i});
        end
        fprintf(f, '\t\t\t\t\t\t\t</substModel>\n');
        fprintf(f, '\t\t\t\t\t\t</siteModel>\n');
        fprintf(f, '\t\t\t\t\t\t<branchRateModel id="StrictClock.c:%s" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.%s"/>\n', geneNames{i}, geneNames{i});
        fprintf(f, '\t\t\t\t\t</distribution>\n');
    end
    fprintf(f, '\t\t\t\t</distribution>\n');

end


fprintf(f, '\t\t\t</distribution>\n');
fprintf(f, '\t\t</distribution>\n');

%% print Operators
if estimateGeneTree
    if isBSSVS
        fprintf(f, '\t\t<operator id="KappaScaler" spec="ScaleOperator" optimise="true" parameter="@kappa" scaleFactor="0.5" weight="10"/>\n');
    else
        for i = 1 : useNrGenes
            fprintf(f, '\t\t<operator id="KappaScaler.%s" spec="ScaleOperator" optimise="true" parameter="@kappa.%s" scaleFactor="0.5" weight="10"/>\n', geneNames{i}, geneNames{i});
            fprintf(f, '\t\t<operator id="FrequenciesExchanger.%s" spec="DeltaExchangeOperator" delta="0.2" weight="10">\n', geneNames{i});
            fprintf(f, '\t\t\t<parameter idref="freqParameter.%s"/>\n', geneNames{i});
            fprintf(f, '\t\t</operator>\n');
        end
    end

    
    fprintf(f, '\t\t<operator id="ClockExchanger" spec="DeltaExchangeOperator" delta="0.2" autoOptimize="true" weight="100">\n');
    for i = 1 : useNrGenes
        fprintf(f, '\t\t\t<parameter idref="clockRate.%s"/>\n', geneNames{i});
    end    
    fprintf(f, '\t\t\t<weightvector id="weightparameter" spec="parameter.IntegerParameter" dimension="%d" estimate="false" lower="0" upper="0">',useNrGenes);
    for i = 1 : useNrGenes-1
        fprintf(f, '%d ', length(genes.(geneNames{i}).sequence{1,1}));
    end
    fprintf(f, '%d', length(genes.(geneNames{useNrGenes}).sequence{1,1}));
    fprintf(f, '</weightvector>\n');
    fprintf(f, '\t\t</operator>\n');


    for i = 1 : useNrGenes
        fprintf(f, '\n');
        fprintf(f, '\t\t<operator id="TreeScaler.t:%s" spec="ScaleOperator" scaleFactor="0.95" tree="@Tree.%s" weight="3.0"/>\n', geneNames{i}, geneNames{i});
        fprintf(f, '\t\t<operator id="TreeRootScaler.t:%s" spec="ScaleOperator" rootOnly="true" scaleFactor="0.7" tree="@Tree.%s" weight="3.0"/>\n', geneNames{i}, geneNames{i});
        fprintf(f, '\t\t<operator id="UniformOperator.t:%s" spec="Uniform" tree="@Tree.%s" weight="30.0"/>\n', geneNames{i}, geneNames{i});
        fprintf(f, '\t\t<operator id="SubtreeSlide.t:%s" spec="SubtreeSlide" size="0.0005" tree="@Tree.%s" weight="30.0"/>\n', geneNames{i}, geneNames{i});
        fprintf(f, '\t\t<operator id="Narrow.t:%s" spec="Exchange" tree="@Tree.%s" weight="3.0"/>\n', geneNames{i}, geneNames{i});
        fprintf(f, '\t\t<operator id="Wide.t:%s" spec="Exchange" isNarrow="false" tree="@Tree.%s" weight="3.0"/>\n', geneNames{i}, geneNames{i});
        fprintf(f, '\t\t<operator id="WilsonBalding.t:%s" spec="WilsonBalding" tree="@Tree.%s" weight="3.0"/>\n', geneNames{i}, geneNames{i});
    end
end

if estimateSpeciesTree
    fprintf(f, '\n');
    
%     fprintf(f, '\t\t<operator id="TreeScaler.t:Species" spec="ScaleOperator" scaleFactor="0.95" tree="@species.tree" weight="5"/>\n');
    fprintf(f, '\t\t<operator id="TreeRootScaler.t:Species" spec="ScaleOperator" rootOnly="true" scaleFactor="0.9" tree="@species.tree" weight="5"/>\n');
    fprintf(f, '\t\t<operator id="UniformOperator.t:Species" spec="Uniform" tree="@species.tree" weight="100"/>\n');
%     fprintf(f, '\t\t<operator id="SubtreeSlide.t:Species" spec="SubtreeSlide" size="0.00005" tree="@species.tree" weight="50"/>\n');
%     fprintf(f, '\t\t<operator id="NodeScaler.t:Species" spec="aim.operators.NodeRe" scaleFactor="0.95" tree="@species.tree" weight="15"/>\n');

    fprintf(f, '\t\t<operator id="Narrow.t:Species" spec="Exchange" tree="@species.tree" weight="15"/>\n');
    fprintf(f, '\t\t<operator id="Wide.t:Species" spec="Exchange" isNarrow="false" tree="@species.tree" weight="5"/>\n');
    fprintf(f, '\t\t<operator id="WilsonBalding.t:Species" spec="WilsonBalding" tree="@species.tree" weight="5"/>\n');

   
    fprintf(f, '\t\t<operator id="coordinatedUniform.t:Species" spec="starbeast2.CoordinatedUniform" speciesTree="@species.tree" weight="5">\n');
    for i = 1 : useNrGenes
        fprintf(f, '\t\t<geneTree idref="Tree.%s"/>\n', geneNames{i});
    end
    fprintf(f, '\t\t</operator>\n');
    
        fprintf(f, '\t\t<operator id="coordinatedExponential.t:Species" spec="starbeast2.CoordinatedExponential" speciesTree="@species.tree" weight="5">\n');
    for i = 1 : useNrGenes
        fprintf(f, '\t\t<geneTree idref="Tree.%s"/>\n', geneNames{i});
    end
    fprintf(f, '\t\t</operator>\n');
    
    fprintf(f, '\t\t<operator id="AllTreeScaler" spec="UpDownOperator" optimise="true" scaleFactor="0.7" weight="2.5">\n');
    fprintf(f, '\t\t\t<up idref="species.tree"/>\n');
    for i = 1 : useNrGenes
        fprintf(f, '\t\t\t<up idref="Tree.%s"/>\n', geneNames{i});
    end
    fprintf(f, '\t\t\t<up idref="NeMean"/>\n');
%     fprintf(f, '\t\t\t<down idref="m"/>\n');
    fprintf(f, '\t\t</operator>\n');
end

% fprintf(f, '\t\t<operator id="migrationExchanger" spec="DeltaExchangeOperator" autoOptimize="true" delta="0.2" weight="30">\n');
% fprintf(f, '\t\t\t<parameter idref="m"/>\n');
% fprintf(f, '\t\t</operator>\n');

fprintf(f, '\t\t<operator id="BirthDeathDifferenceScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@bddiff" weight="10"/>\n');

fprintf(f, '\t\t<operator id="NeScaler" spec="ScaleOperator" scaleFactor="0.95" scaleAll="true" scaleAllIndependently="true" optimise="true" parameter="@Ne" weight="30"/>\n');
if ~estimateGeneTree
    fprintf(f, '\t\t<operator id="mScaler" spec="ScaleOperator" scaleFactor="0.95" scaleAll="true" scaleAllIndependently="true"  optimise="true" parameter="@m" weight="200"/>\n');
    fprintf(f, '\t\t<operator id="NeMeanScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@NeMean" weight="5"/>\n');
    fprintf(f, '\t\t<operator id="eMScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@eM" weight="50"/>\n');
else
    if isBSSVS
        fprintf(f, '\t\t<operator id="mScaler" spec="ScaleOperator" scaleFactor="0.95" scaleAll="true" scaleAllIndependently="true"  optimise="true" parameter="@m" weight="30"/>\n');
        fprintf(f, '\t\t<operator id="NeMeanScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@NeMean" weight="5"/>\n');
        fprintf(f, '\t\t<operator id="eMScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@eM" weight="20"/>\n');
    else
        fprintf(f, '\t\t<operator id="mScaler" spec="ScaleOperator" scaleFactor="0.95" scaleAll="true" scaleAllIndependently="true"  optimise="true" parameter="@m" weight="30"/>\n');
        fprintf(f, '\t\t<operator id="NeMeanScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@NeMean" weight="5"/>\n');
        fprintf(f, '\t\t<operator id="eMScaler" spec="ScaleOperator" scaleFactor="0.8" optimise="true" parameter="@eM" weight="20"/>\n');
    end
end

if isBSSVS
    fprintf(f, '\t\t<operator id="indicatorFlip.s:indicator" spec="BitFlipOperator" parameter="@migration.indicator" uniform="true" weight="20"/>\n');
end




%% print loggers
if estimateGeneTree
    for i = 1 : useNrGenes
        fprintf(f, '\n');
        fprintf(f, '\t\t<logger id="%s.Logger" fileName="%s.%s.trees" logEvery="%d" mode="tree">\n', geneNames{i}, filename, geneNames{i}, printEvery);
        fprintf(f, '\t\t\t<log id="%s.Tree.Logger" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.%s"/>\n', geneNames{i}, geneNames{i});
        fprintf(f, '\t\t</logger>\n');
    end
end
fprintf(f, '\t\t<logger id="species.Logger" fileName="$(filebase).species.trees" logEvery="%d" mode="tree">\n', printEvery);
fprintf(f, '\t\t\t<log id="species.tree.logger" spec="starbeast2.SpeciesTreeLoggerWithGeneFlow" populationModel="@popModelBridge.Species" />\n');
fprintf(f, '\t\t</logger>\n');


fprintf(f, '\t\t<logger id="tracelog" fileName="$(filebase).log" logEvery="%d" model="@posterior" sanitiseHeaders="true" sort="smart">\n', printEvery);
fprintf(f, '\t\t\t<log idref="posterior"/>\n');
fprintf(f, '\t\t\t<log idref="likelihood"/>\n');
fprintf(f, '\t\t\t<log idref="prior"/>\n');
fprintf(f, '\t\t\t<log idref="speciescoalescent"/>\n');
if estimateGeneTree
   fprintf(f, '\t\t\t<log idref="likelihood.GeneTrees"/>\n');
end


for i = 1 : useNrGenes
%     fprintf(f, '\t\t\t<log idref="%s.AIM"/>\n', geneNames{i});    
end
if estimateGeneTree
   for i = 1 : useNrGenes
        fprintf(f, '\t\t\t<log idref="treeLikelihood.%s"/>\n', geneNames{i});
    end
end
if estimateGeneTree
    if isBSSVS
        fprintf(f, '\t\t\t<log idref="kappa"/>\n');
        for i = 1 : useNrGenes
            fprintf(f, '\t\t\t<log idref="clockRate.%s"/>\n', geneNames{i});
        end
    else
        for i = 1 : useNrGenes
            fprintf(f, '\t\t\t<log idref="kappa.%s"/>\n', geneNames{i});
            fprintf(f, '\t\t\t<log idref="freqParameter.%s"/>\n', geneNames{i});        
            fprintf(f, '\t\t\t<log idref="clockRate.%s"/>\n', geneNames{i});
        end
    end
%     fprintf(f, '\t\t\t<log id="NodeHeightLogger" spec="aim.logger.NodeHeightLogger" tree="@species.tree"/>\n');

    for i = 1 : useNrGenes
        fprintf(f, '\t\t\t<log id="Height.%s" spec="beast.evolution.tree.TreeHeightLogger" tree="@Tree.%s"/>\n', geneNames{i}, geneNames{i});
    end
end
fprintf(f, '\t\t\t<log idref="bddiff"/>\n');

fprintf(f, '\t\t\t<log idref="NeMean"/>\n');
fprintf(f, '\t\t\t<log idref="Ne"/>\n');
% fprintf(f, '\t\t\t<log idref="mMax"/>\n');
fprintf(f, '\t\t\t<log idref="eM"/>\n');
fprintf(f, '\t\t\t<log idref="m"/>\n');
if isBSSVS
    fprintf(f, '\t\t\t<log idref="nonZeroMigrationRates"/>\n');
    fprintf(f, '\t\t\t<log idref="migration.indicator"/>\n');
end

fprintf(f, '\t\t</logger>\n');
fprintf(f, '\t\t<logger id="screenlog" logEvery="%d">\n', printEvery/10);
fprintf(f, '\t\t\t<log idref="posterior"/>\n');
fprintf(f, '\t\t</logger>\n');

fprintf(f, '\t</run>\n');
fprintf(f, '</beast>\n');
fclose(f);
end

