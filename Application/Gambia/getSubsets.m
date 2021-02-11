% get the subsets
clear

% define the number of loci
nr_loci = [50 200];
% define the number of different random sets of loci
nr_sets = 1;
% define the number of parallel MCMC chains
nr_reps = 3;

% make a dir for the AIM xmls
% system('rm -r xmls_allequal');
% system('mkdir xmls_allequal');

% temperature scaler for the coupled MCMC
temperature = [0.01, 0.005, 0.0025];
temperature = [0.01, 0.0025];

% length of the coupled MCMC
chainLength = [50000000, 100000000, 50000000];
autodelay = chainLength./[5 2 2];

proportions(1,1) = length(dir('data/maffilter/chr2R/*.fa')); 
proportions(1,2) = length(dir('data/maffilter/chr3L/*.fa')); 
proportions(1,3) = length(dir('data/maffilter/chr3R/*.fa'));
proportions(1,4) = length(dir('data/maffilter/chrX/*.fa'));

max_nr(1,:) =  proportions(1,:);   
% max_nr(1,4) = length(dir('data/maffilter/chrX/*.fa'));
max_nr(2,:) =  [0 0 0 proportions(1,4)];


proportions = [0 0.125 0.125 0.75;
                0 0 0 1;];


chrName = {'all', 'chrX'};

fixed_tree = [true, false];

for i = 1 : size(proportions,1)
    weights = proportions(i,:)/sum(proportions(i,:));
    getSubsetXmls(nr_loci, nr_sets, nr_reps,...
        temperature, chainLength, autodelay, weights, max_nr(i,:),...
        chrName{i}, 3, [50], fixed_tree(i))
end
