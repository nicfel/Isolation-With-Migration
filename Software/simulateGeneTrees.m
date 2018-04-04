% simulates gene trees from species trees

nr_samples = [1 0 0 0 1];
speciation = [0.5 3 5 10]*100000;
daughters = [1 2;
            2 3;
            1 2;
            1 2];
        
timeStep = 0.1;


height = zeros(0,0);
gNe = [1 5 10]; % define the range of gene Ne's
gMig = [0.01 0.05 0.01]; % define the range of gene migration rates

for ne = 1 %: length(gNe)
    for m = 1 %: length(gMig)
        migration_rate = gMig(m); % define the migration rate for the gene trees
        geneNe = gNe(ne); % define the Ne for the species tree
        g = fopen(sprintf('BPP_Ne_%.1f_m_%.2f.log',geneNe,migration_rate),'w');
        fprintf(g, 'Sample\tposterior\tspeceis_height\tgene_height1\n');

        for i = 1 : 1000
            disp(i)
            nr_samples = [2 1 0];
            time = 0;
            while sum(nr_samples)>1
                clear rates
                for a = 1 : length(nr_samples)
                    for b = 1 : length(nr_samples)
                        if a~=b
                            rates(a,b) = migration_rate*nr_samples(a);
                        else
                            rates(a,b) = 1/geneNe*(nr_samples(a))*(nr_samples(a)-1)/2;
                        end
                    end
                end
                rates = rates*timeStep;
                events = poissrnd(rates);        
                while sum(events)>1
                    events = poissrnd(rates);        
                end
                for a= 1 : length(nr_samples)
                    for b = 1 : length(nr_samples)
                        if events(a,b)>0
                            if a~=b
                                nr_samples(a) = nr_samples(a)-1;
                                nr_samples(b) = nr_samples(b)+1;
                            else
                                nr_samples(a) = nr_samples(a)-1;
                            end
                        end                    
                    end
                end
                time = time+timeStep;

                lala = speciation-time;
                spec = find(abs(lala)<timeStep/1000);

                if ~isempty(spec)
                   nr_samples(end+1) = nr_samples(daughters(spec,1)) + nr_samples(daughters(spec,2));
                   nr_samples(daughters(spec,2)) = [];
                   nr_samples(daughters(spec,1)) = [];
                end        
            end
            fprintf(g,'%d\t1\t1\t%f\n',i-1,time);
            
            height(end+1) = time;
        end
    end
end


[y,x]=ksdensity(height);
plot(x,y)