% makes an xml for each locus on the X-Chromosome
clear loci

system('rm -r yule');
system('mkdir yule');
system('mkdir yule/xmls');

loci = dir('data/maffilter/chrX/*.fa');


for i = 1 : length(loci)
    disp(i)
    fasta = fastaread(['data/maffilter/chrX/' loci(i).name]);
    c = 1;
    for k = 1 : length(fasta)
        if isempty(strfind(fasta(k).Header, 'AgamP4')) && ...
                        isempty(strfind(fasta(k).Header, 'AepiE1'))
            newFasta(c) = fasta(k);
            c = c+1;
        end
    end
    
    fname = strrep(loci(i).name, '-','_');
    fname = strrep(fname, '.fa','.xml');

    f = fopen('templateYule.xml');
    g = fopen(['yule/xmls/' fname], 'w');
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'insert_data')
            for k = 1 : length(newFasta)
                fprintf(g, '<sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n',...
                    newFasta(k).Header, newFasta(k).Header, newFasta(k).Sequence);
            end
        elseif contains(line, 'insert_filename')
            fprintf(g, strrep(line, 'insert_filename',  strrep(fname, '.xml','')));
        else
            fprintf(g, line);
        end
    end
    fclose(f);
    fclose(g);
end
