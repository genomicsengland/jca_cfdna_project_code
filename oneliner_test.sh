#get lp_numbers for cfdna samples and print the ones that are also in the report log
for lp_number in $(awk -F'[,\t]' 'FNR==NR{a[$1]++; next} {if($1 in a){print $3}}' ../cfdna_samples.csv ../reports_log.20092017.csv); do
        #print the lp_number
        echo $lp_number;
        #use the lp_number to pull out the tiering json file that you generated and search for geneNames
        awk '$0~/geneName/{print $2}' $(grep $lp_number ../tiering_files_final | awk -F, '{print $2}') | sed 's/[",]//g' | \
        #and only print gene names that are in the thermo_gene_panel
        awk 'FNR==NR{a[$0]++;next} {if($0 in a){print $0}}' - ../thermo_gene_panel; done | \
tr '\n' '\t' | tr '\t\n' ' ' | sed 's/LP/\nLP/g' | awk 'NF>1{print $0}'
