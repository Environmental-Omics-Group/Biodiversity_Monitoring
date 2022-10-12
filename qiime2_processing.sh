# Data Processing with QIIME2
Parts of this script follow the [QIIME2 moving pictures tutorial](https://docs.qiime2.org/2021.2/tutorials/moving-pictures/)
required packages: QIIME2021.2; cutadapt 2.7; pigz 2.4; fastqc 0.11.9; multiqc 1.9; q2-picrust2 2021.2

## Demultiplexing by barcode
 'top' set of primers used for 16SV1 and 16SV4, 'bottom' set of indexes used for rbcL, COI and 18S


```
for top in $(ls ../../lr_raw_niamh/data/*top*_R1_001.fastq.gz | cut -d '/' -f5); do
        echo "now $top demux start"
        sampleid=$(echo ${top} | cut -d '_' -f1)
        namestring=$(echo ${top} | cut -d '_' -f1,2,3)
        cutadapt -e 0.07 --untrimmed-output ../partdemuxed/164_${sampleid}_R1.fastq.gz --untrimmed-paired-output ../partdemuxed/164_${sampleid}_R2.fastq.gz -j 6 -g primers161F=^NNNNNAGAGTTTGATCMTGGCTCAG -o ../demux_fwd/161/00_reads/${sampleid}_R1.fastq.gz -p ../demux_fwd/161/00_reads/${sampleid}_R2.fastq.gz ../../lr_raw_niamh/data/${namestring}_L001_R1_001.fastq.gz  ../../lr_raw_niamh/data/${namestring}_L001_R2_001.fastq.gz
        cutadapt -e 0.07 --untrimmed-output ../partdemuxed/leftover/noprimer_${sampleid}_R1.fastq.gz --untrimmed-paired-output ../partdemuxed/leftover/noprimer_${sampleid}_R2.fastq.gz -j 6 -g primers164F=^NNNNNGTGCCAGCMGCCGCGGTAA -o ../demux_fwd/164/00_reads/${sampleid}_R1.fastq.gz -p ../demux_fwd/164/00_reads/${sampleid}_R2.fastq.gz ../partdemuxed/164_${sampleid}_R1.fastq.gz  ../partdemuxed/164_${sampleid}_R2.fastq.gz

done

for btm in $(ls ../../lr_raw_niamh/data/*bottom*_R1_001.fastq.gz | cut -d '/' -f5); do
        echo "$btm demux start"
        sampleid=$(echo ${btm} | cut -d '_' -f1)
        namestring=$(echo ${btm} | cut -d '_' -f1,2,3)
        cutadapt -e 0.07 --untrimmed-output ../partdemuxed/rbcl_coi_${sampleid}_R1.fastq.gz --untrimmed-paired-output ../partdemuxed/rbcl_coi_${sampleid}_R2.fastq.gz -j 6 -g primers181F=^NNNNNGCTTGTCTCAAAGATTAAGCC -o ../demux_fwd/181/00_reads/${sampleid}_R1.fastq.gz -p ../demux_fwd/181/00_reads/${sampleid}_R2.fastq.gz ../../lr_raw_niamh/data/${namestring}_L001_R1_001.fastq.gz  ../../lr_raw_niamh/data/${namestring}_L001_R2_001.fastq.gz
        
        cutadapt -e 0.07 --untrimmed-output ../partdemuxed/coi_${sampleid}_R1.fastq.gz --untrimmed-paired-output ../partdemuxed/coi_${sampleid}_R2.fastq.gz -j 6 -g primersrbclF=^NNNNNATGCGTTGGAGAGARCGTTTC -o ../demux_fwd/rbcl/00_reads/${sampleid}_R1.fastq.gz -p ../demux_fwd/rbcl/00_reads/${sampleid}_R2.fastq.gz ../partdemuxed/rbcl_coi_${sampleid}_R1.fastq.gz  ../partdemuxed/rbcl_coi_${sampleid}_R2.fastq.gz
   
        cutadapt -e 0.07 --untrimmed-output ../partdemuxed/leftover/noprimer_${sampleid}_R1.fastq.gz --untrimmed-paired-output ../partdemuxed/leftover/noprimer_${sampleid}_R2.fastq.gz -j 6 -g primerscoiF=^NNNNNGGWACWGGWTGAACWGTWTAYCCYCC -o ../demux_fwd/coi/00_reads/${sampleid}_R1.fastq.gz -p ../demux_fwd/coi/00_reads/${sampleid}_R2.fastq.gz ../partdemuxed/coi_${sampleid}_R1.fastq.gz  ../partdemuxed/coi_${sampleid}_R2.fastq.gz

done
```
## Quality checking

```

for amp in 161 164 181 rbcl coi ; do
    mkdir -p ../01_qc/${amp}
    for i in $(ls -1 ../demux_fwd/${amp}/00_reads/*.gz) ; do
        fastqc -o ../01_qc/${amp} -t 6 $i
    done
    multiqc -n ../01_qc/multiqc_report_${amp}.html ../01_qc/${amp}/*
done 
```
## Import to QIIME2

```

for amp in 161 164 181 rbcl coi 
do
    mkdir -p ../demux_fwd/${amp}/01_qiimeimport

    qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path ./manifests/${amp}_manifest.txt \
    --output-path ../demux_fwd/${amp}/01_qiimeimport/paired-end-demux-${amp}.qza \
    --input-format PairedEndFastqManifestPhred33V2

    qiime demux summarize \
    --i-data ../demux_fwd/${amp}/01_qiimeimport/paired-end-demux-${amp}.qza \
    --o-visualization ../demux_fwd/${amp}/01_qiimeimport/paired-end-demux-${amp}.qzv

done
```
## Denoising
Trimming parameters are specific to each amplicon and chosen based on paired-end-demux-${amp}.qzv quality plots, supported by fastqc and multiqc reports

```
for amp in 161 164 181 rbcl coi 
do
    if [ ${amp} = 161 ] ; then
        leftf=0
        leftr=23
        lenf=224
        lenr=233
    elif [ ${amp} = 164 ] ; then
        leftf=2
        leftr=21
        lenf=225
        lenr=247
    elif [ ${amp} = 181 ] ; then
        leftf=0
        leftr=19
        lenf=223
        lenr=238
    elif [ ${amp} = rbcl ] ; then
        leftf=10
        leftr=27
        lenf=224
        lenr=237
    elif [ ${amp} = coi ] ; then
        leftf=0
        leftr=26
        lenf=214
        lenr=232
    fi
    echo "leftf" 
    echo ${leftf}


    qiime dada2 denoise-paired \
    --i-demultiplexed-seqs ../demux_fwd/${amp}/01_qiimeimport/paired-end-demux-${amp}.qza \
    --p-trim-left-f ${leftf} \
    --p-trim-left-r ${leftr} \
    --p-trunc-len-f ${lenf} \
    --p-trunc-len-r ${lenr} \
    --p-n-threads 6 \
    --o-table ../demux_fwd/${amp}/02_DADA2/table_DADA2_${amp}.qza \
    --o-representative-sequences ../demux_fwd/${amp}/02_DADA2/rep-seqs_DADA2_${amp}.qza \
    --o-denoising-stats ../demux_fwd/${amp}/02_DADA2/denoisingstats_DADA2_${amp}.qza \
    --verbose


    qiime metadata tabulate \
    --m-input-file ../demux_fwd/${amp}/02_DADA2/denoisingstats_DADA2_${amp}.qza \
    --o-visualization ../demux_fwd/${amp}/02_DADA2/denoisingstats_DADA2_vis_${amp}.qzv


    qiime feature-table summarize \
    --i-table ../demux_fwd/${amp}/02_DADA2/table_DADA2_${amp}.qza \
    --o-visualization ../demux_fwd/${amp}/02_DADA2/table_${amp}_vis.qzv \
    --m-sample-metadata-file ../metadata.txt

    qiime feature-table tabulate-seqs \
    --i-data ../demux_fwd/${amp}/02_DADA2/rep-seqs_DADA2_${amp}.qza \
    --o-visualization ../demux_fwd/${amp}/02_DADA2/rep-seqs_DADA2_${amp}_vis.qzv
done
```
## Generate a tree

```
for amp in 161 164 181 rbcl coi 
do
    mkdir -p ../demux_fwd/${amp}/03_tree

    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences ../demux_fwd/${amp}/02_DADA2/rep-seqs_DADA2_${amp}.qza \
        --o-alignment ../demux_fwd/${amp}/03_tree/aligned-rep-seqs-${amp}.qza \
        --o-masked-alignment ../demux_fwd/${amp}/03_tree/masked-aligned-rep-seqs-${amp}.qza \
        --o-tree ../demux_fwd/${amp}/03_tree/unrooted-tree-${amp}.qza \
        --o-rooted-tree ../demux_fwd/${amp}/03_tree/rooted-tree-${amp}.qza
    

done
```
##Diversity metrics
Adjust sampling depth per amplicon based on 02_DADA2/table_DADA2_${amp}

```

for amp in 161 164 181 rbcl coi 
do
    
    
    
        
    if [ "${amp}" = "161" ] 
    then
        qiime diversity core-metrics-phylogenetic \
            --i-phylogeny ../demux_fwd/${amp}/03_tree/rooted-tree-${amp}.qza \
            --i-table ../demux_fwd/${amp}/02_DADA2/table_DADA2_${amp}.qza  \
            --p-sampling-depth 10250 \
            --p-n-jobs-or-threads 1 \
            --m-metadata-file ../metadata.txt \
            --output-dir ../demux_fwd/${amp}/04_core-metrics-results

    elif [ "${amp}" = "164" ] 
    then
        qiime diversity core-metrics-phylogenetic \
            --i-phylogeny ../demux_fwd/${amp}/03_tree/rooted-tree-${amp}.qza \
            --i-table ../demux_fwd/${amp}/02_DADA2/table_DADA2_${amp}.qza  \
            --p-sampling-depth 10400 \
            --p-n-jobs-or-threads 1 \
            --m-metadata-file ../metadata.txt \
            --output-dir ../demux_fwd/${amp}/04_core-metrics-results
    elif [ "${amp}" = "181" ] 
    then
        qiime diversity core-metrics-phylogenetic \
            --i-phylogeny ../demux_fwd/${amp}/03_tree/rooted-tree-${amp}.qza \
            --i-table ../demux_fwd/${amp}/02_DADA2/table_DADA2_${amp}.qza  \
            --p-sampling-depth 9070 \
            --p-n-jobs-or-threads 1 \
            --m-metadata-file ../metadata.txt \
            --output-dir ../demux_fwd/${amp}/04_core-metrics-results
    elif [ "${amp}" = "rbcl" ]
    then
        qiime diversity core-metrics-phylogenetic \
            --i-phylogeny ../demux_fwd/${amp}/03_tree/rooted-tree-${amp}.qza \
            --i-table ../demux_fwd/${amp}/02_DADA2/table_DADA2_${amp}.qza  \
            --p-sampling-depth 4650 \
            --p-n-jobs-or-threads 1 \
            --m-metadata-file ../metadata.txt \
            --output-dir ../demux_fwd/${amp}/04_core-metrics-results
    elif [ "${amp}" = "coi" ] 
    then
        qiime diversity core-metrics-phylogenetic \
            --i-phylogeny ../demux_fwd/${amp}/03_tree/rooted-tree-${amp}.qza \
            --i-table ../demux_fwd/${amp}/02_DADA2/table_DADA2_${amp}.qza  \
            --p-sampling-depth 3580 \
            --p-n-jobs-or-threads 1 \
            --m-metadata-file ../metadata.txt \
            --output-dir ../demux_fwd/${amp}/04_core-metrics-results

    fi


done
```
##Taxonomy assignment
16S and 18S on silva138 database obtained from [QIIME2 data resources](https://docs.qiime2.org/2021.2/data-resources/) then trimmed to relevant primer regions
COI reads obtained  at [dereplicated stage](https://forum.qiime2.org/t/building-a-coi-database-from-bold-references/16129) [see also](https://doi.org/10.1002/ece3.6594)
rbcL classifier trained on [diat.database](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/TOMBYZ)

```
for amp in 161 164 181 rbcl coi 
do
    mkdir ../demux_fwd/${amp}/05_taxonomic_assignment
done
        
    qiime feature-classifier classify-sklearn \
    --i-classifier ../classifier/silva-16SV1-classifier.qza \
    --i-reads ../demux_fwd/161/02_DADA2/rep-seqs_DADA2_161.qza \
    --p-n-jobs 4 \
    --o-classification ../demux_fwd/161/05_taxonomic_assignment/taxonomy_161.qza


    qiime feature-classifier classify-sklearn \
    --i-classifier ../classifier/silva-16SV4-classifier.qza \
    --i-reads ../demux_fwd/164/02_DADA2/rep-seqs_DADA2_164.qza \
    --p-n-jobs 4 \
    --o-classification ../demux_fwd/164/05_taxonomic_assignment/taxonomy_164.qza

    qiime feature-classifier classify-sklearn \
    --i-classifier ../classifier/silva-18SV1-classifier.qza \
    --i-reads ../demux_fwd/181/02_DADA2/rep-seqs_DADA2_181.qza \
    --p-n-jobs 4 \
    --o-classification ../demux_fwd/181/05_taxonomic_assignment/taxonomy_181.qza

    
    qiime tools import \  
    --type 'FeatureData[Sequence]' \
    --input-path diat_seqs.fasta \
    --output-path diat_seqs.qza

    qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path diat_taxonomy.txt \
    --output-path diat_taxonomy.qza

    qiime feature-classifier extract-reads \
    --i-sequences diat_seqs.qza \
    --p-f-primer ATGCGTTGGAGAGARCGTTTC \
    --p-r-primer GATCACCTTCTAATTTACCWACAACTG \
    --p-min-length 150 \
    --p-max-length 450 \
    --o-reads diat_trim_seqs.qza

    qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads ../classifier/diat/rbcl-diat-ref-seqs.qza \
    --i-reference-taxonomy ../classifier/diat/diat_taxonomy_rbclonly.qza  \
    --o-classifier ../classifier/diat/diat-rbclonly_trimmed-classifier.qza
    qiime feature-classifier classify-sklearn \
    --i-classifier ../classifier/diat-rbclonly_trimmed-classifier.qza\
    --i-reads ../demux_fwd/rbcl/02_DADA2/rep-seqs_DADA2_rbcl.qza \
    --p-n-jobs 4 \
    --o-classification ../demux_fwd/rbcl/05_taxonomic_assignment/taxonomy_rbcl.qza

    qiime tools import \  
    --type 'FeatureData[Sequence]' \
    --input-path bold_derep1_seqs.fasta \
    --output-path bold_derep1_seqs.qza

    qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path coi_BOLD_taxonomy.txt \
    --output-path bold_derep1_taxa.qza

    qiime feature-classifier extract-reads \
    --i-sequences bold_derep1_seqs.qza \
    --p-f-primer  GGWACWGGWTGAACWGTWTAYCCYCC \
    --p-r-primer  TAAACTTCAGGGTGACCAAAAAATCA \
    --p-min-length 150 \
    --p-max-length 850 \
    --o-reads trim_bold_derep1_seqs.qza

    qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads ../classifier/coi_rescript_bold/trim_bold_derep1_seqs.qza \
    --i-reference-taxonomy ../classifier/coi_rescript_bold/bold_derep1_taxa.qza  \
    --o-classifier ../classifier/coi_trimmed_bold_classifier.qza

        
    qiime feature-classifier classify-sklearn \
    --i-classifier ../classifier/coi_trimmed_bold_classifier.qza \
    --i-reads ../demux_fwd/coi/02_DADA2/rep-seqs_DADA2_coi.qza \
    --p-n-jobs 1 \
    --o-classification ../demux_fwd/coi/05_taxonomic_assignment_trimmed/taxonomy_coi.qza
```

## Overall and pairwise PERMANOVA

```
for amp in 161 164 181 rbcl coi
do
        
    qiime diversity beta-group-significance \
        --i-distance-matrix ../demux_fwd/${amp}/04_core-metrics-results/weighted_unifrac_distance_matrix.qza \
        --m-metadata-file ../metadata.txt \
        --m-metadata-column phase_category \
        --o-visualization ../demux_fwd/${amp}/04_core-metrics-results/weighted-unifrac-beta-group-significance.qzv \
        --p-pairwise


done

for amp in 161 164 181 rbcl coi
do
   for phase1 in Recovery Pesticide Eutrophic
   do 
    for phase2 in Pesticide Eutrophic SemiPristine 
    do
    if [ "$phase1" = "$phase2" ] 
            then
                continue      # prevent duplication
            elif [ "$phase1" = "Eutrophic" ] && [ "$phase2" = "Pesticide" ] 
            then
                continue
            fi
            echo ${phase1}
            echo ${phase2}

            qiime diversity adonis \
                --i-distance-matrix ../demux_fwd/${amp}/04_core-metrics-results/weighted_unifrac_distance_matrix.qza \
                --m-metadata-file ../metadata.txt \
                --p-formula "phase_category" \
                --o-visualization ../demux_fwd/${amp}/04_core-metrics-results/weighted-unifrac-adonis.qzv 
        done
    done
done
```   
## Picrust

```
for amp in 161 164; do
    mkdir -p ${amp}/picrust_pipeline
    mkdir -p ${amp}/ko_phase
    mkdir -p ${amp}/ko_phase_pseud
    mkdir -p ${amp}/ancom_ko_phase
    mkdir -p ${amp}/ancom_ko_phase_out
    mkdir -p ${amp}/ancom_ko_data
    mkdir -p ${amp}/low_abund_ko_rm_phase/


    qiime picrust2 full-pipeline --i-table ${amp}/rarefied_table.qza 
    --i-seq ../${amp}/02_DADA2/rep-seqs_DADA2_${amp}.qza \
    --output-dir ${amp}/picrust_pipeline \
    --p-placement-tool epa-ng \
    --p-threads 3 \
    --p-hsp-method mp \
    --p-max-nsti 2 \
    --verbose

    for phase1 in Recovery Pesticide Eutrophic  ; do
        for phase2 in  Pesticide Eutrophic SemiPristine ; do
        
            if [ "$phase1" = "$phase2" ] 
            then
                continue     
            elif [ "$phase1" = "Eutrophic" ] && [ "$phase2" = "Pesticide" ] 
            then
                continue
            fi
            echo ${phase1}
            echo ${phase2}

            qiime feature-table filter-samples \
                --i-table ${amp}/picrust_pipeline/ko_metagenome.qza \
                --m-metadata-file ../../metadata.txt \
                --p-where "[phase_category] IN ('${phase1}', '${phase2}')" \
                --o-filtered-table ${amp}/ko_phase/ko_${phase1}_${phase2}.qza
            

            qiime feature-table filter-features \
                --i-table  ${amp}/ko_phase/ko_${phase1}_${phase2}.qza \
                --p-min-frequency 50 \
                --o-filtered-table ${amp}/low_abund_ko_rm_phase/ko_metagenome_min50.qza
            
            
            qiime composition add-pseudocount \
                --i-table ${amp}/ko_phase/ko_${phase1}_${phase2}.qza \
                --o-composition-table ${amp}/ko_phase_pseud/ko_pseud_${phase1}_${phase2}.qza
            
            qiime composition ancom \
                --i-table ${amp}/ko_phase_pseud/ko_pseud_${phase1}_${phase2}.qza \
                --m-metadata-file ../../metadata.txt  \
                --m-metadata-column phase_category \
                --o-visualization ${amp}/ancom_ko_phase/ancom_ko_${phase1}_${phase2}.qzv
            
            qiime tools export \
                --input-path ${amp}/ancom_ko_phase/ancom_ko_${phase1}_${phase2}.qzv \
                --output-path ${amp}/ancom_ko_phase_out/ancom_ko_${phase1}_${phase2} 
                
            cp ${amp}/ancom_ko_phase_out/ancom_ko_${phase1}_${phase2}/data.tsv ${amp}/ancom_ko_data/clr_${phase1}_${phase2}.tsv
            cp ${amp}/ancom_ko_phase_out/ancom_ko_${phase1}_${phase2}/ancom.tsv ${amp}/ancom_ko_data/ancom_${phase1}_${phase2}.tsv

            
        
         done
    done
done
```
