process LIFTOVER_1 {
    label 'QC'
    scratch false
    publishDir "$params.OutDir/liftover",mode: "copy"

    input:
        val(clsa_plink)
        path(sqc)
        path(mqc)

    output:
        tuple path("clsa_gen_v3_sort_chr.bed"), path("clsa_gen_v3_sort_chr.bim"), path("clsa_gen_v3_sort_chr.fam"), emit : bfiles_sorted
    
    script:
    plink_file = files(clsa_plink)[0]
    plink_prefix = "${plink_file.getParent()}/${plink_file.getSimpleName()}" //path + prefix
    """
    awk 'NR > 1 {if(\$7 == 1){print \$1"\t"\$1}}' $sqc > remove_sqc.txt

    awk 'NR > 1 {if(\$8 == 1 || \$9 == 1 || \$10 == 1 || \$11 == 1 || \$12 == 1 || \$13 == 1){print \$1"\t"\$1}}' clsa_mqc_v3.txt > exclude_mqc.txt

    plink --bfile $plink_prefix --remove remove_sqc.txt --exclude exclude_mqc.txt --make-bed --out clsa_gen_v3_filtered

    plink --bfile clsa_gen_v3_filtered --make-bed --out clsa_gen_v3_sort
    plink --bfile clsa_gen_v3_sort --merge-x --make-bed --out clsa_gen_v3_sort_merged
    plink --bfile clsa_gen_v3_sort_merged --make-bed --output-chr chrMT --out clsa_gen_v3_sort_chr
    """
}

process BIM_TO_BED {
    label 'python_script'
    executor 'local'
    publishDir "$params.OutDir/liftover",mode: "copy"

    input:
        tuple path(bed), path(bim), path(fam)

    output:
        path("liftover_redo_UCSC_BED")

    script:
    """
#!/usr/bin/python
with open("${bim}", 'r') as bim, open("liftover_redo_UCSC_BED", 'w') as bed:
    for line in bim:
        fields = line.strip().split()
        chrom = fields[0]
        chrom_start = int(fields[3]) - 1
        chrom_end = int(fields[3])
        rsid = fields[1]
        
        bed_fields = [str(chrom), str(chrom_start), str(chrom_end), str(rsid)]
        bed.write('\\t'.join(bed_fields) + '\\n')
    """
}

process LIFTOVER_2 {
    label 'QC'
    scratch false
    publishDir "$params.OutDir/liftover",mode: "copy"

    input:
        val(ucsc_bed)
        path(chain)
        tuple path(bed), path(bim), path(fam)

    output:
        tuple path("clsa_gen_v3_Hg38.bed"), path("clsa_gen_v3_Hg38.bim"), path("clsa_gen_v3_Hg38.fam")
    
    script:
    bfiles_prefix = bed.getSimpleName() // prefix
    """
    liftOver $ucsc_bed $chain lifted.bed lifted.unmapped
    plink --bfile $bfiles_prefix --update-map lifted.bed 3 4 --exclude lifted.unmapped --make-bed --keep-allele-order --output-chr chrMT --out clsa_gen_v3_Hg38
    """
}

process PLINK2REFERENCE {
    label 'python_script'
    executor 'local'
    publishDir "$params.OutDir/liftover",mode: "copy"

    input:
        tuple path(bed), path(bim), path(fam)
        path(fasta)

    output:
        path("clsa_gen_v3_Hg38_ra.remove.txt"), emit : remove
        path("clsa_gen_v3_Hg38_ra.strand_flip.txt"), emit : strand_flip
        path("clsa_gen_v3_Hg38_ra.force_a1.txt"), emit: force_a1

    script:
    """
#!/usr/bin/env python3

import pysam

def strand_flip(a):
    return { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A' }[a] 

n_total_variants = 0
n_non_snps = 0
n_palindromic = 0
n_flip_strand = 0
n_force_ref_allele = 0
n_no_ref_match = 0 

with open('$bim', 'rt') as ibim, pysam.FastaFile('$fasta') as ifasta, \
    open('clsa_gen_v3_Hg38_ra.remove.txt', 'wt') as oremove, open('clsa_gen_v3_Hg38_ra.strand_flip.txt', 'wt') as oflip, open('clsa_gen_v3_Hg38_ra.force_a1.txt', 'wt') as oforce:

    fasta_chroms = set(list(ifasta.references))
    for line in ibim:
        fields = line.rstrip().split()
        chrom, varid, pos, a1, a2 = fields[0], fields[1], int(fields[3]), fields[4], fields[5]
        
        n_total_variants += 1

        if chrom not in fasta_chroms:
            chrom = chrom[3:] if chrom.startswith('chr') else f'chr{chrom}'
            if chrom not in fasta_chroms:
                print(f'Warning: skipping chromosome {fields[0]} because it is not in FASTA file.')
                continue

        if a1 not in {'A', 'C', 'G', 'T'} or a2 not in {'A', 'C', 'G', 'T'}:
            oremove.write(f'{varid}\\n')
            n_non_snps += 1
            continue

        if (a1 in {'A', 'T'} and a2 in {'A', 'T'}) or (a1 in {'C', 'G'} and a2 in {'C', 'G'}):
            oremove.write(f'{varid}\\n')
            n_palindromic += 1
            continue

        ref_base = None
        for base in ifasta.fetch(chrom, pos - 1, pos):
            ref_base = base

        if ref_base == a2:
            n_force_ref_allele += 1
            oforce.write(f'{varid}\\t{ref_base}\\n')
        elif ref_base != a1:
            flipped_a1 = strand_flip(a1)
            flipped_a2 = strand_flip(a2)
            if ref_base == flipped_a2:
                n_force_ref_allele += 1
                oforce.write(f'{varid}\\t{ref_base}\\n')
                n_flip_strand += 1
                oflip.write(f'{varid}\\n')
            elif ref_base == flipped_a1:
                n_flip_strand += 1
                oflip.write(f'{varid}\\n')
            else:
                n_no_ref_match += 1
                oremove.write(f'{varid}\\n')
    """
}

process LIFTOVER_3 {
    label 'QC'
    scratch false
    publishDir "$params.OutDir/liftover",mode: "copy"

    input:
        path(remove)
        path(strand_flip)
        path(force_a1)
        tuple path(bed), path(bim), path(fam)

    output:
        tuple path("clsa_gen_v3_Hg38_ra.bed"), path("clsa_gen_v3_Hg38_ra.bim"), path("clsa_gen_v3_Hg38_ra.fam")
    
    script: 
    plink_prefix = bed.getSimpleName() //prefix
    """
    plink --bfile $plink_prefix --exclude $remove --make-bed --out temp1
    plink --bfile temp1 --flip $strand_flip --make-bed --out temp2
    plink --bfile temp2 --a1-allele $force_a1 --make-bed --out clsa_gen_v3_Hg38_ra
    """
}

process GENOTYPE_QC_1 {
    label 'QC'
    scratch false
    publishDir "$params.OutDir/genotype_qc",mode: "copy"

    input:
        tuple path("clsa_gen_v3_Hg38_ra.bed"), path("clsa_gen_v3_Hg38_ra.bim"), path("clsa_gen_v3_Hg38_ra.fam")
        val(anderson_bed_gz)
        val(dust_masker_gz)
        val(sqc)

    output:
        path("R_check.het"), emit : r_check
        tuple path("clsa_gen_v3_1_23_Hg38_6.bed"), path("clsa_gen_v3_1_23_Hg38_6.bim"), path("clsa_gen_v3_1_23_Hg38_6.fam"), emit : plink_files
        path("non_euro.txt"), emit: non_euro

    script:
    """
    plink --bfile clsa_gen_v3_Hg38_ra --chr 1-22, X --make-bed --allow-extra-chr --out clsa_gen_v3_1_23_Hg38_x
    awk '{if (NR > 1) {if (\$5 != 4) print \$1,\$1}}' $sqc > non_euro.txt
    
    plink --bfile clsa_gen_v3_1_23_Hg38_x --remove non_euro.txt --make-bed -out clsa_gen_v3_1_23_Hg38_euro
    
    gzip -d $anderson_bed_gz | awk '{print \$0"\t"1}' > HG38.Anderson2010.set
    gzip -d $dust_masker_gz | awk '{print \$0"\t"1}' > HG38.dust_masker.set
    
    plink --bfile clsa_gen_v3_1_23_Hg38_euro --exclude 'range' HG38.Anderson2010.set --make-bed --out clsa_gen_v3_1_23_Hg38_euro_clean_1 
    plink --bfile clsa_gen_v3_1_23_Hg38_euro_clean_1 --exclude 'range' HG38.dust_masker.set --make-bed --out clsa_gen_v3_1_23_Hg38_euro_clean  

    plink --bfile clsa_gen_v3_1_23_Hg38_euro_clean --geno 0.1 --make-bed --out clsa_gen_v3_1_23_Hg38_1
    plink --bfile clsa_gen_v3_1_23_Hg38_1 --mind 0.1 --make-bed --out clsa_gen_v3_1_23_Hg38_2
    plink --bfile clsa_gen_v3_1_23_Hg38_2 --geno 0.01 --make-bed --out clsa_gen_v3_1_23_Hg38_3
    plink --bfile clsa_gen_v3_1_23_Hg38_3 --mind 0.01 --make-bed --out clsa_gen_v3_1_23_Hg38_4
    
    plink --bfile clsa_gen_v3_1_23_Hg38_4 --maf 0.05 --make-bed --out clsa_gen_v3_1_23_Hg38_5
    
    plink --bfile clsa_gen_v3_1_23_Hg38_5 --indep-pairwise 1000 100 0.9 --out indepSNP
    plink --bfile clsa_gen_v3_1_23_Hg38_5 --extract indepSNP.prune.in --make-bed --out clsa_gen_v3_1_23_Hg38_6
    
    plink --bfile clsa_gen_v3_1_23_Hg38_6 --het --out R_check
    """
}

process CALC_HET_OUTLIERS {
    label 'r_script'
    scratch false
    publishDir "$params.OutDir/genotype_qc",mode: "copy"

    input:
        val(r_check)

    output:
        path("fail-het-qc.txt")

    shell:
    '''
    #!/usr/bin/env Rscript

    het <- read.table("!{r_check}", head=TRUE)
    het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
    het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-5*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+5*sd(het$HET_RATE)));
    het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
    write.table(het_fail, "fail-het-qc.txt", row.names=FALSE)
    '''
}

process GENOTYPE_QC_2 {
    label 'QC'
    scratch false
    publishDir "$params.OutDir/genotype_qc",mode: "copy"

    input:
        tuple path("clsa_gen_v3_1_23_Hg38_6.bed"), path("clsa_gen_v3_1_23_Hg38_6.bim"), path("clsa_gen_v3_1_23_Hg38_6.fam")
        val(fail_het)
        val(mqc)
        val(sqc)

    output:
        tuple path("clsa_gen_v3_1_23_Hg38_10.bed"), path("clsa_gen_v3_1_23_Hg38_10.bim"), path("clsa_gen_v3_1_23_Hg38_10.fam"), emit : plink_files
        path("plink.sexcheck"), emit: plink_sexcheck
    
    script:
    """
    sed 's/"// g' $fail_het | awk '{print\$1, \$2}'> het_fail_ind.txt
    plink --bfile clsa_gen_v3_1_23_Hg38_6 --remove het_fail_ind.txt --make-bed --out clsa_gen_v3_1_23_Hg38_7

    plink --bfile clsa_gen_v3_1_23_Hg38_7 --check-sex 0.4 0.8 
    grep "PROBLEM" plink.sexcheck | awk '{print\$1,\$2}'> sex_discrepancy.txt
    plink --bfile clsa_gen_v3_1_23_Hg38_7 --remove sex_discrepancy.txt --make-bed --out clsa_gen_v3_1_23_Hg38_8 

    plink --bfile clsa_gen_v3_1_23_Hg38_8 --hwe 1e-15 --make-bed --out clsa_gen_v3_1_23_Hg38_9

    awk '{if (\$9 == 1){print \$1}}' $mqc > no_pass_hwe.txt
    
    plink --bfile clsa_gen_v3_1_23_Hg38_9 --exclude no_pass_hwe.txt --make-bed --out clsa_gen_v3_1_23_Hg38_10
    awk '{if (\$7 == 1){print \$1 "\t" \$1}}' $sqc > no_pass_het.txt

    plink --bfile clsa_gen_v3_1_23_Hg38_9 --remove no_pass_het.txt --make-bed --out clsa_gen_v3_1_23_Hg38_10
    """
}

process QUANT_PHENOTYPE_QC {
    label 'r_script'
    scratch false
    publishDir "$params.OutDir/quant_phenotype_qc",mode: "copy"

    input:
        tuple val(bed), val(bim), val(fam)
        val(sexcheck)
        val(sqc)
        val(base_pheno)

    output:
        path("phenotypes_bin_10.txt"), emit: phenotypes
        path("covariates.txt"), emit: covariates
        path("covariates_ss.txt"), emit: covariates_ss

    shell:
    '''
#!/usr/bin/env Rscript

library(dplyr)

CoP6_Baseline <- read.csv("/home/justb11/projects/def-gsarah/clsa/2104010_UdM_SGagliano_Baseline/2104010_UdM_SGagliano_CoP6_Baseline.csv", header = T, row.names = 1)
CoP6_Baseline$ADM_GWAS3_COM <- as.character(CoP6_Baseline$ADM_GWAS3_COM)

fam <- read.table("/scratch/justb11/Full_Pipeline/work/94/5adc6a0ffcd36e7075d306f96012ea/clsa_gen_v3_1_23_Hg38_10.fam", header = F)

CoP6_Edited <- CoP6_Baseline[!is.na(match(CoP6_Baseline$ADM_GWAS3_COM, fam$V2)),]

sqc <- read.table("/home/justb11/projects/def-gsarah/clsa/clsa_sqc_v3.txt", header = T)
sqc_subset <- sqc[,c("ADM_GWAS_COM", "batch", "ePC1", "ePC2", "ePC3", "ePC4", "ePC5", "ePC6", "ePC7", "ePC8", "ePC9", "ePC10")]
names(sqc_subset)[1] <- "ADM_GWAS3_COM"

CoP6_Cov <- CoP6_Edited[,c("ADM_GWAS3_COM", "ADM_DCS_AGE_COM")]

sexcheck <- read.table("/scratch/justb11/Full_Pipeline/work/94/5adc6a0ffcd36e7075d306f96012ea/plink.sexcheck", header = T)
sexcheck_ok <- sexcheck[sexcheck$STATUS == "OK",]
CoP6_Cov$BSEX <- sexcheck_ok$SNPSEX[match(CoP6_Cov$ADM_GWAS3_COM, sexcheck_ok$IID)]

CoP6_Cov_Full <- merge(CoP6_Cov, sqc_subset, by = c("ADM_GWAS3_COM"))

CoP6_Pheno <- CoP6_Edited[,c("ADM_GWAS3_COM","HWT_DBMI_COM", "VA_ETDRS_BOTH_RSLT_COM", "HGT_HEIGHT_M_COM", "BP_SYSTOLIC_ALL_AVG_COM", "GS_TRIAL1_MAX_COM",
                "FAS_F_SCORE_COM", "CRT_MRT_CORRANS_COM", "BP_PULSE_ALL_AVG_COM", "BLD_CREAT_COM", "TUG_TIME_COM", "TON_IOPG_R_COM",
                "WHC_RATIO_COM" , "DXA_WB_T_S_BMD_COM", "DXA_OI_TOTAL_PERCENT_FAT_COM")]

CoP6_Pheno$HWT_DBMI_COM[CoP6_Pheno$HWT_DBMI_COM %in% c(999.96, 999.99)] <- NA #no values are actually this..
CoP6_Pheno$VA_ETDRS_BOTH_RSLT_COM[CoP6_Pheno$VA_ETDRS_BOTH_RSLT_COM %in% c(-88.8)] <- NA
CoP6_Pheno$FAS_F_SCORE_COM[CoP6_Pheno$FAS_F_SCORE_COM %in% c(701)] <- NA #no values are actually this..
CoP6_Pheno$BLD_CREAT_COM[CoP6_Pheno$BLD_CREAT_COM %in% c(-8888)] <- NA
CoP6_Pheno$TUG_TIME_COM[CoP6_Pheno$TUG_TIME_COM %in% c(-88)] <- NA

CoP6_SDs <- sapply(CoP6_Pheno, function(x) sd(x, na.rm = T))
CoP6_means <- sapply(CoP6_Pheno, function(x) mean(x, na.rm = T))

CoP6_NoOutliers <- CoP6_Pheno

for (i in 2:ncol(CoP6_NoOutliers)){
    CoP6_NoOutliers[[i]][CoP6_NoOutliers[[i]] > (CoP6_means[[i]] + (CoP6_SDs[[i]] * 4)) | CoP6_NoOutliers[[i]] < (CoP6_means[[i]] - (CoP6_SDs[[i]] * 4))] <- NA
}

FID <- fam$V1[match(CoP6_Cov_Full$ADM_GWAS3_COM, fam$V2)]
CoP6_Cov_Full <- cbind(FID, CoP6_Cov_Full)
CoP6_Cov_Full <- CoP6_Cov_Full[,-3]

FID <- fam$V1[match(CoP6_NoOutliers$ADM_GWAS3_COM, fam$V2)]
CoP6_NoOutliers <- cbind(FID, CoP6_NoOutliers)

names(CoP6_Cov_Full)[2] <- "IID"
names(CoP6_NoOutliers)[2] <- "IID"

CoP6_Cov_ss <- CoP6_Cov_Full[,!(names(CoP6_Cov_Full) == "BSEX")]

write.table(CoP6_NoOutliers, "phenotypes_bin_10.txt", row.names=F, sep="	", quote = F)
write.table(CoP6_Cov_Full, "covariates.txt", row.names=F, sep="	", quote = F)
write.table(CoP6_Cov_ss, "covariates_ss.txt", row.names=F, sep="	", quote = F)
    '''
}

process SPLIT_PHENO{
    label 'QC'
    scratch false
    publishDir "$params.OutDir/quant_phenotype_qc",mode: "copy"

    input:
        path(phenos)

    output:
        path("phenotype_*.txt")

    script:
    """ 
    num_cols=\$(awk '{print NF; exit}' $phenos)

    for (( i=3; i<=\$num_cols; i++ ))
    do
        awk -v i=\${i} '{print \$1 "\t" \$2 "\t" \$i}' $phenos > phenotype_\${i}.txt
    done
    """
}

process REGENIE_STEP1 {
    label 'regenie_step1'
    scratch false
    publishDir "$params.OutDir/regenie_step1", mode: "copy"

    input:
        tuple val(bed), val(bim), val(fam)
        val(phenotypes)
        val(covariates)

    output:
        path("fit_bin_out_*.loco")
        path("fit_bin_out_pred.list"), emit: fit_pred

    script:
    plink_prefix = "${bed.getParent()}/${bed.getSimpleName()}" //path + prefix
    """

    regenie \
    --step 1 \
    --bed $plink_prefix \
    --covarFile $covariates \
    --phenoFile $phenotypes \
    --bsize 100 \
    --threads 8 \
    --out fit_bin_out
    """
}

process REGENIE_STEP1_MALE {
    label 'regenie_step1'
    scratch false
    publishDir "$params.OutDir/regenie_step1_male", mode: "copy"

    input:
        tuple val(bed), val(bim), val(fam)
        val(phenotypes)
        val(covariates)

    output:
        path("fit_bin_out_*.loco")
        path("fit_bin_out_pred.list"), emit: fit_pred

    script:
    plink_prefix = "${bed.getParent()}/${bed.getSimpleName()}" //path + prefix
    """
    unique_values=\$(cut -f3 $phenotypes | sort | uniq | wc -l)
    if ( \$unique_values -eq 2 || \$unique_values -eq 3 ) ; then
        binary="--bt --cc12 --force-qt"
    else
        binary=" "
    fi

    regenie \
    --step 1 \
    --bed $plink_prefix \
    --covarFile $covariates \
    --phenoFile $phenotypes \
    --sex-specific 'male' \
    --bsize 100 \
    --threads 8 \
    --out fit_bin_out
    """
}

process REGENIE_STEP1_FEMALE {
    label 'regenie_step1'
    scratch false
    publishDir "$params.OutDir/regenie_step1_female", mode: "copy"

    input:
        tuple val(bed), val(bim), val(fam)
        val(phenotypes)
        val(covariates)

    output:
        path("fit_bin_out_*.loco")
        path("fit_bin_out_pred.list"), emit: fit_pred

    script:
    plink_prefix = "${bed.getParent()}/${bed.getSimpleName()}" //path + prefix
    """
    regenie \
    --step 1 \
    --bed $plink_prefix \
    --covarFile $covariates \
    --phenoFile $phenotypes \
    --sex-specific 'female' \
    --bsize 100 \
    --threads 8 \
    --out fit_bin_out
    """
}

process SUBSET_BGEN {
    label 'PLINK2_Subset'
    scratch false
    publishDir "$params.OutDir/bgen", mode: "copy"

    input:
        each(bgen)
        path(sample)
        path(non_euro)

    output:
        path("clsa_imp_*_euro.pgen")

    script:
    chrom = (bgen.getBaseName() =~ /(?<=clsa_imp_)(.*?)(?=_v3)/)[0][0]
    """
    plink2 \
    --bgen $bgen ref-first \
    --sample $sample \
    --remove $non_euro \
    --make-pgen \
    --out clsa_imp_${chrom}_euro
    """
}

process ADD_SEX_TO_SAMPLE{
    label 'QC'
    scratch false
    publishDir "$params.OutDir/bgen", mode: "copy"

    input:
        path(sample_nosex)
        path(clsa_plink)

    output:
        val("clsa_imp_v3_withsex.sample")

    //take info from clsa_gen_v3.fam file
    script:
    fam = "${clsa_plink.getParent()}/clsa_gen_v3.fam"
    //remove sample file header.
    """
    head -n 1 $sample_nosex > clsa_imp_v3_withsex.sample
    join -j 1 -o 1.1,1.2,1.3,2.5 <(tail -n +2 $sample_nosex | sort -k1) <(sort -k1 $fam) >> clsa_imp_v3_withsex.sample
    """
}

process REGENIE_STEP2 {
    label 'regenie_step2'
    scratch false
    publishDir "$params.OutDir/regenie_step2", mode: "copy"

    input:
        each(pgen)
        each(phenotype_num)
        val(phenotypes)
        val(covariates)
        val(fit_pred)

    output:
        path("results*.regenie")

    script:
    phenotype = phenotypes[0]
    pgen_prefix = "${pgen.getParent()}/${pgen.getSimpleName()}" //path + prefix
    pgen_simple = pgen.getSimpleName()
    chrom_num = (pgen.getBaseName() =~ /(?<=clsa_imp_)(.*?)(?=_euro)/)[0][0]
    phenos_path = "${phenotype.getParent()}" //path
    """
    unique_values=\$(cut -f3 $phenotype | sort | uniq | wc -l)
    if ( \$unique_values -eq 2 || \$unique_values -eq 3 ) ; then
        binary="--bt --cc12"
    else
        binary=" "
    fi

    mkdir -p step_2_results_pheno_${phenotype_num}
    regenie --step 2 \
        --pgen $pgen_prefix \
        --covarFile $covariates \
        --phenoFile $phenos_path/phenotype_${phenotype_num}.txt \
        \$binary \
        --bsize 200 \
        --minMAC 20 \
        --ref-first \
        --threads 8 \
        --chr $chrom_num \
        --pred $fit_pred \
        --out results
    """
}

process REGENIE_STEP2_MALE {
    label 'regenie_step2'
    scratch false
    publishDir "$params.OutDir/regenie_step2_male", mode: "copy"

    input:
        each(pgen)
        each(phenotype_num)
        val(phenotypes)
        val(covariates)
        val(fit_pred)
        path(done)

    output:
        path("results*.regenie")

    script:
    phenotype = phenotypes[0]
    pgen_prefix = "${pgen.getParent()}/${pgen.getSimpleName()}" //path + prefix
    pgen_simple = pgen.getSimpleName()
    chrom_num = (pgen.getBaseName() =~ /(?<=clsa_imp_)(.*?)(?=_euro)/)[0][0]
    phenos_path = "${phenotype.getParent()}" //path
    """
    unique_values=\$(cut -f3 $phenotype | sort | uniq | wc -l)
    if ( \$unique_values -eq 2 || \$unique_values -eq 3 ) ; then
        binary="--bt --cc12"
    else
        binary=" "
    fi

    mkdir -p step_2_results_pheno_${phenotype_num}
    regenie --step 2 \
        --pgen $pgen_prefix \
        --covarFile $covariates \
        --phenoFile $phenos_path/phenotype_${phenotype_num}.txt \
        \$binary \
        --sex-specific 'male' \
        --bsize 200 \
        --minMAC 20 \
        --ref-first \
        --threads 8 \
        --chr $chrom_num \
        --pred $fit_pred \
        --out results
    """
}

process REGENIE_STEP2_FEMALE {
    label 'regenie_step2'
    scratch false
    publishDir "$params.OutDir/regenie_step2_female", mode: "copy"

    input:
        each(pgen)
        each(phenotype_num)
        val(phenotypes)
        val(covariates)
        val(fit_pred)
        path(done)

    output:
        path("results*.regenie")

    script:
    phenotype = phenotypes[0]
    pgen_prefix = "${pgen.getParent()}/${pgen.getSimpleName()}" //path + prefix
    pgen_simple = pgen.getSimpleName()
    chrom_num = (pgen.getBaseName() =~ /(?<=clsa_imp_)(.*?)(?=_euro)/)[0][0]
    phenos_path = "${phenotype.getParent()}" //path
    """
    unique_values=\$(cut -f3 $phenotype | sort | uniq | wc -l)
    if ( \$unique_values -eq 2 || \$unique_values -eq 3 ) ; then
        binary="--bt --cc12 --force-qt"
    else
        binary=" "
    fi
    mkdir -p step_2_results_pheno_${phenotype_num}
    regenie --step 2 \
        --pgen $pgen_prefix \
        --covarFile $covariates \
        --phenoFile $phenos_path/phenotype_${phenotype_num}.txt \
        \$binary \
        --sex-specific 'female' \
        --bsize 200 \
        --minMAC 20 \
        --ref-first \
        --threads 8 \
        --chr $chrom_num \
        --pred $fit_pred \
        --out results
    """
}

process CREATE_PHENO_LIST_BOTH {
    label 'QC'
    scratch false

    input:
        each(step2_results)

    output:
        path("pheno-list_unclosed.json")

    script:
    file_name = (step2_results.getName() =~ /(?<=results_)(.*?)(?=.regenie)/)[0][0]
    """
    echo '{"assoc_files" : ["$step2_results"], "phenocode" : "$file_name" }' > pheno-list_unclosed.json
    """
}

process CREATE_PHENO_LIST_MALE {
    label 'QC'
    scratch false

    input:
        each(step2_results)

    output:
        path("pheno-list_unclosed.json")

    script:
    file_name = (step2_results.getName() =~ /(?<=results_)(.*?)(?=.regenie)/)[0][0]
    """
    echo '{"assoc_files" : ["$step2_results"], "phenocode" : "$file_name" }' > pheno-list_unclosed.json
    """
}

process CREATE_PHENO_LIST_FEMALE {
    label 'QC'
    scratch false

    input:
        each(step2_results)

    output:
        path("pheno-list_unclosed.json")

    script:
    file_name = (step2_results.getName() =~ /(?<=results_)(.*?)(?=.regenie)/)[0][0]
    """
    echo '{"assoc_files" : ["$step2_results"], "phenocode" : "$file_name" }' > pheno-list_unclosed.json
    """
}

process CLOSE_JSON_BOTH {
    label 'QC'
    scratch false
    publishDir "$params.OutDir/pheweb_$sex",mode: "copy"

    input:
        path(pheno_list_unclosed)
        val(sex)

    output:
        path("pheno-list.json")

    script:
    """
    sed '1s/^/[/;\$!s/\$/,/;\$s/\$/]/' $pheno_list_unclosed > pheno-list.json
    """
}

process CLOSE_JSON_FEMALE {
    label 'QC'
    scratch false
    publishDir "$params.OutDir/pheweb_$sex",mode: "copy"

    input:
        path(pheno_list_unclosed)
        val(sex)

    output:
        path("pheno-list.json")

    script:
    """
    sed '1s/^/[/;\$!s/\$/,/;\$s/\$/]/' $pheno_list_unclosed > pheno-list.json
    """
}

process CLOSE_JSON_MALE {
    label 'QC'
    scratch false
    publishDir "$params.OutDir/pheweb_$sex",mode: "copy"

    input:
        path(pheno_list_unclosed)
        val(sex)

    output:
        path("pheno-list.json")

    script:
    """
    sed '1s/^/[/;\$!s/\$/,/;\$s/\$/]/' $pheno_list_unclosed > pheno-list.json
    """
}

process LOG10P_TO_P_SORT_BOTH {
    label 'r_script'
    scratch false
    publishDir "$params.OutDir/pheweb_$sex",mode: "copy"

    input:
        each(step2_results)
        val(sex)

    output:
        path("results*.regenie")

    shell:
    final_output = step2_results.getName()
    '''
#!/usr/bin/env Rscript
library(dplyr)

df <- read.table("!{step2_results}", header=TRUE)
df <- df %>%
    mutate(P = 10 ** -(LOG10P)) %>%
    arrange(CHROM, GENPOS)

write.table(df, file="!{final_output}", row.names=FALSE, sep="\t", quote = FALSE)
    '''
}

process LOG10P_TO_P_SORT_MALE {
    label 'r_script'
    scratch false
    publishDir "$params.OutDir/pheweb_$sex",mode: "copy"

    input:
        each(step2_results)
        val(sex)

    output:
        path("results*.regenie")

    shell:
    final_output = step2_results.getName()
    '''
#!/usr/bin/env Rscript
library(dplyr)

df <- read.table("!{step2_results}", header=TRUE)
df <- df %>%
    mutate(P = 10 ** -(LOG10P)) %>%
    arrange(CHROM, GENPOS)

write.table(df, file="!{final_output}", row.names=FALSE, sep="\t", quote = FALSE)
    '''
}

process LOG10P_TO_P_SORT_FEMALE {
    label 'r_script'
    scratch false
    publishDir "$params.OutDir/pheweb_$sex",mode: "copy"

    input:
        each(step2_results)
        val(sex)

    output:
        path("results*.regenie")

    shell:
    final_output = step2_results.getName()
    '''
#!/usr/bin/env Rscript
library(dplyr)

df <- read.table("!{step2_results}", header=TRUE)
df <- df %>%
    mutate(P = 10 ** -(LOG10P)) %>%
    arrange(CHROM, GENPOS)

write.table(df, file="!{final_output}", row.names=FALSE, sep="\t", quote = FALSE)
    '''
}

process CREATE_CONFIG_PY {
    label 'QC'
    scratch false
    publishDir "$params.OutDir/pheweb",mode: "copy"

    output:
        path("config.py")

    script:
    """
    echo hg_build_number = 38 > config.py
    echo -e "field_aliases = {'CHROM': 'chrom', 'GENPOS': 'pos', 'ALLELE0' : 'ref', 'ALLELE1': 'alt', 'BETA': 'beta', 'SE': 'sebeta', 'P': 'pval', 'N': 'num_samples'}" >> config.py
    """
}

workflow {

    //liftover
    bfiles_sorted = LIFTOVER_1(params.clsa_plink, params.mqc, params.sqc)
    ucsc_bed = BIM_TO_BED(bfiles_sorted)
    bfiles_hg38 = LIFTOVER_2(ucsc_bed, params.chain, bfiles_sorted)
    p2r_output = PLINK2REFERENCE(bfiles_hg38, params.fasta)
    lifted = LIFTOVER_3(p2r_output.remove, p2r_output.strand_flip, p2r_output.force_a1, bfiles_hg38)

    //QC scripts    
    genotype_QC_1 = GENOTYPE_QC_1(lifted, params.anderson_bed_gz, params.dust_masker_gz, params.mqc) //run the genotype QC
    fail_het = CALC_HET_OUTLIERS(genotype_QC_1.r_check)
    genotypeQCed = GENOTYPE_QC_2(genotype_QC_1.plink_files, fail_het, params.mqc, params.sqc) //run the genotype QC
    phenotypeQCed = QUANT_PHENOTYPE_QC(genotypeQCed.plink_files, genotypeQCed.plink_sexcheck, params.sqc, params.base_phenotypes) //create formatted phenotype and covariate file for regenie
    split_phenos = SPLIT_PHENO(phenotypeQCed.phenotypes)
    sample_file = ADD_SEX_TO_SAMPLE(params.clsa_samples, params.clsa_plink)

    //step 1 regenie
    step1 = REGENIE_STEP1(genotypeQCed.plink_files, phenotypeQCed.phenotypes, phenotypeQCed.covariates)
    step1_male = REGENIE_STEP1_MALE(genotypeQCed.plink_files, phenotypeQCed.phenotypes, phenotypeQCed.covariates_ss)
    step1_female = REGENIE_STEP1_FEMALE(genotypeQCed.plink_files, phenotypeQCed.phenotypes, phenotypeQCed.covariates_ss)

    //preperation for step 2
    bgens = Channel.fromPath(params.clsa_bgen)
    euro_pgens = SUBSET_BGEN(bgens, sample_file, genotype_QC_1.non_euro)

    //need to find a way to NOT hardcode this..
    phenotypes_cols = [3,4,5,6,7,8,9,10,11,12,13,14,15,16]

    //step 2 regenie
    //merge phenotypes into one file
    step2_results_both = REGENIE_STEP2(euro_pgens, phenotypes_cols, split_phenos, phenotypeQCed.covariates, step1.fit_pred).collectFile(storeDir: "$params.OutDir/regenie_step2/merged_both", keepHeader : true, sort : 'true')
    step2_results_male = REGENIE_STEP2_MALE(euro_pgens, phenotypes_cols, split_phenos, phenotypeQCed.covariates_ss, step1_male.fit_pred, step2_results_both).collectFile(storeDir: "$params.OutDir/regenie_step2/merged_male", keepHeader : true, sort : 'true')
    step2_results_female = REGENIE_STEP2_FEMALE(euro_pgens, phenotypes_cols, split_phenos, phenotypeQCed.covariates_ss, step1_female.fit_pred, step2_results_male).collectFile(storeDir: "$params.OutDir/regenie_step2/merged_male", keepHeader : true, sort : 'true')

    //transform logp into p using Rscript
    pheweb_ready_both = LOG10P_TO_P_SORT_BOTH(step2_results_both, 'both')
    pheweb_ready_male = LOG10P_TO_P_SORT_MALE(step2_results_male, 'male')
    pheweb_ready_female = LOG10P_TO_P_SORT_FEMALE(step2_results_female, 'female')

    //create pheno-list and config file for pheweb
    CREATE_CONFIG_PY()

    pheno_list_unclosed_both = CREATE_PHENO_LIST_BOTH(pheweb_ready_both).collectFile(storeDir: "$params.OutDir/pheweb/pheno-list_both")
    pheno_list_unclosed_male = CREATE_PHENO_LIST_MALE(pheweb_ready_male).collectFile(storeDir: "$params.OutDir/pheweb/pheno-list_male")
    pheno_list_unclosed_female = CREATE_PHENO_LIST_FEMALE(pheweb_ready_female).collectFile(storeDir: "$params.OutDir/pheweb/pheno-list_female")
    
    CLOSE_JSON_BOTH(pheno_list_unclosed_both, 'both')
    CLOSE_JSON_MALE(pheno_list_unclosed_male, 'male')
    CLOSE_JSON_FEMALE(pheno_list_unclosed_female, 'female')

}
