"""
PipelineExome.py - Tasks for variant calling
============================================

Reference
---------

"""
# Import modules
import os
import CGATCore.IOTools as IOTools
from CGATCore import Pipeline as P
import CGATCore.Experiment as E
from CGATCore.Pipeline import cluster_runnable
import numpy as np
import pandas as pd
import CGATCore.CSV as csv
import CGAT.VCF as VCF
import collections
import re
from urllib.request import urlopen
import itertools
from bs4 import BeautifulSoup, NavigableString
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as R
import copy

# Set PARAMS in calling module
#PARAMS = {}


def getGATKOptions():
    return "-l mem_free=1.4G"


def makeSoup(address):
    sock = urlopen(address)
    htmlSource = sock.read()
    soup = BeautifulSoup(htmlSource)
    return soup

##############################################################################


def GATKReadGroups(infile, outfile, genome,
                   library="unknown", platform="Illumina",
                   platform_unit="1", track="unknown",
                   tmp_dir='.',cluster_options = None):
    '''Reorders BAM according to reference fasta and adds read groups'''

    if track == 'unknown':
        track = P.snip(os.path.basename(infile), ".bam")
    tmpdir_gatk = P.get_temp_dir(tmp_dir)
    if cluster_options is None:
        cluster_options = {}
    #job_options = getGATKOptions()
    job_threads = 3

    statement = '''picard ReorderSam
                    INPUT=%(infile)s
                    OUTPUT=%(tmpdir_gatk)s/%(track)s.reordered.bam
                    REFERENCE=%(genome)s
                    ALLOW_INCOMPLETE_DICT_CONCORDANCE=true
                    VALIDATION_STRINGENCY=SILENT ;''' % locals()

    statement += '''samtools index %(tmpdir_gatk)s/%(track)s.reordered.bam ;
                 ''' % locals()

    statement += '''picard AddOrReplaceReadGroups
                    INPUT=%(tmpdir_gatk)s/%(track)s.reordered.bam
                    OUTPUT=%(outfile)s
                    RGLB=%(track)s
                    RGPL=%(platform)s
                    RGPU=%(platform_unit)s
                    RGSM=%(track)s
                    VALIDATION_STRINGENCY=SILENT ;''' % locals()

    statement += '''samtools index %(outfile)s ;
                 ''' % locals()
    statement += '''rm -rf %(tmpdir_gatk)s ;''' % locals()

    P.run(statement, **cluster_options)

##############################################################################


def GATKIndelRealign(infile, outfile, genome, intervals, padding, threads=4):
    '''Realigns BAMs around indels using GATK'''
    # Note that this is no longer recommanded if use HaplotypeCaller or MuTect2

    intervalfile = outfile.replace(".bam", ".intervals")
    job_options = getGATKOptions()
    job_threads = 3

    statement = '''GenomeAnalysisTK
                    -T RealignerTargetCreator
                    -o %(intervalfile)s
                    --num_threads %(threads)s
                    -R %(genome)s
                    -L %(intervals)s
                    -ip %(padding)s
                    -I %(infile)s; ''' % locals()

    statement += '''GenomeAnalysisTK
                    -T IndelRealigner
                    -o %(outfile)s
                    -R %(genome)s
                    -I %(infile)s
                    -targetIntervals %(intervalfile)s;''' % locals()
    P.run(statement)

##############################################################################


def GATKBaseRecal(infile, outfile, genome,
                intervals, padding, dbsnp, solid_options="",
                gatk_version=3, gatk=None, gatk_options='-Xmx4g',
                cluster_options=None):
    '''Recalibrates base quality scores using GATK'''

    track = P.snip(os.path.basename(infile), ".bam")
    tmpdir_gatk = P.get_temp_dir('.')
    job_options = getGATKOptions()
    job_threads = 3
    if cluster_options is None:
        cluster_options = {}

    if gatk_version == 3:
        statement = '''GenomeAnalysisTK
                        -T BaseRecalibrator
                        --out %(tmpdir_gatk)s/%(track)s.recal.grp
                        -R %(genome)s
                        -L %(intervals)s
                        -ip %(padding)s
                        -I %(infile)s
                        --knownSites %(dbsnp)s %(solid_options)s ;
                        ''' % locals()

        statement += '''GenomeAnalysisTK
                        -T PrintReads -o %(outfile)s
                        -BQSR %(tmpdir_gatk)s/%(track)s.recal.grp
                        -R %(genome)s
                        -I %(infile)s ;
                        ''' % locals()

        statement += '''rm -rf %(tmpdir_gatk)s ;''' % locals()
    elif gatk_version == 4:
        if gatk is None:
            errormsg = 'Currently you need to specify where to find GATK4'
            raise ValueError(errormsg)
        statement = '''{gatk} 
                        BaseRecalibrator
                        -I {infile}
                        -R {genome}
                        --known-sites {dbsnp} {solid_options}
                        -O {tmpdir_gatk}/{track}.recal.table ;
                        '''.format(**locals())
        statement += '''{gatk} 
                        ApplyBQSR
                        -R {genome}
                        -I {infile}
                        -bqsr {tmpdir_gatk}/{track}.recal.table
                        -O {outfile}
                        '''.format(**locals())
    else:
        errormsg = 'GATK version has to be either 3 or 4'
        raise ValueError(errormsg)
    P.run(statement, **cluster_options)

##############################################################################


def haplotypeCaller(infile, outfile, genome,
                    dbsnp, intervals, padding, options,
                    gatk_version=3, gatk=None, gatk_options='-Xmx4g',
                    cluster_options=None, wgs=False):
    '''Call SNVs and indels using GATK HaplotypeCaller in all members of a
    family together'''
    job_options = getGATKOptions()
    job_threads = 3
    if cluster_options is None:
        cluster_options = {}

    if gatk_version == 3:
        # no wgs here. not interested anyway
        statement = '''GenomeAnalysisTK
                        -T HaplotypeCaller
                        -ERC GVCF
                        -variant_index_type LINEAR
                        -variant_index_parameter 128000
                        -o %(outfile)s
                        -R %(genome)s
                        -I %(infile)s
                        --dbsnp %(dbsnp)s
                        -L %(intervals)s
                        -ip %(padding)s
                        %(options)s''' % locals()
    elif gatk_version == 4:
        statement = '''{gatk} --java-options {gatk_options}
                        HaplotypeCaller
                        -R {genome}
                        -I {infile}
                        -O {outfile}
                        -ERC GVCF
                        -L {intervals}
                        --dbsnp {dbsnp}'''.format(**locals())
        if wgs:
            statement += '''
                        --pcr-indel-model NONE
                        '''.format(**locals()) 
        else:
            statement += '''
                        -ip {paddings}
                        '''.format(**locals())
    P.run(statement, **cluster_options)


##############################################################################


def GenomicsDBImport(infiles, outfolder, chrom, gatk,
                    gatk_options='-Xmx4g', cluster_options=None):
    '''Move all GVCF files to GenomicsDB, one chromosome at a time'''
    # !!!Only work with GATK4!!!
    if cluster_options is None:
        cluster_options = {}
    statement = '{gatk} --java-options {gatk_options} GenomicsDBImport '.format(**locals())
    for infile in infiles:
        statement += '-V {} '.format(infile)
    statement += '''--genomicsdb-workspace-path {outfolder}
                    -L {chrom}
                    '''.format(**locals())
    P.run(statement, **cluster_options)                    


##############################################################################


def CombineGVCFs(inputfiles, outfile, genome, gatk, gatk_options='-Xmx4g'):
    '''Combine gvcfs for genotypeGVCFs (gatk4 only)'''
    # not tested!
    statement = '''{gatk} --java-options {gatk_options}
                CombineGVCFs
                -R {genome}
                -O {outfile}
                '''.format(**locals())
    for infile in inputfiles:
        statement += '-V {} '.format(infile)

    p.run(statement)



##############################################################################


def genotypeGVCFs(inputfiles, outfile, genome, gatk,
                gatk_version=3, gatk_options='-Xmx4g', From='',
                cluster_options = None):
    '''Joint genotyping of all samples together'''
    # Can get it from GVCFs or GenomicsDB
    # GenomicsDB only applies to gatk4
    # From has to be ('', 'gendb://')
    # If use gatk4 and From == '' (which is gvcf),
    #  the input has to be from CombineGVCFs
    # From = '' not tested for gatk4
    if cluster_options is None:
        cluster_options = {}

    job_options = getGATKOptions()
    job_threads = 4

    if gatk_version == 3:
        statement = '''GenomeAnalysisTK
                        -T GenotypeGVCFs
                        -o %(outfile)s
                        -R %(genome)s
                        --variant %(inputfiles)s''' % locals()
    elif gatk_version == 4:
        statement = '''{gatk} --java-options {gatk_options} 
                    GenotypeGVCFs
                    -V {From}{inputfiles}
                    -R {genome}
                    -O {outfile}
                    '''.format(**locals())
    P.run(statement, job_threads=job_threads, **cluster_options)


##############################################################################


def mutectSNPCaller(infile, outfile, mutect_log, genome, cosmic,
                    dbsnp, call_stats_out, job_memory, job_threads,
                    quality=20, max_alt_qual=150, max_alt=5,
                    max_fraction=0.05, tumor_LOD=6.3, strand_LOD=2,
                    normal_panel=None,
                    infile_matched=None,
                    gatk_key=None,
                    artifact=False):
    '''Call SNVs using Broad's muTect'''
    # TS. this is currently CGAT specific. How to generalise?

    job_memory_java = job_memory.lower()
    job_threads = 1

    statement = '''module load apps/java/jre1.6.0_26;
                   java -Xmx%(job_memory_java)s -jar
                   /ifs/apps/bio/muTect-1.1.4/muTect-1.1.4.jar
                   --analysis_type MuTect
                   --reference_sequence %(genome)s
                   --cosmic %(cosmic)s
                   --dbsnp %(dbsnp)s
                   --input_file:tumor %(infile)s
                   --out %(call_stats_out)s
                   --enable_extended_output
                   --vcf %(outfile)s
                ''' % locals()
    if artifact:
        statement += ''' --artifact_detection_mode '''

    if infile_matched:
        statement += '''--min_qscore %(quality)s
                        --gap_events_threshold 2
                        --max_alt_alleles_in_normal_qscore_sum %(max_alt_qual)s
                        --max_alt_alleles_in_normal_count %(max_alt)s
                        --max_alt_allele_in_normal_fraction %(max_fraction)s
                        --tumor_lod %(tumor_LOD)s
                        --input_file:normal %(infile_matched)s ''' % locals()
    if normal_panel:
        statement += ''' --normal_panel %(normal_panel)s ''' % locals()

    if gatk_key:
        statement += " -et NO_ET -K %(gatk_key)s " % locals()

    statement += " > %(mutect_log)s " % locals()

    P.run(statement)

##############################################################################


def strelkaINDELCaller(infile_control, infile_tumor, outfile, genome, config,
                       outdir, job_memory, job_threads):
    '''Call INDELs using Strelka'''

    statement = '''
    rm -rf %(outdir)s;
    /ifs/apps/bio/strelka-1.0.14/bin/configureStrelkaWorkflow.pl
    --normal=%(infile_control)s  --tumor=%(infile_tumor)s
    --ref=%(genome)s  --config=%(config)s  --output-dir=%(outdir)s;
    make -j %(job_threads)s -C %(outdir)s''' % locals()

    P.run(statement)

##############################################################################


def variantAnnotator(vcffile, bamlist, outfile, genome,
                     dbsnp, annotations="", snpeff_file=""):
    '''Annotate variant file using GATK VariantAnnotator'''
    job_options = getGATKOptions()
    job_threads = 3

    if "--useAllAnnotations" in annotations:
        anno = "--useAllAnnotations"
    elif annotations:
        anno = " -A " + " -A ".join(annotations)
    else:
        anno = ""

    statement = '''GenomeAnalysisTK -T VariantAnnotator
                    -R %(genome)s
                    -I %(bamlist)s
                    -A SnpEff
                    --snpEffFile %(snpeff_file)s
                    -o %(outfile)s
                    --variant %(vcffile)s
                    -L %(vcffile)s
                    --dbsnp %(dbsnp)s
                    %(anno)s''' % locals()
    P.run(statement)

##############################################################################


def variantRecalibrator(infile, outfile, genome, mode, dbsnp=None,
                        kgenomes=None, hapmap=None, omni=None, mills=None,
                        gatk_version=3, gatk=None, gatk_options='-Xmx4g',
                        wgs=False, cluster_options=None):
    '''Create variant recalibration file'''
    # If wgs == True, include DP as a train feature
    job_options = getGATKOptions()
    job_threads = 3
    if cluster_options is None:
        cluster_options = {}

    track = P.snip(outfile, ".recal")
    if mode == 'SNP' and gatk_version == 3:
        statement = '''GenomeAnalysisTK -T VariantRecalibrator
        -R %(genome)s
        -input %(infile)s
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %(hapmap)s
        -resource:omni,known=false,training=true,truth=true,prior=12.0 %(omni)s
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 %(dbsnp)s
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 %(kgenomes)s
        -an QD -an SOR -an MQRankSum
        -an ReadPosRankSum -an FS -an MQ
        --maxGaussians 4
        -mode %(mode)s
        -recalFile %(outfile)s
        -tranchesFile %(track)s.tranches
        -rscriptFile %(track)s.plots.R ''' % locals()
    elif mode == 'INDEL' and gatk_version == 3:
        statement = '''GenomeAnalysisTK -T VariantRecalibrator
        -R %(genome)s
        -input %(infile)s
        -resource:mills,known=true,training=true,truth=true,prior=12.0 %(mills)s
        -an QD -an MQRankSum
        -an ReadPosRankSum -an FS -an MQ
        --maxGaussians 4
        --minNumBadVariants 5000
        -mode %(mode)s
        -recalFile %(outfile)s
        -tranchesFile %(track)s.tranches
        -rscriptFile %(track)s.plots.R ''' % locals()
    elif mode == 'SNP' and gatk_version == 4:
        statement = '''{gatk} --java-options {gatk_options}
            VariantRecalibrator
            -R {genome}
            -V {infile}
            --resource hapmap,known=false,training=true,truth=true,prior=15.0:{hapmap}
            --resource omni,known=false,training=true,truth=true,prior=12.0:{omni}
            --resource 1000G,known=false,training=true,truth=false,prior=10.0:{kgenomes}
            --resource dbsnp,known=true,training=false,truth=false,prior=2.0:{dbsnp}
            -an QD 
            -an FS 
            -an SOR 
            -an MQ
            -an MQRankSum 
            -an ReadPosRankSum 
            -mode SNP
            --max-gaussians 4
            -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0
            -O {outfile}
            --tranches-file {track}.tranches
            --rscript-file {track}.plots.R
            '''.format(**locals())
    elif mode == 'INDEL' and gatk_version == 4:
        statement = '''{gatk} --java-options {gatk_options}
            VariantRecalibrator
            -R {genome}
            -V {infile}
            --resource mills,known=true,training=true,truth=true,prior=12.0:{mills}
            -an QD
            -an FS
            -an SOR
            -an MQ
            -an MQRankSum
            -an ReadPosRankSum
            --max-gaussians 4
            -mode INDEL
            -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0
            -O {outfile}
            --tranches-file {track}.tranches
            --rscript-file {track}.plots.R
            '''.format(**locals())
                
    if wgs:
        statement += ' -an DP '
    P.run(statement, **cluster_options)

##############################################################################


def applyVariantRecalibration(vcf, recal, tranches, outfile, genome, mode,
                            gatk_version=3, gatk=None, gatk_options='-Xmx4g',
                            cluster_options=None):
    '''Perform variant quality score recalibration using GATK '''
    job_options = getGATKOptions()
    job_threads = 3
    if cluster_options is None:
        cluster_options = {}

    if gatk_version == 3:
        statement = '''GenomeAnalysisTK -T ApplyRecalibration
        -R %(genome)s
        -input:VCF %(vcf)s
        -recalFile %(recal)s
        -tranchesFile %(tranches)s
        --ts_filter_level 99.0
        -mode %(mode)s
        -o %(outfile)s ''' % locals()
    elif gatk_version == 4:
        statement = '''{gatk} --java-options {gatk_options}
                    ApplyVQSR
                    -R {genome}
                    -V {vcf}
                    -ts-filter-level 99.0
                    -mode {mode}
                    -O {outfile}
                    --tranches-file {tranches}
                    --recal-file {recal}
                    '''.format(**locals())
                    
    P.run(statement, **cluster_options)

##############################################################################


def vcfToTable(infile, outfile, genome, columns):
    '''Converts vcf to tab-delimited file'''
    job_options = getGATKOptions()
    job_threads = 3

    statement = '''GenomeAnalysisTK -T VariantsToTable
                   -R %(genome)s
                   -V %(infile)s
                   --showFiltered
                   --allowMissingData
                   %(columns)s
                   -o %(outfile)s''' % locals()
    P.run(statement)

##############################################################################


def selectVariants(infile, outfile, genome, select):
    '''Filter de novo variants based on provided jexl expression'''
    statement = '''GenomeAnalysisTK -T SelectVariants
                    -R %(genome)s
                    --variant %(infile)s
                    -select '%(select)s'
                    -log %(outfile)s.log
                    -o %(outfile)s''' % locals()
    P.run(statement)

##############################################################################


def buildSelectStatementfromPed(filter_type, pedfile, template):
    '''Build a select statement from a template and a pedigree file'''
    pedigree = csv.DictReader(
        IOTools.open_file(pedfile), delimiter='\t', fieldnames=[
            'family', 'sample', 'father', 'mother', 'sex', 'status'])
    affecteds = []
    unaffecteds = []
    parents = []
    select = None
    # loop over pedigree file and establish relationships
    for row in pedigree:
        if row['status'] == '2':
            if filter_type == "denovo":
                father = row['father']
                mother = row['mother']
                proband = row['sample']
            elif filter_type == "dominant" or filter_type == "recessive":
                affecteds += [row['sample']]
            if filter_type == "recessive":
                parents += [row['father'], row['mother']]
        if row['status'] == '1':
            if filter_type == "dominant":
                unaffecteds += [row['sample']]
            elif filter_type == "recessive":
                if row['sample'] not in parents:
                    unaffecteds += [row['sample']]

    # Build select statement from template
    if filter_type == "denovo":
        select = template.replace("father", father)
        select = select.replace("mother", mother)
        select = select.replace("proband", proband)
    elif filter_type == "dominant":
        affecteds_exp = '").getPL().1==0&&vc.getGenotype("'.join(affecteds)
        if len(unaffecteds) == 0:
            unaffecteds_exp = ''
        else:
            unaffecteds_exp = '&&vc.getGenotype("' + \
                ('").isHomRef()&&vc.getGenotype("'.join(unaffecteds)) + \
                '").isHomRef()'
        select = template.replace("affecteds_exp", affecteds_exp)
        select = select.replace("unaffecteds_exp", unaffecteds_exp)
    elif filter_type == "recessive":
        affecteds_exp = '").getPL().2==0&&vc.getGenotype("'.join(affecteds)
        unaffecteds_exp = '").getPL().2!=0&&vc.getGenotype("'.join(unaffecteds)
        if len(parents) == 0:
            parents_exp = ''
        else:
            parents_exp = '&&vc.getGenotype("' + \
                ('").getPL().1==0&&vc.getGenotype("'.join(parents)) + \
                '").getPL().1==0'
        select = template.replace("affecteds_exp", affecteds_exp)
        select = select.replace("unaffecteds_exp", unaffecteds_exp)
        select = select.replace("parents_exp", parents_exp)

    return select

##############################################################################


def guessSex(infile, outfile):
    '''Guess the sex of a sample based on ratio of reads
    per megabase of sequence on X and Y'''
    statement = '''calc `samtools idxstats %(infile)s
                    | grep 'X'
                    | awk '{print $3/($2/1000000)}'`
                    /`samtools idxstats %(infile)s | grep 'Y'
                    | awk '{print $3/($2/1000000)}'`
                    | tr -d " " | tr "=" "\\t" | tr "/" "\\t"
                    > %(outfile)s'''
    P.run(statement)

##############################################################################


def filterMutect(infile, outfile, logfile,
                 control_id, tumour_id,
                 min_t_alt, min_n_depth,
                 max_n_alt_freq, min_t_alt_freq,
                 min_ratio):
    ''' filter the mutect snps'''

    reasons = collections.Counter()

    def comp(base):
        '''return complementary base'''
        comp_dict = {"C": "G", "G": "C", "A": "T", "T": "A"}
        return comp_dict[base]

    with IOTools.open_file(outfile, "w") as outf:
        with IOTools.open_file(infile, "r") as inf:
            for line in inf.readlines():
                # need to find location of control and tumor columns
                if line.startswith('#CHROM'):
                    columns = line.split("\t")
                    for x in range(0, len(columns)):
                        if control_id in columns[x]:
                            control_col = x
                        elif tumour_id in columns[x]:
                            tumor_col = x

                if line.startswith('#'):
                    # write out all comment lines
                    outf.write(line)

                else:
                    values = line.split("\t")

                    if values[6] == "PASS":
                        t_values = values[tumor_col].split(":")
                        t_ref, t_alt = list(
                            map(float, (t_values[2].split(","))))
                        t_depth = t_alt + t_ref
                        n_values = values[control_col].split(":")
                        n_ref, n_alt = list(
                            map(float, (n_values[2].split(","))))
                        n_depth = n_alt + n_ref
                        np.seterr(divide='ignore')

                        t_freq = np.divide(t_alt, t_depth)
                        n_freq = np.divide(n_alt, n_depth)

                        # filter
                        if not t_alt > min_t_alt:
                            reasons["Low_tumour_alt_count"] += 1
                            continue

                        if not t_freq >= min_t_alt_freq:
                            reasons["Low_tumour_alt_freq"] += 1
                            continue

                        if not n_depth >= min_n_depth:
                            reasons["Low_normal_depth"] += 1
                            continue

                        if not n_freq <= max_n_alt_freq:
                            reasons["high_normal_alt_freq"] += 1
                            continue

                        if (np.divide(t_freq, n_freq) >= min_ratio or n_freq == 0):
                            outf.write(line)
                    else:
                        reasons["Mutect_reject"] += 1

    with IOTools.open_file(logfile, "w") as outf:
        outf.write("%s\n" % "\t".join(("reason", "count")))
        for reason in reasons:
            outf.write("%s\t%i\n" % (reason, reasons[reason]))


# the following two functions should be generalised
# currently they operate only on mutect output
def compileMutationalSignature(infiles, outfiles):
    '''takes a list of mutect output files and compiles per sample mutation
    signatures'''

    delim = ":"

    def lookup(b1, b2):
        '''return lookup key for a pair of bases'''
        return(b1 + delim + b2)

    def breakKey(key):
        '''take a lookup key and return the elements'''
        return key.split(delim)

    def comp(base):
        '''return complementary base'''
        comp_dict = {"C": "G", "G": "C", "A": "T", "T": "A"}
        return comp_dict[base]

    def getID(infile):
        return P.snip(os.path.basename(infile),
                      ".mutect.snp.annotated.filtered.vcf")

    outfile1 = IOTools.open_file(outfiles[0], "w")
    mutations = ["C:T", "C:A", "C:G", "A:C", "A:T", "A:G"]

    outfile1.write("%s\t%s\t%s\t%s\t%s\n" % ("patient_id", "base_change",
                                             "ref", "alt", "frequency"))
    patient_freq = {}

    for infile in infiles:
        patient_id = getID(infile)
        mut_dict = {}
        for comb in mutations:
            mut_dict[comb] = 0

        with IOTools.open_file(infile, "r") as f:
            for line in f.readlines():
                if line.startswith('#'):
                    continue
                values = line.split("\t")
                key = lookup(values[3], values[4])
                if key in mut_dict:
                    mut_dict[key] += 1
                else:
                    comp_key = lookup(
                        comp(values[3]), comp(values[4]))
                    mut_dict[comp_key] += 1

        patient_freq[patient_id] = mut_dict

    for mutation in mutations:
        base1, base2 = breakKey(mutation)
        for infile in infiles:
            patient_id = getID(infile)
            outfile1.write("%s\t%s\t%s\t%s\t%s\n" % (patient_id, mutation,
                                                     base1, base2,
                                                     patient_freq[patient_id]
                                                     [mutation]))
    outfile1.close()

    outfile2 = IOTools.open_file(outfiles[1], "w")
    outfile2.write("%s\t%s\n" % ("patient_id",
                                 "\t".join(mutations)))
    for infile in infiles:
        patient_id = getID(infile)
        frequencies = "\t".join(map(str, [patient_freq[patient_id][x]
                                          for x in mutations]))
        outfile2.write("%s\t%s\n" % (patient_id, frequencies))
    outfile2.close()


##############################################################################
# utility functions for pipeline_exome_cancer
##############################################################################

@cluster_runnable
def parseMutectCallStats(infile, outfile):
    '''take the call stats outfile from mutect and summarise the
    reasons for variant rejection'''

    single_dict = collections.defaultdict(int)
    combinations_dict = collections.defaultdict(int)

    with IOTools.open_file(infile, "rb") as infile:
        lines = infile.readlines()
        for i, line in enumerate(lines):
            if i < 2:
                continue
            values = line.strip().split("\t")
            judgement, justification = (values[-1], values[-2])
            if judgement == "REJECT":
                reasons = justification.split(",")
                if len(reasons) == 1:
                    single_dict[reasons[0]] += 1
                else:
                    for reason in reasons:
                        combinations_dict[reasons[0]] += 1

    df = pd.DataFrame([single_dict, combinations_dict])

    df = df.transpose()
    df.columns = ["single", "combination"]
    df = df.sort("single", ascending=False)
    df.index.name = "justification"
    df.to_csv(outfile, header=True, index=True, sep="\t")


@cluster_runnable
def defineEBioStudies(cancer_types, outfile):
    ''' For the cancer types specified in pipeline.yml, identify the
    relevent studies in eBio '''

    cancer_types = cancer_types.split(",")

    cancer_studies_url = "http://www.cbioportal.org/webservice.do?cmd=getCancerStudies"
    genetic_profiles_url = "http://www.cbioportal.org/webservice.do?cmd=getGeneticProfiles"

    type2study_dict = collections.defaultdict(list)
    study2table_dict = collections.defaultdict(list)

    soup = makeSoup(cancer_studies_url)
    for lines in soup.body:
        if isinstance(lines, NavigableString):
            for line in lines.split("\n"):
                if "cancer_study_id" not in line:
                    values = str(line).strip().split("\t")
                    if len(values) > 1:
                        cancer_type = values[1].split(" (")[0]
                        if cancer_type in cancer_types:
                            study = re.sub(".\n", "", values[0])
                            type2study_dict[cancer_type].append(study)

    for study_list in list(type2study_dict.values()):
        for study in study_list:
            soup = makeSoup(genetic_profiles_url + "&cancer_study_id=" + study)
            lines = str(soup.body).split('\n')
            for line in lines:
                values = line.strip().split("\t")
                if len(values) > 1:
                    if values[1] == "Mutations":
                        genetic_profile = values[0]
                        study2table_dict[study] = genetic_profile

    outf = IOTools.open_file(outfile, "w")

    for cancer_type, study_id in type2study_dict.items():
        for study in study_id:
            table_id = study2table_dict[study]
            outf.write("%s\t%s\t%s\n" % (cancer_type, study, table_id))

    outf.close()


@cluster_runnable
def extractEBioinfo(eBio_ids, vcfs, outfile):
    '''find the number of mutations identitified in previous studies (eBio_ids)
    for the mutated genes in the vcfs'''

    genes = set()

    for vcf in vcfs:
        infile = VCF.VCFFile(IOTools.open_file(vcf))
        for vcf_entry in infile:
            # assumes all vcf entries without "REJECT" are "PASS"
            if vcf_entry.filter != "REJECT":
                info_entries = vcf_entry.info.split(";")
                for entry in info_entries:
                    if "SNPEFF_GENE_NAME" in entry:
                        genes.update((entry.split("=")[1],))

    eBio_ids = IOTools.open_file(eBio_ids, "r")

    tissue_counts = collections.defaultdict(
        lambda: collections.defaultdict(
            lambda: collections.defaultdict(int)))

    def chunks(l, n):
        ''' Yield successive n-sized chunks from l '''
        for i in range(0, len(l), n):
            yield l[i:i + n]

    # delete me
    E.info("number of genes: %i" % len(list(genes)))

    for line in eBio_ids:
        tissue, study, table = line.strip().split("\t")

        n = 0

        for i in range(0, len(list(genes)), 250):

            genes_chunk = list(genes)[i:i + 250]

            # TS sporadic error when querying with a single gene at a time
            # "urllib2.URLError: <urlopen error [Errno 110] Connection timed out>"
            # max URL length appears to be 8200 characters,
            # try doing 250 genes at a time?

            gene_list = "+".join(list(genes_chunk))

            n += len(genes_chunk)

            E.info("number of genes processed: %i" % n)

            url = ("http://www.cbioportal.org/webservice.do?cmd=getProfileData&"
                   "case_set_id=%(study)s_all&genetic_profile_id=%(table)s&"
                   "gene_list=%(gene_list)s" % locals())

            df = pd.io.parsers.read_csv(
                url, comment="#", sep="\t", index_col=0)

            for gene in genes_chunk:

                tmp_df = df[df['COMMON'] == gene]

                # check dataframe contains data!
                if tmp_df.shape[0] != 0:
                    # seem to be having issues with gene set containing duplicates!
                    # --> dataframe with repeated instances of gene after selection
                    # so splice to first row and recreate dataframe from series
                    if tmp_df.shape[0] > 1:
                        tmp_df = pd.DataFrame(tmp_df.iloc[0]).T

                    tissue_counts[tissue][gene]["total"] += tmp_df.shape[1] - 2
                    tissue_counts[tissue][gene][
                        "mutations"] += int(tmp_df.count(1)) - 1

    out = IOTools.open_file(outfile, "w")

    tissues = list(tissue_counts.keys())

    out.write("gene\t%s\n" % "\t".join([
        "%s_frequency" % x.replace(" ", "_") for x in tissues]))

    for gene in genes:
        freq_values = []
        for tissue in tissues:
            total = tissue_counts[tissue][gene]["total"]
            mutations = tissue_counts[tissue][gene]["mutations"]
            freq_values.append(round(np.divide(float(mutations), total), 4))

        out.write("%s\t%s\n" % (gene, "\t".join(map(str, freq_values))))

    out.close()


def intersectionHeatmap(infiles, outfile):
    ''' calculate the intersection between the infiles and plot'''

    pandas2ri.activate()

    name2genes = {}
    df = pd.DataFrame(columns=["id_1", "id_2", "intersection", "perc"])

    ix = 0
    for inf in infiles:

        name = P.snip(os.path.basename(inf)).split(".")[0]
        name = name.replace(".", "_")

        with IOTools.open_file(inf, "r") as f:
            genes = set()

            for line in f:
                if line[0] == "#":
                    continue

                values = line.strip().split("\t")
                info = values[7].split(";")

                for x in info:
                    if x.split("=")[0] == "SNPEFF_GENE_NAME":
                        gene_name = x.split("=")[1]
                        break

                # if no gene name found, line is skipped
                if gene_name:
                    genes.update((gene_name,))

        name2genes[name] = genes
        df.loc[ix] = [name, name, len(genes), 1.0]
        ix += 1

    for pair in itertools.permutations(list(name2genes.keys()), 2):
        id_1, id_2 = pair
        intersection = len(name2genes[id_1].intersection(name2genes[id_2]))
        not_intersecting = len(
            name2genes[id_1].symmetric_difference(name2genes[id_2]))
        intersection_perc = float(intersection) / (intersection +
                                                   not_intersecting)

        df.loc[ix] = [id_1, id_2, intersection, intersection_perc]
        ix += 1

    variant = os.path.basename(outfile).replace(
        "overlap_", "").replace("_heatmap.png", "")

    plotIntersectionHeatmap = R('''
    function(df){
    library(ggplot2)
    m_txt = element_text(size=15)
    m_txt_90 = element_text(size=15, angle=90, vjust=0.5, hjust=1)
    l_txt = element_text(size=20)

    p = ggplot(df, aes(id_1, id_2, fill=100*perc)) +
    geom_tile() +
    geom_text(aes(label=intersection), size=3) +
    scale_fill_gradient(name="Intersection (%%)", limits=c(0,100),
                       low="yellow", high="dodgerblue4") +
    theme(axis.text.x = m_txt_90, axis.text.y = m_txt,
          legend.text = m_txt, legend.title = m_txt,
          aspect.ratio=1) +
    xlab("") + ylab("") +
    ggtitle("%(variant)s")

    ggsave("%(outfile)s", width=10, height=10)
    }''' % locals())

    plotIntersectionHeatmap(df)


@cluster_runnable
def filterQuality(infile, qualstr, qualfilter, outfiles):
    '''
    Filter variants based on quality.  Columns to filter on and
    how they should be filtered can be specified in the pipeline.yml.
    Currently only implemented to filter numeric columns.  "." is assumed
    to mean pass.
    '''
    columns = IOTools.open_file(infile).readline()
    columns = columns.split("\t")
    qualparams = qualstr.split(",")
    qualdict = dict()
    fdict = dict()
    for param in qualparams:
        param = param.split("'")

        # column to filter on
        col = param[0]
        # string of >, <, >= or <= depending how the column should
        # be filtered
        lessmore = param[1]

        # score to filter by
        score = float(param[2])

        assert col in columns, "column %s not in variant table" % col

        ind = columns.index(col)
        i = 0
        iset = set([0, 1])
        with IOTools.open_file(infile) as input:
            for line in input:
                # rows one and two are headers
                if i > 1:
                    line = line.strip().split("\t")
                    if line[ind] == ".":
                        iset.add(i)
                    elif lessmore == ">":
                        if float(line[ind]) > score:
                            iset.add(i)
                    elif lessmore == ">=":
                        if float(line[ind]) > score:
                            iset.add(i)
                    elif lessmore == "<":
                        if float(line[ind]) < score:
                            iset.add(i)
                    elif lessmore == "<=":
                        if float(line[ind]) <= score:
                            iset.add(i)
                    if i not in iset:
                        fdict.setdefault(i, [])
                        fdict[i].append("%s=%s" % (col, line[ind]))
                i += 1
        qualdict[col] = iset
    if qualfilter == "all":
        allqual = set.intersection(*list(qualdict.values()))
    elif qualfilter == "any":
        allqual = set.union(*list(qualdict.values()))
    i = 0
    out = IOTools.open_file(outfiles[0], "w")
    out2 = IOTools.open_file(outfiles[1], "w")
    with IOTools.open_file(infile) as input:
        for line in input:
            if i in allqual:
                out.write(line)
            else:
                line = line.strip()
                out2.write("%s\t%s\n" % (line, ",".join(fdict[i])))
            i += 1
    out.close()
    out2.close()


@cluster_runnable
def FilterExacCols(infile, exac_suffs, exac_thresh):
    '''
    Returns a set of line indices indicating lines where either of the alleles
    called have a frequency of greater that exac_thresh in any of the
    populations specified as exac_suffs.
    Where no data is available an allele frequency of -1 is used.

    Exac provide data as AC_xxx and AN_xxx where AC is the allele count
    - the number of times the allele has been called
    - and AN is chromosome count - the number of
    samples in which the allele could have been called - in population xxx.
    AC / AN = allele frequecy.

    exac_suffs are any columns where an AC_xxx and AN_xxx column is provided
    in the VCF, e.g. Adj will calculate allele frequency from the AC_Adj
    and AN_Adj columns

    '''
    # read columns from the input VCF
    exac_suffs = exac_suffs.split(",")
    cols = IOTools.open_file(infile).readline().strip().split("\t")
    nD = dict()
    afdict = dict()
    for e in exac_suffs:
        # find the columns with the appropriate information
        # Allele count
        AC_i = cols.index("AC_%s" % (e))
        # Allele Number
        AN_i = cols.index("AN_%s" % (e))
        # Genotype
        GT_i = cols.index('GT')
        nlist = set()
        n = 0
        AFS = []
        with IOTools.open_file(infile) as input:
            for line in input:
                if n > 1:
                    line = line.strip().split("\t")
                    # At multi-allelic sites, comma delimited AC and AN values
                    # are provided
                    # "." and "NA" indicate no data here - this is represented
                    # as an AF of -1
                    AC = line[AC_i].replace(".", "-1").replace(
                        "NA", "-1").split(",")
                    AN = line[AN_i].replace(".", "1").replace(
                        "NA", "1").split(",")
                    AC = np.array([float(a) for a in AC])
                    AN = np.array([float(a) for a in AN])
                    AF = AC / AN
                    AF2 = [af if af > 0 else 0 for af in AF]
                    AF = np.insert(AF, 0, (1 - sum(AF2)))

                    # Chromosome count is usually the same for all minor
                    # alleles (but not always)
                    # If it is not the same the AC and AN lists should be the
                    # same length
                    # Otherwise AN will have length 1
                    if len(AC) != len(AN):
                        AN = [AN] * len(AC)

                    # Record the genotype called in this sample for this SNP
                    GT = line[GT_i]
                    GT = GT.replace(".", '0')
                    GT = GT.split("/")
                    GT[0], GT[1] = int(GT[0]), int(GT[1])

                    # If the variant is not in ExAC the ExAC columns show "."
                    # but the site
                    # may still have been called as multi allelic
                    # - use -1 for all frequencies
                    # in this case
                    if max(GT) > (len(AF) - 1):
                        AF = np.array([-1] * (max(GT) + 1))

                    AF1 = AF[GT[0]]
                    AF2 = AF[GT[1]]
                    AFS.append((AF1, AF2))
                    # Remember where both allele frequencies are
                    # greater than exac_thresh
                    if AF1 >= exac_thresh and AF2 >= exac_thresh:
                        nlist.add(n)
                else:
                    AFS.append(('NA', 'NA'))
                n += 1
        afdict[e] = AFS
        nD[e] = nlist

    ns = set.union(*list(nD.values()))
    return afdict, ns


@cluster_runnable
def FilterFreqCols(infile, thresh, fcols):
    '''
    Returns a set of line indices indicating lines where either of the alleles
    called have a frequency of less than thresh in all of the columns specified
    in fcols.
    No information - assigned allele frequency of -1.
    '''
    fcols = fcols.split(",")
    # read the column headings from the variant table
    cols = IOTools.open_file(infile).readline().strip().split("\t")
    # store allele frequency columns
    AFdict = dict()
    # store low frequency indices
    nD = dict()
    for col in fcols:
        ind = cols.index(col)
        GT_i = cols.index('GT')
        n = 0
        nlist = set()
        AFS = []
        with IOTools.open_file(infile) as input:
            for line in input:
                if n > 1:
                    line = line.strip().split("\t")
                    GT = line[GT_i].replace(".", "0").split("/")
                    af = line[ind].split(",")
                    AF = []
                    # where the allele frequency is not numeric
                    # "." or "NA" use -1 to indicate no data
                    for a in af:
                        try:
                            AF.append(float(a))
                        except:
                            AF.append(float(-1))
                    AF2 = [l if l > 0 else 0 for l in AF]
                    AF = np.array(AF)
                    AF = np.insert(AF, 0, 1 - sum(AF2))
                    GT[0] = int(GT[0])
                    GT[1] = int(GT[1])
                    # If the variant is not in database the column shows "."
                    # but the site
                    # may still have been called as multi allelic
                    # - use -1 for all frequencies
                    # in this case
                    if max(GT[0], GT[1]) > (len(AF) - 1):
                        AF = [float(-1)] * (max(GT[0], GT[1]) + 1)
                    AF1 = AF[GT[0]]
                    AF2 = AF[GT[1]]
                    if AF1 >= thresh and AF2 >= thresh:
                        nlist.add(n)
                    AFS.append((AF1, AF2))
                else:
                    AFS.append(('NA', 'NA'))
                n += 1
        AFdict[col] = AFS
        nD[col] = nlist

    ns = set.union(*list(nD.values()))
    return AFdict, ns


def WriteFreqFiltered(infile, exacdict, exacinds, otherdict, otherinds,
                      outfiles):
    '''
    Writes the output of the frequency filtering steps to file, including
    the new columns showing calculated allele frequency for the allele
    in this specific sample.
    '''
    x = 0
    out = IOTools.open_file(outfiles[0], "w")
    out2 = IOTools.open_file(outfiles[1], "w")

    exaccols = list(exacdict.keys())
    othercols = list(otherdict.keys())

    # column names for the new columns
    exacnewcols = ["%s_calc" % c for c in exaccols]
    othernewcols = ["%s_calc" % c for c in othercols]

    with IOTools.open_file(infile) as infile:
        for line in infile:
            line = line.strip()
            if x <= 1:
                # write the column names
                out.write("%s\t%s\t%s\n" % (line, "\t".join(exacnewcols),
                                            "\t".join(othernewcols)))
            else:
                freqs = []
                for key in exaccols:
                    col = exacdict[key]
                    freqs.append(col[x])
                for key in othercols:
                    col = otherdict[key]
                    freqs.append(col[x])
                freqs = [str(freq) for freq in freqs]

                if x not in exacinds and x not in otherinds:
                    out.write("%s\t%s\n" % (line, "\t".join(freqs)))
                else:
                    out2.write("%s\t%s\n" % (line, "\t".join(freqs)))
            x += 1
    out.close()
    out2.close()


@cluster_runnable
def filterRarity(infile, exac, freqs, thresh, outfiles):
    '''
    Filter out variants which are common in any of the exac or other
    population datasets as specified in the pipeline.yml.
    '''
    exacdict, exacinds = FilterExacCols(infile, exac, thresh)
    otherdict, otherinds = FilterFreqCols(infile, thresh, freqs)
    WriteFreqFiltered(infile, exacdict, exacinds, otherdict,
                      otherinds, outfiles)


@cluster_runnable
def filterDamage(infile, damagestr, outfiles):
    '''
    Filter variants which have not been assessed as damaging by any
    of the specified tools.
    Tools and thresholds can be specified in the pipeline.yml.

    Does not account for multiple alt alleles - if any ALT allele has
    been assessed as damaging with any tool the variant is kept,
    regardless of if this is the allele called in the sample.

    '''
    damaging = damagestr.split(",")
    cols = IOTools.open_file(infile).readline().strip().split("\t")

    D = dict()
    # parses the "damage string" from the pipeline.yml
    # this should be formatted as COLUMN|result1-result2-...,COLUMN|result1...
    # where variants with any of these results in this column will
    # be retained
    for d in damaging:
        d = d.split("|")
        col = d[0]
        res = d[1].split("-")
        i = cols.index(col)
        D[col] = ((res, i))

    x = 0
    out = IOTools.open_file(outfiles[0], "w")
    out2 = IOTools.open_file(outfiles[1], "w")
    with IOTools.open_file(infile) as input:
        for line in input:
            if x > 1:
                # grep for specific strings within this column of this
                # line of the input file
                line = line.strip().split("\t")
                isdamaging = 0
                for key in D:
                    res, i = D[key]
                    current = line[i]
                    for r in res:
                        if re.search(r, current):
                            isdamaging = 1
                if isdamaging == 1:
                    out.write("%s\n" % "\t".join(line))
                else:
                    out2.write("%s\n" % "\t".join(line))
            else:
                out.write(line)
            x += 1
    out.close()
    out2.close()


@cluster_runnable
def filterFamily(infile, infile2, outfiles):
    '''
    Filter variants according to the output of calculateFamily -
    only variants shared by both members of a family will be kept.
    '''
    cps1 = set()
    cps2 = set()

    # make a list of variants in infile1
    with IOTools.open_file(infile) as input:
        for line in input:
            line = line.strip().split("\t")
            chrom = line[0]
            pos = line[1]
            cp = "%s_%s" % (chrom, pos)
            cps1.add(cp)

    # make a list of variants in infile2
    with IOTools.open_file(infile2) as input:
        for line in input:
            line = line.strip().split("\t")
            chrom = line[0]
            pos = line[1]
            cp = "%s_%s" % (chrom, pos)
            cps2.add(cp)

    # only variants in both are of interest
    cps = cps1 & cps2

    out = IOTools.open_file(outfiles[0], "w")
    out2 = IOTools.open_file(outfiles[1], "w")
    with IOTools.open_file(infile) as input:
        for line in input:
            line = line.strip().split("\t")
            if "%s_%s" % (line[0], line[1]) in cps:
                out.write("%s\n" % "\t".join(line))
            else:
                out2.write("%s\n" % "\t".join(line))
    out.close()
    out2.close()


@cluster_runnable
def CleanVariantTables(genes, variants, cols, outfile):
    variants = pd.read_csv(variants, sep="\t")
    variants = variants.drop(0)

    vp1 = copy.copy(variants[['CHROM',
                              'POS', 'QUAL', 'ID', 'REF1', 'ALT', 'GT']])
    alleles = vp1['REF1'].str.cat(
        vp1['ALT'].str.strip(), sep=",").str.split(",")

    vp1['GT'] = vp1['GT'].str.replace(".", "0")
    inds1 = vp1['GT'].str.get(0).astype(int).values
    inds2 = vp1['GT'].str.get(-1).astype(int).values
    x = 0
    a1s = []
    a2s = []
    gts = []
    homhet = []
    for allele in alleles:
        i1 = int(inds1[x])
        i2 = int(inds2[x])
        a1 = allele[i1]
        a2 = allele[i2]
        a1s.append(a1)
        a2s.append(a2)
        if a1 == a2:
            homhet.append("HOM")
        else:
            homhet.append("HET")
        gts.append("%s%s" % (a1, a2))
        x += 1
    vp1['HOMHET'] = homhet
    vp1['Allele1'] = a1s
    vp1['Allele2'] = a2s
    vp1['Genotype'] = gts
    vp1 = vp1.drop(['REF1', 'ALT', 'GT'], 1)
    vp1[cols] = copy.copy(variants[cols])

    Ls = []
    for gene in [line.strip()
                 for line in IOTools.open_file(genes[0]).readlines()]:
        cp = []
        with IOTools.open_file(genes[1]) as infile:
            for line in infile:
                r = re.search(gene, line)
                if r:
                    line = line.strip().split("\t")
                    chrom = line[0]
                    pos = line[1]
                    cp.append("%s_%s" % (chrom, pos))
        cp = set(cp)
        for c in cp:
            Ls.append((gene, c.split("_")))
    df = pd.DataFrame(Ls)
    df['CHROM'] = df[1].str.get(0)
    df['POS'] = df[1].str.get(1)
    df = df.drop(1, 1)
    df.columns = ['gene', 'CHROM', 'POS']
    variants = vp1.merge(df, 'left')
    variants.to_csv(outfile, sep="\t")

###############################################################################

# population freq and cadd annotation
from CGAT import gnomad_utils, bravo_utils, kaviar_utils, utils

def pop_annotate(line_cache, variant_cache, header, fields, outf, PARAMS):
    # annotate
    # gnomad
    if 'gnomad_path' in PARAMS['pop_freqs']:
        gnomads = gnomad_utils.overall_freqs(
                list(variant_cache.keys()),
                PARAMS['pop_freqs']['gnomad_path']
        )
        for variant in variant_cache:
            af = gnomads[variant]['gnomad_af']
            if af is None:
                af = ''
            hom_f = gnomads[variant]['gnomad_hom_f']
            if hom_f is None:
                hom_f = ''
            variant_cache[variant]['gnomad_af'] = af
                     
            variant_cache[variant]['gnomad_hom_f'] = hom_f
                    
    # bravo
    if 'bravo_vcf' in PARAMS['pop_freqs']:
        bravos = bravo_utils.bravo(
                list(variant_cache.keys()),
                PARAMS['pop_freqs']['bravo_vcf']
        )
        for variant in variant_cache:
            if variant not in bravos:
                variant_cache[variant]['bravo_af'] = ''
                variant_cache[variant]['bravo_hom_f'] = ''
            else:
                variant_cache[variant]['bravo_af'] = \
                    bravos[variant]['af']
                variant_cache[variant]['bravo_hom_f'] = \
                    bravos[variant]['Hom']*2 / bravos[variant]['an']
    # kaviar
    if 'kaviar_vcf' in PARAMS['pop_freqs']:
        kaviars = kaviar_utils.kaviar(
                list(variant_cache.keys()),
                PARAMS['pop_freqs']['kaviar_vcf']
        )
        for variant in variant_cache:
            if variant not in kaviars:
                variant_cache[variant]['kaviar_af'] = ''
            else:
                variant_cache[variant]['kaviar_af'] = \
                    kaviars[variant]['af']

    for line in line_cache:
        record = dict(zip(header, line.rstrip().split('\t')))
        INFO = record['INFO']
        pop_info = []
        for alt in record['ALT'].split(','):
            v_id = '-'.join([
                record['CHROM'],
                record['POS'],
                record['REF'],
                alt,
            ])
            pop_info.append('|'.join(
                [alt] + [str(variant_cache[v_id][f]) for f in fields]
            ))
        pop_info = 'PF=' + ','.join(pop_info)
        new_INFO = ';'.join([INFO, pop_info])
        record['INFO'] = new_INFO
        new_line = '\t'.join([record[h] for h in header]) + '\n'
        outf.write(new_line)

def construct_pop_info(fields):
    info_line = (
            '##INFO=<ID=PF,Number=.,Type=String,Description="'
            'Population frequency such as gnomAD and Bravo. '
            'Format: Allele'
    )
    # add fields
    info_line = '|'.join([info_line] + fields)
    # closing
    info_line += '">\n'
    return info_line

def get_cadd_dict(cadd_file):
    result = {}
    with IOTools.open_file(cadd_file,'r') as inf:
        for line in inf:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            v_id = '-'.join(row[:4])
            result[v_id] = row[-1]
    return result

def cadd_annotate(line_cache, cadd_dict, header, outf):
    # annotate
    
    for line in line_cache:
        record = dict(zip(header, line.rstrip().split('\t')))
        INFO = record['INFO']
        cadd_info = []
        for alt in record['ALT'].split(','):
            if not set(alt).issubset('ATCGatcg'):
                continue
            v_id = '-'.join([
                record['CHROM'],
                record['POS'],
                record['REF'],
                alt,
            ])

            cleaned_id = utils.clean_variant(v_id)
            if cleaned_id not in cadd_dict:
                cleaned_id = utils.clean_variant(v_id, left=False)
            cadd_info.append(cadd_dict[cleaned_id])
        cadd_info = 'CADD=' + ','.join(cadd_info)
        new_INFO = ';'.join([INFO, cadd_info])
        record['INFO'] = new_INFO
        new_line = '\t'.join([record[h] for h in header]) + '\n'
        outf.write(new_line)

@cluster_runnable
def add_cadd(infile, outfile, cadd_file, PARAMS):
    '''
    add cadd given cadd.vcf.gz
    '''
    info_line = (
            '##INFO=<ID=CADD,Number=A,Type=Float,Description="'
            'CADD phred score">\n'
    )
    # no need to use tabix since cadd file is small
    cadd_dict = get_cadd_dict(cadd_file)
    line_cache = []
    variant_cache = {}
    header = []
    chrom = os.path.basename(infile).split('.')[1][3:]
    min_pos,max_pos = None,None
    with IOTools.open_file(infile, 'r') as inf, \
            IOTools.open_file(outfile, 'w') as outf:
        for line in inf:
            if line.startswith('##'):
                outf.write(line)
            elif line.startswith('#'):
                # add an INFO line
                outf.write(info_line)
                outf.write(line)
                header = line[1:].rstrip().split('\t')
            else:
                # write to cache
                if not set(line.split('\t')[3]).issubset(set('ATCG')):
                    continue
                if min_pos is None:
                    min_pos = int(line.split('\t')[1])
                line_cache.append(line)

                if len(line_cache) >= PARAMS['pop_freqs']['cache_size']:
                    # get cadd_dict
                    max_pos = int(line.split('\t')[1])
                    # dump cache
                    cadd_annotate(line_cache, cadd_dict, header, outf)
                    min_pos = max_pos
                    # empty caches
                    line_cache = []

        # finally dump cache
        if line_cache:
            max_pos = int(line.split('\t')[1])
            cadd_annotate(line_cache, cadd_dict, header, outf)
            line_cache = []

    
@cluster_runnable
def add_pop_freqs(infile, outfile, PARAMS):
    '''
    add pop freqs, such as gnomad, bravo and kaviar
    shove everything in INFO=<ID=PF.
    '''
    fields = []
    if 'gnomad_path' in PARAMS['pop_freqs']:
        fields.extend(['gnomad_af','gnomad_hom_f'])
    if 'bravo_vcf' in PARAMS['pop_freqs']:
        fields.extend(['bravo_af','bravo_hom_f'])
    if 'kaviar_vcf' in PARAMS['pop_freqs']:
        fields.append('kaviar_af')
    line_cache = []
    variant_cache = {}
    header = []
    with IOTools.open_file(infile, 'r') as inf, \
            IOTools.open_file(outfile, 'w') as outf:
        for line in inf:
            if line.startswith('##'):
                outf.write(line)
            elif line.startswith('#'):
                # add an INFO line
                info_line = construct_pop_info(fields)
                outf.write(info_line)
                outf.write(line)
                header = line[1:].rstrip().split('\t')
            else:
                # write to cache
                line_cache.append(line)
                row = line.split('\t')
                for alt in row[4].split(','):
                    v_id = '-'.join([row[0],row[1],row[3],alt])
                    variant_cache[v_id]={}

                if len(line_cache) >= PARAMS['pop_freqs']['cache_size']:
                    # dump cache
                    pop_annotate(line_cache, variant_cache, header, fields, outf, PARAMS)
                    # empty caches
                    line_cache = []
                    variant_cache = {}

        # finally dump cache
        if line_cache:
            pop_annotate(line_cache, variant_cache, header, fields, outf, PARAMS)
            line_cache = []
            variant_cache = {}
