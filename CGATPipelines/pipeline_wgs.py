# load modules
from ruffus import transform, mkdir, follows, merge, regex, suffix, \
    jobs_limit, files, collate, add_inputs, formatter, \
    active_if, subdivide
from ruffus.combinatorics import permutations
import sys
import os
import csv
import glob
import copy
import re
import shutil
import decimal
import numpy as np
import pandas as pd
import CGATCore.Experiment as E
import CGATCore.IOTools as IOTools
from CGATCore import Pipeline as P
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGATPipelines.PipelineExome as PipelineExome
import CGATPipelines.PipelineExomeAncestry as PipelineExomeAncestry
from CGATPipelines.Report import run_report


# Use fastq/get_fastq.py to soft link 
#  all fastq files before use this script

PARAMS = P.get_parameters(
    ["{}/pipeline.yml".format(os.path.splitext(__file__)[0]), "pipeline.yml"])

INPUT_FORMATS = ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz")
INPUT_FORMATS = [os.path.join(PARAMS['fastq_path'], i) for i in INPUT_FORMATS]

# initiate pipeline from output of the pipeline_readqc.
REGEX_FORMATS = regex(r"{}/trimmed-(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)"\
        .format(PARAMS['fastq_path']))

CLUSTER_OPTIONS = dict(
    job_memory="8G",
    job_options='-l h_rt=240:0:0 -l h=!allen-606-19'
)
###############################################################################
###############################################################################
###############################################################################
# Alignment to a reference genome


@follows( mkdir("{}/bam".format(PARAMS['output_folder'])) )
@transform(
    INPUT_FORMATS,
    REGEX_FORMATS,
    r"{}/bam/\1.bam".format(PARAMS['output_folder'])
)
def mapReads(infiles, outfile):
    '''Map reads to the genome using BWA-MEM (output=SAM), convert to BAM,
    sort and index BAM file'''
    track = P.snip(os.path.basename(outfile), ".bam")
    m = PipelineMapping.BWAMEM(remove_unique=PARAMS["bwa_remove_non_unique"])
    statement = m.build((infiles,), outfile)
    P.run(statement, job_memory="16G", job_options='-l h_rt=240:0:0', job_threads=PARAMS["bwa_threads"])


###############################################################################
###############################################################################
###############################################################################
# Post-alignment QC


@transform(
    mapReads,
    regex(r"{}/bam/(\S+).bam".format(PARAMS['output_folder'])),
    r"{}/bam/\1.picard_stats".format(PARAMS['output_folder'])
)
def PicardAlignStats(infile, outfile):
    '''Run Picard CollectMultipleMetrics on each BAM file'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"]
    cluster_options = copy.copy(CLUSTER_OPTIONS)
    cluster_options['job_memory']='20G'
    PipelineMappingQC.buildPicardAlignmentStats(
            infile,
            outfile,
            genome,
            cluster_options,
    )


###############################################################################
###############################################################################
###############################################################################
# GATK


@follows(PicardAlignStats, mkdir("{}/gatk".format(PARAMS['output_folder'])))
@transform(
        mapReads,
        regex(r"{}/bam/(\S+).bam".format(PARAMS['output_folder'])),
        r"{}/gatk/\1.readgroups.bam".format(PARAMS['output_folder'])
)
def GATKReadGroups(infile, outfile):
    '''Reorders BAM according to reference fasta and adds read groups using
    GATK'''
    '''Reorders BAM according to reference fasta and add read groups using
    SAMtools, realigns around indels and recalibrates base quality
    scores using GATK

    '''

    track = re.sub(r'-\w+-\w+\.bam', '', os.path.basename(infile))
    tmpdir_gatk = P.get_temp_dir('.')
    job_threads = PARAMS["gatk_threads"]
    library = PARAMS["readgroup"]["library"]
    platform = PARAMS["readgroup"]["platform"]
    platform_unit = PARAMS["readgroup"]["platform_unit"]
    genome = PARAMS["human_ref"]
    PipelineExome.GATKReadGroups(infile, outfile, genome,
                                 library, platform,
                                 platform_unit, track,
                                 PARAMS['tmpdir'], CLUSTER_OPTIONS)
    IOTools.zap_file(infile)

###############################################################################
###############################################################################
###############################################################################
# Recalibration, 

@transform(GATKReadGroups,
           regex(r"{}/gatk/(\S+).readgroups.bam".format(PARAMS['output_folder'])),
           r"{}/gatk/\1.bqsr.bam".format(PARAMS['output_folder']))
def GATKBaseRecal(infile, outfile):
    '''recalibrates base quality scores using GATK'''
    intrack = P.snip(os.path.basename(infile), ".bam")
    outtrack = P.snip(os.path.basename(outfile), ".bam")
    dbsnp = PARAMS["gatk_dbsnp"]
    solid_options = PARAMS["gatk_solid_options"]
    genome = PARAMS["human_ref"]
    intervals = None
    padding = None
    PipelineExome.GATKBaseRecal(infile, outfile, genome, intervals,
                                padding, dbsnp, solid_options,
                                gatk_version=4, gatk=PARAMS['gatk']['executable'],
                                cluster_options = CLUSTER_OPTIONS)
    IOTools.zap_file(infile)

###############################################################################
###############################################################################
###############################################################################
# Coverage of targetted area


@transform(GATKBaseRecal,
           regex(r"{}/gatk/(\S+).bqsr.bam".format(PARAMS['output_folder'])),
           r"{}/gatk/\1.cov".format(PARAMS['output_folder']))
def buildCoverageStats(infile, outfile):
    '''Generate coverage statistics for regions of interest from a bed
    file using Picard'''
    picard_opts = '-Xmx9G -XX:+UseParNewGC -XX:+UseConcMarkSweepGC' % locals()
    job_threads = 3


    human_ref = PARAMS['human_ref']
    statement = '''picard %(picard_opts)s CollectWgsMetrics
    I=%(infile)s
    O=%(outfile)s
    R=%(human_ref)s''' % locals()
    P.run(
        statement,
        job_memory="16G",
        job_options='-l h_rt=240:0:0',
        job_threads=3
    )


# Whole genome takes a lot of time! subdivide into chroms
###############################################################################
###############################################################################
###############################################################################
# Variant Calling
# Whole genome takes a lot of time! subdivide into chroms
# Due to the way subdivide works, need to subdivide first, then do haplotypecall

@follows(
    buildCoverageStats,
    mkdir("{}/sentinel".format(PARAMS['output_folder']))
)
@subdivide(
    GATKBaseRecal,
    regex(r"\S+/(\S+).bqsr.bam"),
    r"{}/sentinel/\1.haplotype.chr*".format(PARAMS['output_folder']),
    r"{}/sentinel/\1.haplotype.chr".format(PARAMS['output_folder']),
)
def subdivideForHaplotypeCaller(infile, outfiles, stem):
    for chrom in PARAMS['chromosomes']:
        outfile = stem+chrom
        IOTools.touch_file(outfile)

@transform(
        subdivideForHaplotypeCaller,
        regex(r"\S+/(\S+).haplotype.chr(\w+)"),
        r"{}/variants/\1.haplotypeCaller.chr\2.g.vcf.gz".format(PARAMS['output_folder']),
        r"\1",
        r"\2",
)
def haplotypeCaller(sentinel, outfile, sample, intervals):
    '''Call SNVs and indels using GATK HaplotypeCaller in individuals'''
    genome = PARAMS["human_ref"]
    job_threads = PARAMS["gatk_threads"]
    dbsnp = PARAMS["gatk_dbsnp"]
    padding = None
    options = PARAMS["gatk_hc_options"]
    cluster_options = copy.copy(CLUSTER_OPTIONS)
    #cluster_options['job_memory']="12G"
    infile = '{}/gatk/{}.bqsr.bam'.format(PARAMS['output_folder'], sample)
    PipelineExome.haplotypeCaller(infile, outfile, genome, dbsnp,
                                  intervals, padding, options, wgs=True,
                                  gatk_version=4,
                                  gatk=PARAMS['gatk']['executable'],
                                  cluster_options=cluster_options)

###############################################################################
# import to database before merging gvcfs
@collate(
        haplotypeCaller,
        regex(r"{}/variants/\S+.haplotypeCaller.chr(\w+).g.vcf.gz".format(PARAMS['output_folder'])),
        r"{}/variants/gatkImportSentinel.chr\1".format(PARAMS['output_folder']),
        r"\1"
)
def gatkDBImport(infiles, outfile, chrom):
    outfolder = '{}/variants/dbi_{}'.format(PARAMS['output_folder'],chrom)
    PipelineExome.GenomicsDBImport(infiles, outfolder, chrom, PARAMS['gatk']['executable'],
                    cluster_options=CLUSTER_OPTIONS)
    IOTools.touch_file(outfile)


###############################################################################
# joined call gvcfs
@transform(
        gatkDBImport,
        regex(r"{}/variants/gatkImportSentinel.chr(\w+)".format(PARAMS['output_folder'])),
        r"{}/variants/\1.vcf.gz".format(PARAMS['output_folder']),
        r"\1"
        )
def genotypeGVCFs(infile, outfile, chrom):
    '''genotypeGVCFS'''
    # get dbi folder
    infolder = '{}/variants/dbi_{}'.format(PARAMS['output_folder'], chrom)
    PipelineExome.genotypeGVCFs(
        infolder,
        outfile,
        PARAMS['human_ref'],
        PARAMS['gatk']['executable'],
        gatk_version=4,
        From = 'gendb://',
        cluster_options = CLUSTER_OPTIONS,
    )
    # delete DB
    # need implementation

###############################################################################

@transform(
    genotypeGVCFs,
    regex(r"(\S+)/(\w+).vcf.gz"),
    r"\1/\2.snp_vqsr.recal",
    r"\1/\2.snp_vqsr.tranches",
    r"\2",
)
def variantRecalibratorSnps(infile, outfile, outfile2, chrom):
    '''Create variant recalibration file'''
    # skip MT
    if chrom == 'MT':
        IOTools.touch_file(outfile)
        IOTools.touch_file(outfile2)
        return None

    cluster_options = copy.copy(CLUSTER_OPTIONS)
    cluster_options['job_memory'] = '12G'
    genome = PARAMS["human_ref"]
    dbsnp = PARAMS["gatk_dbsnp"]
    job_threads = PARAMS["gatk_threads"]
    track = P.snip(outfile, ".recal")
    kgenomes = PARAMS["gatk_kgenomes"]
    hapmap = PARAMS["gatk_hapmap"]
    omni = PARAMS["gatk_omni"]
    mode = 'SNP'
    gatk = PARAMS['gatk']['executable']
    PipelineExome.variantRecalibrator(infile, outfile, genome, mode, dbsnp,
                                      kgenomes, hapmap, omni, gatk=gatk, gatk_version=4,
                                      cluster_options=cluster_options)

###############################################################################

@follows(variantRecalibratorSnps)
@transform(
    genotypeGVCFs,
    regex(r"(\S+)/(\w+).vcf.gz"),
    add_inputs(r"\1/\2.snp_vqsr.recal",
              r"\1/\2.snp_vqsr.tranches"),
    r"\1/\2.snp_vqsr.vcf.gz",
    r"\2",
)
def applyVariantRecalibrationSnps(infiles, outfile, chrom):
    '''Perform variant quality score recalibration using GATK '''
    vcf, recal, tranches = infiles
    if chrom == 'MT':
        shutil.copyfile(vcf, outfile)
        shutil.copyfile(vcf+'.tbi', outfile+'.tbi')
        return None

    genome = PARAMS["human_ref"]
    mode = 'SNP'
    gatk = PARAMS['gatk']['executable']
    PipelineExome.applyVariantRecalibration(vcf, recal, tranches, outfile,
                                            genome, mode, gatk=gatk, gatk_version=4,
                                            cluster_options=CLUSTER_OPTIONS)

###############################################################################
# Indel recalibration


@transform(
    applyVariantRecalibrationSnps,
    regex(r"(\S+)/(\w+).snp_vqsr.vcf.gz"),
    r"\1/\2.vqsr.recal",
    r"\1/\2.vqsr.tranches",
    r"\2",
)
def variantRecalibratorIndels(infile, outfile, outfile2, chrom):
    '''Create variant recalibration file'''
    # skip MT
    if chrom == 'MT':
        IOTools.touch_file(outfile)
        IOTools.touch_file(outfile2)
        return None
    cluster_options = copy.copy(CLUSTER_OPTIONS)
    cluster_options['job_memory'] = '10G'
    genome = PARAMS["human_ref"]
    job_threads = PARAMS["gatk_threads"]
    track = P.snip(outfile, ".recal")
    mills = PARAMS["gatk_mills"]
    dbsnp = PARAMS["gatk_dbsnp"]
    mode = 'INDEL'
    gatk = PARAMS['gatk']['executable']
    PipelineExome.variantRecalibrator(infile, outfile, genome, mode, gatk=gatk,
                                      gatk_version=4, mills=mills,
                                      cluster_options=cluster_options)

###############################################################################


@follows(variantRecalibratorIndels)
@transform(
    applyVariantRecalibrationSnps,
    regex(r"(\S+)/(\w+).snp_vqsr.vcf.gz"),
    add_inputs(r"\1/\2.vqsr.recal",
              r"\1/\2.vqsr.tranches"),
    r"\1/\2.vqsr.vcf.gz",
    r"\2",
)
def applyVariantRecalibrationIndels(infiles, outfile, chrom):
    '''Perform variant quality score recalibration using GATK '''
    vcf, recal, tranches = infiles
    if chrom == 'MT':
        shutil.copy(vcf, outfile)
        shutil.copy(vcf+'.tbi', outfile+'.tbi')
        return None
    genome = PARAMS["human_ref"]
    mode = 'INDEL'
    gatk = PARAMS['gatk']['executable']
    PipelineExome.applyVariantRecalibration(vcf, recal, tranches,
                                            outfile, genome, mode,
                                            gatk=gatk, gatk_version=4,
                                            cluster_options=CLUSTER_OPTIONS)

###############################################################################

@follows(mkdir('{}/vep'.format(PARAMS['output_folder'])))
@transform(
    applyVariantRecalibrationIndels,
    regex(r"(\S+)/variants/(\w+).vqsr.vcf.gz"),
    r"\1/vep/\2.vep.vcf.gz",
    r"\2.vep.vcf"
)
def vep(infile, outfile, tmpfile):
    ''' vep '''
    vep_dict = PARAMS['annotation']['vep']
    tmpfile = os.path.join(PARAMS['tmpdir'], tmpfile)
    cluster_options = copy.copy(CLUSTER_OPTIONS)
    cluster_options['job_memory'] = '9G'
    statement = '''
        {executable} -i {infile} {options} --dir {dir}
        -a {build} --fork {fork} -o {tmpfile};
        bgzip -c {tmpfile} > {outfile};
        tabix -p vcf {outfile};
        rm {tmpfile}
    '''.format(**vep_dict, **locals())
    P.run(statement, **cluster_options)


###############################################################################

@follows(mkdir('{}/cadd'.format(PARAMS['tmpdir'])))
@transform(
    vep,
    regex(r"\S+/(\w+).vep.vcf.gz"),
    r"{}/cadd/\1.cadd.gz".format(PARAMS['tmpdir']),
)
def cadd(infile, outfile):
    cadd = PARAMS['cadd']['executable']
    statement = '''{cadd} {infile} {outfile};
        tabix -p vcf {outfile}'''.format(**locals())
    P.run(statement, job_memory="10G")

'''population annotation'''
'''doesn't work on MT'''
@follows(mkdir('{}/pop'.format(PARAMS['output_folder'])))
@transform(
    vep,
    regex('\S+/(\w+).vep.vcf.gz'),
    r'{}/pop/\1.vcf.gz'.format(PARAMS['output_folder']),
    r'\1.vcf',
    r'\1',
)
def addPop(infile, outfile, basename, chrom):
    if chrom == 'MT':
        #IOTools.shadow_file(outfile)
        shutil.copyfile(infile, outfile)
        return None
    temp_dir = P.get_temp_dir(PARAMS['tmpdir'])
    temp_outfile = os.path.join(temp_dir, basename)
    PipelineExome.add_pop_freqs(infile, temp_outfile, PARAMS, submit=True, job_memory="10G", job_options='-l h_rt=240:0:0')
    statement = 'bgzip -c {} >{}'.format(temp_outfile, outfile) 
    P.run(statement, job_memory="2G", job_options='-l h_rt=240:0:0')
    os.unlink(temp_outfile)
    os.rmdir(temp_dir)
    #IOTools.zap_file(infile)

@follows(
    mkdir('{}/anno'.format(PARAMS['output_folder'])),
    cadd
)
@transform(
    addPop,
    regex('\S+/(\w+).vcf.gz'),
    r'{}/anno/\1.vcf.gz'.format(PARAMS['output_folder']),
    r'\1.vcf',
    r'\1',
)
def addCadd(infile, outfile, basename, chrom):
    cadd_file = '{}/cadd/{}.cadd.gz'.format(PARAMS['tmpdir'],chrom)
    temp_dir = P.get_temp_dir(PARAMS['tmpdir'])
    temp_outfile = os.path.join(temp_dir, basename)
    PipelineExome.add_cadd(infile, temp_outfile, cadd_file, PARAMS, submit=True, job_memory="10G", job_options='-l h_rt=240:0:0')
    statement = 'bgzip -c {temp_outfile} > {outfile}; tabix -p vcf {outfile}'.format(**locals()) 
    P.run(statement, job_memory="2G", job_options='-l h_rt=240:0:0')
    os.unlink(temp_outfile)
    os.rmdir(temp_dir)
    #IOTools.zap_file(infile)

###############################################################################
###############################################################################
# GATK CNV (beta)
# follow https://gatkforums.broadinstitute.org/dsde/discussion/11682/
# subject to change since this is Beta.
# not yet implemented!!! Don't use the follow code

@transform(
    genotypeGVCFs,
    regex(r"\S+/(\S+).vcf.gz"),
    r"\1.interval_list",
    r"\1",
)
def generateIntervalList(infile, outfile, chrom):
    '''Generate interval_list for cnv inference'''
    # get window size
    window_size = PARAMS['gatk']['cnv_window_size'].get(
            chrom,
            PARAMS['gatk']['cnv_window_size']['default']
    )
    # get contig length
    contig_length = 0
    with IOTools.open_file(infile) as inf:
        for line in inf:
            if line.startswith('##contig=<ID={}'.format(chrom)):
                contig_length = int(line.split('length=')[1].split('>')[0])
            break
    if contig_length == 0:
        msg = 'Contig is of length 0: {}'.format(chrom)
        raise ValueError(msg)
    with open(outfile, 'w') as outf:
        intervals = list(np.arange(1, contig_length, window_size))
        intervals.append(contig_length)
        for ind, end in enumerate(intervals):
            if not ind:
                continue
            start = intervas[ind-1]
            line = '{chrom}:{start}-{end}\n'.format(**locals())
            out.write(line)
    

@follows(
    generateIntervalList,
    mkdir("{}/cnv_temp".format(PARAMS['output_folder']))
)
@transform(
    subdivideForHaplotypeCaller,
    regex(r"\S+/(\S+).haplotype.chr(\w+)"),
    r"{}/gatk/\1.bqsr.bam",
    r"{}/cnv_temp/\1.counts.chr\2.hdf5".format(PARAMS['output_folder']),
    r"\2",
)
def CollectReadCounts(sentinel, infile, outfile, chrom):
    '''collect read counts from bam'''
    gatk_options = '-Xmx1G'
    gatk = PARAMS['gatk']['executable']
    interval_file = '{}.interval_list'.format(chrom)
    statement = '''
        {gatk} --java-options {gatk_options} CollectReadCounts 
        -I {infile}
        -L {interval_file}
        --interval-merging-rule OVERLAPPING_ONLY
        -O {outfile}
    '''.format(**locals())
    P.run(statement, **CLUSTER_OPTIONS)
####
# Note that the following step relies on gcnvkernel, which relies on pymc3 3.1
# If you go to gatk4.0's path and install gcnvkernel from there,
# and conda install pymc3, pymc3 will be greater than 3.1, and gcnvkernel will
# panic. pip install pymc3==3.1 instead
# Don't use gatk4.0's yml to update your env! it will mess up everything dear to you!!
@collate(
    CollectReadCounts,
    regex(r"(\S+)/\S+.counts.(\w+).hdf5"),
    r"\1/\2.DGCP",
    r"\1",
    r"\2",
)
def DetermineGermlineContigPloidy(infiles, outfile, output_folder, chrom):
    '''Call Ploidy'''
    '''failed on MT, maybe sampled reads are too low?'''
    cluster_options = copy.copy(CLUSTER_OPTIONS)
    cluster_options['job_memory']="7G"
    gatk_options = '-Xmx1G'
    gatk = PARAMS['gatk']['executable']
    ploidy_prior = PARAMS['gatk']['contig_ploidy_prior']
    statement = '''
        {gatk} --java-options {gatk_options} DetermineGermlineContigPloidy 
    '''.format(**locals())

    # comment out the following block when not testing!
    if chrom == 'chrMT':
        IOTools.touch_file(outfile)
        return None

    for infile in infiles:
        statement += '-I {} '.format(infile)
    statement += '''
        --contig-ploidy-priors {ploidy_prior}
        --output {output_folder}
        --output-prefix {chrom}
    '''.format(**locals())
    P.run(statement, **cluster_options)
    IOTools.touch_file(outfile)

@follows(
    mkdir('{}/cnv'.format(PARAMS['output_folder'])),
    DetermineGermlineContigPloidy,
)
@collate(
    CollectReadCounts,
    regex(r"\S+/\S+.counts.chr(\w+).hdf5"),
    r"{}/cnv/chr\1.GCC".format(PARAMS['output_folder']),
    r"{}/cnv".format(PARAMS['output_folder']),
    r"\1",
)
def GermlineCNVCaller(infiles, outfile, out_dir, chrom):
    '''Germline CNV caller'''
    '''takes inputs from both ploidycaller and count'''
    # note that test doesn't work on MT
    cluster_options = copy.copy(CLUSTER_OPTIONS)
    cluster_options['job_memory']="12G"
    gatk_options = '-Xmx1G'
    gatk = PARAMS['gatk']['executable']
    ploidy_prior = PARAMS['gatk']['contig_ploidy_prior']
    if chrom == 'MT':
        IOTools.touch_file(outfile)
        return None
    ploidyCallPath = '{}/cnv_temp/chr{}-calls'.format(PARAMS['output_folder'], chrom)
    output_dir = '{}/chr{}'.format(out_dir, chrom)
    statement = '''
        {gatk} --java-options {gatk_options} GermlineCNVCaller
        --run-mode COHORT
        --contig-ploidy-calls {ploidyCallPath}
        --output {out_dir}
        --output-prefix {chrom}
    '''.format(**locals())
    for infile in infiles:
        statement += '-I {} '.format(infile)
    P.run(statement, **cluster_options)
    IOTools.touch_file(outfile)



def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
