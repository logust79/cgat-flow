"""====================
ReadQc pipeline
====================


The readqc pipeline imports unmapped reads from one or more input
files and performs basic quality control steps. The pipeline performs
also read pre-processing such as quality trimming or adaptor removal.

Quality metrics are based on the fastqc tools, see
http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/ for further
details.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning`
on general information how to use CGAT pipelines.

When pre-processing reads before mapping, the workflow of the pipeline
is as follows:

1. Run the ``full`` target to perform initial QC on the raw data. Then
   build the report (``build_report`` target).

2. Inspect the output to decide if and what kind of pre-processing is
   required.

3. Edit the configuration file ``pipeline.yml`` to activate
   pre-processing and parameterize it appropriately. Note that
   parameters can be set on a per-sample basis.

4. Rerun the ``full`` target and ``build_report`` targets. The data
   will now be processed and additional QC will be performed on the
   processed data. Note that all the processed data will be found in
   the :file:`processed.dir` directory.


.. note::

   If you set the option auto_remove, you will ned to run the
   pipeline at least once without any pre-processors.

Configuration
-------------

See :file:`pipeline.yml` for setting configuration values affecting
the workflow (pre-processing or no pre-processing) and options for
various pre-processing tools.

Input
-----

Reads are imported by placing files or linking to files in the :term:
`working directory`.

The default file format assumes the following convention:

   <sample>.<suffix>

The ``suffix`` determines the file type. The following suffixes/file
types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the :file:
   `fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq.2.gz
   Paired-end reads in fastq format.
   The two fastq files must be sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input files.
   Thus it might be difficult to mix different formats.

Pipeline output
----------------

The major output is a set of HTML pages and plots reporting on the quality of
the sequence archive.

Example
=======

Example data is available at

https://www.cgat.org/downloads/public/cgatpipelines/pipeline_test_data/test_readqc.tgz

To run the example, simply unpack and untar::

   wget -qO- https://www.cgat.org/downloads/public/cgatpipelines/
             pipeline_test_data/test_readqc.tgz | tar -xvz
   cd test_readqc
   python <srcdir>/pipeline_readqc.py make full

Code
====

To add a new pre-processing tool, the following changes are required:

1. Add a new tool wrapper to :module:`PipelinePreprocess`. Derive
   the wrapper from :class:`PipelinePreprocess.ProcessTools`. Make
   sure to add a corresponding entry in the Requirements section of
   the module.

2. Add the tool to the task :func:`processReads` in this module.

Requirements:

"""

# import ruffus
from ruffus import transform, merge, follows, mkdir, regex, suffix, \
    jobs_limit, subdivide, collate, active_if, originate

# import useful standard python modules
import sys
import os
import re
import shutil
import sqlite3
import glob

# import modules from the CGAT code collection
import CGATCore.Experiment as E
import CGATPipelines.PipelineMapping as PipelineMapping
from CGATCore import Pipeline as P
import CGATPipelines.PipelineReadqc as PipelineReadqc
import CGATPipelines.PipelinePreprocess as PipelinePreprocess
import CGATCore.IOTools as IOTools
from CGATPipelines.Report import run_report

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# define input files and preprocessing steps
# list of acceptable input formats
INPUT_FORMATS = ["*.fastq.1.gz", "*.fastq.gz",
                 "*.sra", "*.csfasta.gz", "*.remote"]

# Regular expression to extract a track from an input file. Does not preserve
# a directory as part of the track.
REGEX_TRACK = r"([^/]+).(fastq.1.gz|fastq.gz|sra|csfasta.gz|remote)"

# Regular expression to extract a track from both processed and unprocessed
# files
REGEX_TRACK_BOTH = \
    r"({}/processed.dir/)*([^/]+)\.(fastq.1.gz|fastq.gz|sra|csfasta.gz|remote)"\
    .format(PARAMS['exportdir'])

SEQUENCEFILES_REGEX = r"(\S+).(?P<suffix>fastq.1.gz|fastq.gz|sra|csfasta.gz|remote)"

# some cluster options are not passed properly, 
#  so I wrote this to come around the issue
CLUSTER_OPTIONS = dict(
    job_options = '-l h_rt=240:0:0 -l h=!allen-606-19',
    job_memory = '16G',
)

def connect():
    '''
    Setup a connection to an sqlite database
    '''

    dbh = sqlite3.connect(PARAMS['database'])
    return dbh


@transform(INPUT_FORMATS, regex("(.*)"), r"\1")
def unprocessReads(infiles, outfiles):
    """dummy task - no processing of reads."""


# if preprocess tools are specified, preprocessing is done on output that has
# already been generated in the first run
if PARAMS.get("preprocessors", None):
    if PARAMS["auto_remove"]:
        # check if fastqc has been run
        for x in IOTools.flatten([glob.glob(y) for y in INPUT_FORMATS]):
            f = re.match(REGEX_TRACK, x).group(1) + ".fastqc"
            if not os.path.exists(f):
                raise ValueError(
                    "file %s missing, "
                    "you need to run the pipeline once before "
                    "specifying 'auto_remove'" % f)

        @follows(mkdir("fasta.dir"))
        @transform(unprocessReads,
                   regex(SEQUENCEFILES_REGEX),
                   r"fasta.dir/\1.fasta")
        def makeAdaptorFasta(infile, outfile):
            '''Make a single fasta file for each sample of all contaminant adaptor
            sequences for removal
            '''

            print(infile)
            print(REGEX_TRACK)

            PipelinePreprocess.makeAdaptorFasta(
                infile=infile,
                outfile=outfile,
                track=re.match(REGEX_TRACK, infile).groups()[0],
                dbh=connect(),
                contaminants_file=PARAMS['contaminants_path'])

        @merge(makeAdaptorFasta, "contaminants.fasta")
        def aggregateAdaptors(infiles, outfile):
            '''
            Collate fasta files into a single contaminants file for
            adapter removal.
            '''
            tempfile = P.get_temp_filename()
            infiles = " ".join(infiles)

            statement = """
            cat %(infiles)s | fastx_reverse_complement > %(tempfile)s;
            cat %(tempfile)s %(infiles)s | fastx_collapser > %(outfile)s;
            rm -f %(tempfile)s
            """
            P.run(statement, **CLUSTER_OPTIONS)

    else:
        @follows(mkdir("fasta.dir"))
        @transform(INPUT_FORMATS,
                   regex(SEQUENCEFILES_REGEX),
                   r"fasta.dir/\1.fasta")
        def aggregateAdaptors(infile, outfile):
            IOTools.touch_file(outfile)

    @follows(mkdir("{}/processed.dir".format(PARAMS['exportdir'])),
             aggregateAdaptors)
    @subdivide(INPUT_FORMATS,
               regex(SEQUENCEFILES_REGEX),
               r"{}/processed.dir/trimmed-\1.fastq*.gz".format(PARAMS['exportdir']))
    def processReads(infile, outfiles):
        '''process reads from .fastq and other sequence files.
        '''
        trimmomatic_options = PARAMS["trimmomatic_options"]

        if PARAMS["auto_remove"]:
            trimmomatic_options = " ILLUMINACLIP:%s:%s:%s:%s:%s:%s " % (
                "contaminants.fasta",
                PARAMS["trimmomatic_mismatches"],
                PARAMS["trimmomatic_p_thresh"],
                PARAMS["trimmomatic_c_thresh"],
                PARAMS["trimmomatic_min_adapter_len"],
                PARAMS["trimmomatic_keep_both_reads"]) + trimmomatic_options

        elif PARAMS["trimmomatic_adapter"]:
            trimmomatic_options = " ILLUMINACLIP:%s:%s:%s:%s:%s:%s " % (
                PARAMS["trimmomatic_adapter"],
                PARAMS["trimmomatic_mismatches"],
                PARAMS["trimmomatic_p_thresh"],
                PARAMS["trimmomatic_c_thresh"],
                PARAMS["trimmomatic_min_adapter_len"],
                PARAMS["trimmomatic_keep_both_reads"]) + trimmomatic_options

        job_threads = PARAMS["threads"]
        job_memory = "12G"

        track = re.match(REGEX_TRACK, infile).groups()[0]

        m = PipelinePreprocess.MasterProcessor(
            save=PARAMS["save"],
            summarize=PARAMS["summarize"],
            threads=PARAMS["threads"],
            qual_format=PARAMS['qual_format'])

        for tool in P.as_list(PARAMS["preprocessors"]):

            if tool == "fastx_trimmer":
                m.add(PipelinePreprocess.FastxTrimmer(
                    PARAMS["fastx_trimmer_options"],
                    threads=PARAMS["threads"]))
            elif tool == "trimmomatic":
                m.add(PipelinePreprocess.Trimmomatic(
                    trimmomatic_options,
                    threads=PARAMS["threads"]))
            elif tool == "sickle":
                m.add(PipelinePreprocess.Sickle(
                    PARAMS["sickle_options"],
                    threads=PARAMS["threads"]))
            elif tool == "trimgalore":
                m.add(PipelinePreprocess.Trimgalore(
                    PARAMS["trimgalore_options"],
                    threads=PARAMS["threads"]))
            elif tool == "flash":
                m.add(PipelinePreprocess.Flash(
                    PARAMS["flash_options"],
                    threads=PARAMS["threads"]))
            elif tool == "reversecomplement":
                m.add(PipelinePreprocess.ReverseComplement(
                    PARAMS["reversecomplement_options"]))
            elif tool == "pandaseq":
                m.add(PipelinePreprocess.Pandaseq(
                    PARAMS["pandaseq_options"],
                    threads=PARAMS["threads"]))
            elif tool == "cutadapt":
                cutadapt_options = PARAMS["cutadapt_options"]
                if PARAMS["auto_remove"]:
                    cutadapt_options += " -a file:contaminants.fasta "
                m.add(PipelinePreprocess.Cutadapt(
                    cutadapt_options,
                    threads=PARAMS["threads"],
                    untrimmed=PARAMS['cutadapt_reroute_untrimmed'],
                    process_paired=PARAMS["cutadapt_process_paired"]))
            else:
                raise NotImplementedError("tool '%s' not implemented" % tool)

        statement = m.build((infile,), "{}/processed.dir/trimmed-".format(PARAMS['exportdir']), track)
        P.run(
                statement,
                job_memory=job_memory,
                job_threads=job_threads,
                job_options=CLUSTER_OPTIONS['job_options'])

else:
    @follows(mkdir("{}/processed.dir".format(PARAMS['exportdir'])))
    def processReads():
        """dummy task - no processing of reads."""


@active_if(PARAMS["general_reconcile"] == 1)
@follows(mkdir("{}/reconciled.dir".format(PARAMS['exportdir'])))
@transform(processReads, regex(
    r"{}/processed.dir\/trimmed-(.*)\.fastq\.1\.gz".format(PARAMS['exportdir'])),
    r"{}/reconciled.dir/trimmed-\1.fastq.1.gz".format(PARAMS['exportdir']))
def reconcileReads(infile, outfile):
    if PARAMS["general_reconcile"] == 1:
        in1 = infile
        in2 = infile.replace(".fastq.1.gz", ".fastq.2.gz")
        outfile = outfile.replace(".fastq.1.gz",  "")
        job_threads = PARAMS["threads"]
        job_memory = "8G"
        statement = """cgat fastqs2fastqs
            --method=reconcile
            --output-filename-pattern=%(outfile)s.fastq.%%s.gz
            %(in1)s %(in2)s"""

        P.run(statement, **CLUSTER_OPTIONS)


@follows(reconcileReads)
@follows(mkdir(PARAMS["exportdir"]),
         mkdir(os.path.join(PARAMS["exportdir"], "fastqc")))
@transform((unprocessReads, processReads),
           regex(REGEX_TRACK),
           r"\1.fastqc")
def runFastqc(infiles, outfile):
    '''run Fastqc on each input file.

    convert sra files to fastq and check mapping qualities are in
    solexa format.  Perform quality control checks on reads from
    .fastq files.
    '''
    # MM: only pass the contaminants file list if requested by user,
    # do not make this the default behaviour
    if PARAMS['use_custom_contaminants']:
        m = PipelineMapping.FastQc(nogroup=PARAMS["readqc_no_group"],
                                   outdir=PARAMS["exportdir"] + "/fastqc",
                                   contaminants=PARAMS['contaminants_path'],
                                   qual_format=PARAMS['qual_format'])
    else:
        m = PipelineMapping.FastQc(nogroup=PARAMS["readqc_no_group"],
                                   outdir=PARAMS["exportdir"] + "/fastqc",
                                   qual_format=PARAMS['qual_format'])

    if PARAMS["general_reconcile"] == 1:
        infiles = infiles.replace("{}/processed.dir/trimmed".format(PARAMS['exportdir']),
                                  "{}/reconciled.dir/trimmed".format(PARAMS['exportdir']))

    statement = m.build((infiles,), outfile)
    P.run(statement, **CLUSTER_OPTIONS)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(runFastqc, suffix(".fastqc"), "_fastqc.load")
def loadFastqc(infile, outfile):
    '''load FASTQC stats into database.'''
    track = P.snip(infile, ".fastqc")
    filename = os.path.join(
        PARAMS["exportdir"], "fastqc", track + "*_fastqc", "fastqc_data.txt")
    PipelineReadqc.loadFastqc(filename,
                              database_url=PARAMS["database"]["url"])
    IOTools.touch_file(outfile)


@follows(mkdir(PARAMS["exportdir"]),
         mkdir(os.path.join(PARAMS["exportdir"], "fastq_screen")))
@active_if(PARAMS["fastq_screen_run"] == 1)
@transform((unprocessReads, processReads),
           regex(REGEX_TRACK_BOTH),
           r"%s/fastq_screen/\2.fastqscreen" % PARAMS['exportdir'])
def runFastqScreen(infiles, outfile):
    '''run FastqScreen on input files.'''

    # variables required for statement built by FastqScreen()
    tempdir = P.get_temp_dir(".")
    outdir = os.path.join(PARAMS["exportdir"], "fastq_screen")

    # configure job_threads with fastq_screen_options from PARAMS
    job_threads = re.findall(r'--threads \d+', PARAMS['fastq_screen_options'])
    if len(job_threads) != 1:
        raise ValueError("Wrong number of threads for fastq_screen")

    job_threads = int(re.sub(r'--threads ', '', job_threads[0]))

    # Create fastq_screen config file in temp directory
    # using parameters from Pipeline.yml
    with IOTools.open_file(os.path.join(tempdir, "fastq_screen.conf"),
                           "w") as f:
        for i, k in list(PARAMS.items()):
            if i.startswith("fastq_screen_database"):
                f.write("DATABASE\t%s\t%s\n" % (i[22:], k))

    m = PipelineMapping.FastqScreen()
    statement = m.build((infiles,), outfile)
    P.run(statement, job_memory="8G", job_options='-l h_rt=240:0:0')
    shutil.rmtree(tempdir)
    IOTools.touch_file(outfile)


@active_if(PARAMS["fastq_screen_run"] == 1)
@merge(runFastqc, ["fastqc_basic_statistics.tsv"])
def summarizeFastqc(infiles, outfiles):
    all_files = []

    for infile in infiles:
        track = P.snip(infile, ".fastqc")
        all_files.extend(glob.glob(
            os.path.join(
                PARAMS["exportdir"],
                "fastqc",
                track + "*_fastqc", "fastqc_data.txt")))

    dfs = PipelineReadqc.read_fastqc(
        all_files,
        track_regex="([^/]+).fastq.(\d+)")
    for key, df in dfs.items():
        fn = re.sub("basic_statistics", key, outfiles[0])
        E.info("writing to {}".format(fn))
        with IOTools.open_file(fn, "w") as outf:
            df.to_csv(outf, sep="\t", index=True)


@merge(runFastqScreen,
       ["fastqscreen_summary.tsv", "fastcscreen_details.tsv"])
def summarizeFastqScreen(infiles, outfiles):
    all_files = []
    for infile in infiles:
        all_files.extend(glob.glob(IOTools.snip(infile, "screen") + "*_screen.txt"))
    if len(all_files) == 0:
        E.warn("no fastqcscreen results to concatenate")
        for x in outfiles:
            IOTools.touch_file(x)
        return
    df_summary, df_details = PipelineReadqc.read_fastq_screen(
        all_files,
        track_regex="([^/]+).fastq.(\d+)")
    df_summary.to_csv(outfiles[0], sep="\t", index=True)
    df_details.to_csv(outfiles[1], sep="\t", index=True)


@merge(runFastqc, "status_summary.tsv.gz")
def buildFastQCSummaryStatus(infiles, outfile):
    '''load fastqc status summaries into a single table.'''
    exportdir = os.path.join(PARAMS["exportdir"], "fastqc")
    PipelineReadqc.buildFastQCSummaryStatus(infiles, outfile, exportdir)


@follows(loadFastqc)
@merge(runFastqc, "basic_statistics_summary.tsv.gz")
def buildFastQCSummaryBasicStatistics(infiles, outfile):
    '''load fastqc summaries into a single table.'''
    exportdir = os.path.join(PARAMS["exportdir"], "fastqc")
    PipelineReadqc.buildFastQCSummaryBasicStatistics(infiles, outfile,
                                                     exportdir)


@follows(mkdir("{}/experiment.dir".format(PARAMS['exportdir'])), loadFastqc)
@collate(runFastqc,
         regex("({}/processed.dir/)*(.*)-([^-]*).fastqc".format(PARAMS['exportdir'])),
         r"{}/experiment.dir/\2_per_sequence_quality.tsv".format(PARAMS['exportdir']))
def buildExperimentLevelReadQuality(infiles, outfile):
    """
    Collate per sequence read qualities for all replicates per experiment.
    Replicates are the last part of a filename, eg. Experiment-R1,
    Experiment-R2, etc.
    """
    exportdir = os.path.join(PARAMS["exportdir"], "fastqc")
    PipelineReadqc.buildExperimentReadQuality(infiles, outfile, exportdir)


@collate(buildExperimentLevelReadQuality,
         regex("(.+)/(.+)_per_sequence_quality.tsv"),
         r"\1/experiment_per_sequence_quality.tsv")
def combineExperimentLevelReadQualities(infiles, outfile):
    """
    Combine summaries of read quality for different experiments
    """
    infiles = " ".join(infiles)
    statement = ("cgat combine_tables "
                 "  --log=%(outfile)s.log "
                 "  --regex-filename='.+/(.+)_per_sequence_quality.tsv' "
                 "%(infiles)s"
                 "> %(outfile)s")
    P.run(statement, **CLUSTER_OPTIONS)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(combineExperimentLevelReadQualities,
           regex(".+/(.+).tsv"),
           r"\1.load")
def loadExperimentLevelReadQualities(infile, outfile):
    P.load(infile, outfile)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform((buildFastQCSummaryStatus, buildFastQCSummaryBasicStatistics),
           suffix(".tsv.gz"), ".load")
def loadFastqcSummary(infile, outfile):
    P.load(infile, outfile, options="--add-index=track")


@follows(loadFastqc,
         loadFastqcSummary,
         loadExperimentLevelReadQualities,
         summarizeFastqScreen,
         summarizeFastqc)
def full():
    pass


@follows(mkdir("MultiQC_report.dir"))
@originate("MultiQC_report.dir/multiqc_report.html")
def renderMultiqc(infile):
    '''build mulitqc report'''

    statement = (
        "export LANG=en_GB.UTF-8 && "
        "export LC_ALL=en_GB.UTF-8 && "
        "multiqc . -f && "
        "mv multiqc_report.html MultiQC_report.dir/")

    P.run(statement, **CLUSTER_OPTIONS)


@follows(renderMultiqc)
@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting documentation build process from scratch")
    run_report(clean=True)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
