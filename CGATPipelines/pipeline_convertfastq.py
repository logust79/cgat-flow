# The aim of this small ruffus task is to convert bam files to fastq files

from ruffus import transform, merge, follows, mkdir, regex, suffix, \
jobs_limit, subdivide, collate, active_if, originate

import CGATPipelines.Pipeline as P
import sys
import os

P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'paired_end': False})

PARAMS = P.PARAMS


INFILES = "*.bam"

@transform(INFILES,
           suffix(".bam"),
           r"\1_1.fastq.gz")
def bam2fastq(infile, outfile):

    outfile = infile.strip(".bam")

    statement = '''cat %(infile)s | cgat bam2fastq %(outfile)s_1.fastq.gz %(outfile)s_2.fastq.gz '''

    P.run()

@follows(bam2fastq)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
