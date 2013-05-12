from collections import Counter, defaultdict, namedtuple
from itertools import groupby
import logging

import annotation
import gff
import pysam
import pyximport; pyximport.install()

import intervaltree


logging.basicConfig(level=logging.INFO)


def read_samfile(samfile, paired=False):
    '''
    expects pysam.Samfile

    if paired=True, assumes alignments are ordered by read
    '''

    def read():
        for alignment in samfile:
            reference_ID = ''
            if alignment.tid != -1:
                reference_ID = samfile.getrname(alignment.tid)
            start = alignment.pos + 1
            end = alignment.aend
            read_ID = alignment.qname
            tags = dict(alignment.tags)
            mismatches = tags.get('NM', 0)
            strand = '-' if alignment.is_reverse else '+'
            yield Alignment(reference_ID, start, end, read_ID,
                            mismatches, strand)

    def read_paired():
        alignments = read()
        while True:
            a = alignments.next()
            b = alignments.next()
            yield PairedAlignment(a, b)

    if paired:
        return read_paired()
    else:
        return read()


class Alignment(object):
    def __init__(self, reference_ID, start, end, read_ID,
                 mismatches, strand):

        self.reference_ID = reference_ID
        self.start = start
        self.end = end
        self.read_ID = read_ID
        self.mismatches = mismatches
        self.strand = strand


class PairedAlignment(object):
    def __init__(self, a, b):
        self.a = a
        self.b = b

    @property
    def start(self):
        if self.a.start < self.b.start:
            return self.a.start
        else:
            return self.b.start

    @property
    def end(self):
        if self.b.end > self.a.end:
            return self.b.end
        else:
            return self.a.end

    @property
    def reference_ID(self):
        return self.a.reference_ID

    @property
    def read_ID(self):
        return self.a.read_ID

    @property
    def mismatches(self):
        return self.a.mismatches + self.b.mismatches


def filter_read_alignments_by_reference(alignments, exclude):
    '''assumes alignments are sorted by read'''

    for read_alignments in group_alignments_by_read(alignments):
        read_alignments = list(read_alignments)
        excluded = False

        for al in read_alignments:
            if al.reference_ID in exclude:
                excluded = True
                break

        if not excluded:
            for al in read_alignments:
                yield al


def group_alignments_by_read(alignments):
    '''assumes alignments are sorted by read'''
    return (als for ID, als in groupby(alignments, lambda al: al.read_ID))


GeneAlignment = namedtuple('GeneAlignment', 'gene alignment')

def GeneAlignmentFinder(genes, transcripts):
    gene_index = RegionIndex(genes)
    transcripts_by_ID = {t.ID: t for t in transcripts}

    def find(alignments):
        for alignment in alignments:

            transcript = transcripts_by_ID.get(alignment.reference_ID)

            # Did it align to a transcript?
            if transcript:
                yield GeneAlignment(transcript.gene, alignment)

            else:
                # Find which gene(s) this alignment overlaps by position.
                for gene in gene_index.find(alignment.reference_ID, alignment.strand,
                                            alignment.start, alignment.end):

                    yield GeneAlignment(gene, alignment)
    return find


def count_genes(gene_alignments):
    '''assumes alignments are sorted by read'''

    gene_counts = Counter()

    # Group gene alignments by read ID
    grouped = groupby(gene_alignments, lambda gal: gal.alignment.read_ID)
    for read_ID, group in grouped:

        # Get the unique set of gene(s) that this read aligned to
        genes = {gal.gene for gal in group}

        # Only count reads that align unambiguously
        # i.e. there should be only one gene in the unique set
        if len(genes) == 1:
            gene = genes.pop()
            gene_counts[gene] += 1

    return gene_counts


class RegionIndex(object):
    def __init__(self, regions):
        self.index = defaultdict(intervaltree.IntervalTree)

        for region in regions:
            key = (region.reference.name, region.strand)
            self.index[key].insert(region.start, region.end, region)

    def find(self, reference_name, start, end, strand):
        key = (reference_name, strand)
        return self.index[key].find(start, end)


class Transcript(annotation.Transcript):
    def __init__(self, ID):
        super(Transcript, self).__init__()
        self.ID = ID


class AnnotationBuilder(annotation.AnnotationBuilder):

    aliases = {
        'reference': ['chromosome'],
        'gene': ['pseudogene', 'transposable_element_gene'],
        'transcript': ['mRNA', 'snRNA', 'rRNA', 'snoRNA', 'mRNA_TE_gene',
                       'miRNA', 'tRNA', 'ncRNA', 'pseudogenic_transcript'],
        'exon': ['pseudogenic_exon'],
    }

    def reference(self, record):
        ref = annotation.Reference(record.seqid, record.end)
        for child in self.handle(record.children):
            child.reference = ref

        return ref

    def transcript(self, record):
        t = Transcript(record.attributes['ID'])
        for child in self.handle(record.children):
            if isinstance(child, annotation.Exon):
                child.transcript = t
        return t


build_annotation = AnnotationBuilder()


if __name__ == '__main__':

    #samfile = pysam.Samfile('alignments/bowtie2_alignments.sam')
    samfile = pysam.Samfile('test.sam', 'r')
    alignments = read_samfile(samfile)

    logging.info('Loading annotation')

    with open('TAIR10_GFF3_genes.gff') as fh:
        records = gff.Reader(fh)
        logging.info('Building tree')
        tree = gff.GFFTree(records)
        logging.info('Transforming to annotation')
        anno = build_annotation(tree.roots)

    logging.info('Done loading annotation')

    # Filter alignments by mismatch
    MISMATCH_THRESHOLD = 2
    alignments = (al for al in alignments if al.mismatches <= MISMATCH_THRESHOLD)

    # Filter read alignments that aligned to rRNA
    #alignments = filter_read_alignments_by_reference(alignments, rRNA_IDs)

    # Find which gene(s) each alignment corresponds to
    gene_alignment_finder = GeneAlignmentFinder(anno.genes, anno.transcripts)
    gene_alignments = gene_alignment_finder(alignments)
        
    logging.info('Making gene counts')

    gene_counts = count_genes(gene_alignments)
    for gene, counts in gene_counts.items():
        print gene, counts
