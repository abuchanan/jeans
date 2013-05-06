from collections import Counter

from nose.tools import eq_, ok_
import pysam

import genecounter as gc


def test_pysam_positions():
    '''
    It's important to get positions right, and in the biology
    world it's constantly confusing; some people use 1-based closed
    intervals, others use 0-based half-open.

    These tests figure out what pysam does.
    '''

    samfile = pysam.Samfile('test.sam', 'r')
    a = samfile.next()
    # In the file, the position column is "2784862"
    # the cigar column is "85M"

    # Note that this is 1 less than the value in the SAM file.
    # Pysam apparently subtracts 1 to use a 0-based interval.
    eq_(a.pos, 2784861)
    # This is the position above plus the match length (85),
    # so this is a 0-based, half-open interval.
    eq_(a.aend, 2784946)


def test_PairedAlignment():
    a = gc.Alignment('refid', 25, 50, 'readid', 3, True)
    b = gc.Alignment('refid', 70, 95, 'readid', 1, False)
    c = gc.PairedAlignment(a, b)

    eq_(c.start, 25)
    eq_(c.end, 95)
    eq_(c.reference_ID, 'refid')
    eq_(c.read_ID, 'readid')
    eq_(c.mismatches, 4)

    # note the positions are backwards now...
    d = gc.Alignment('refid', 70, 95, 'readid', 1, False)
    e = gc.Alignment('refid', 25, 50, 'readid', 3, True)
    f = gc.PairedAlignment(d, e)

    # ...but the start and end are the same
    eq_(f.start, 25)
    eq_(f.end, 95)


def test_read_sam():
    samfile = pysam.Samfile('test.sam', 'r')
    alignments = gc.read_samfile(samfile)

    alignments = list(alignments)
    eq_(len(alignments), 22)

    a = alignments[0]
    eq_(a.reference_ID, 'gi|30407130|emb|AL954747.1|')
    eq_(a.read_ID, 'DB775P1:240:D1TE2ACXX:4:1216:12799:92206')
    eq_(a.strand, '+')

    # See test_pysam_position() above. Pysam uses a 0-based, half-open
    # interval. I want to transfer back to a 1-based, closed interval.
    eq_(a.start, 2784862)
    eq_(a.end, 2784946)
    eq_(a.mismatches, 6)


def test_read_sam_paired():
    samfile = pysam.Samfile('test.sam', 'r')
    alignments = gc.read_samfile(samfile, paired=True)

    alignments = list(alignments)
    eq_(len(alignments), 11)

    a = alignments[0]
    eq_(a.reference_ID, 'gi|30407130|emb|AL954747.1|')
    eq_(a.read_ID, 'DB775P1:240:D1TE2ACXX:4:1216:12799:92206')

    eq_(a.start, 2784862)
    eq_(a.end, 2785034)
    eq_(a.mismatches, 11)


def test_filter_alignments_by_reference():

    class Alignment(object):
        def __init__(self, reference_ID, read_ID):
            self.reference_ID = reference_ID
            self.read_ID = read_ID


    a = Alignment('included', 'read1')
    b = Alignment('excluded_one', 'read1')

    c = Alignment('included', 'read2')
    d = Alignment('included', 'read2')

    e = Alignment('included', 'read3')
    f = Alignment('excluded_two', 'read3')

    alignments = [a, b, c, d, e, f]
    exclude = ['excluded_one', 'excluded_two']

    res = gc.filter_read_alignments_by_reference(alignments, exclude)
    res = list(res)

    eq_(res, [c, d])


def test_group_alignments_by_read():

    class Alignment(object):
        def __init__(self, read_ID):
            self.read_ID = read_ID

    a = Alignment('read1')
    b = Alignment('read1')
    c = Alignment('read2')
    d = Alignment('read2')
    e = Alignment('read2')
    f = Alignment('read3')

    alignments = [a, b, c, d, e, f]

    res = gc.group_alignments_by_read(alignments)
    res = [list(group) for group in res]

    eq_(res, [[a, b], [c, d, e], [f]])


def test_GeneAlignmentFinder():

    class Region(object):
        def __init__(self, reference_ID, start, end, strand):
            self.reference_ID = reference_ID
            self.start = start
            self.end = end
            self.strand = strand

        def __repr__(self):
            return 'Region({}, {}, {}, {})'.format(self.reference_ID, self.start, self.end, self.strand)

    class Transcript(object):
        def __init__(self, ID, start, end, strand, gene):
            self.ID = ID
            self.start = start
            self.end = end
            self.strand = strand
            self.gene = gene


    genes = [
        Region('ref_1', 10, 30, '+'),
        Region('ref_2', 10, 30, '+'),
        Region('ref_2', 20, 30, '+'),
        Region('ref_3', 20, 30, '-'),
        Region('ref_3', 20, 30, '+'),
    ]

    transcripts = [
        Transcript('transcript_1', 30, 40, '+', genes[0]),
    ]

    alignments = [
        Region('ref_1', 10, 20, '+'),
        Region('ref_2', 25, 40, '+'),

        # start, end, and strand don't really matter for transcripts
        Region('transcript_1', -1, -1, 'dont care'),

        Region('ref_3', 25, 40, '-'),
    ]

    finder = gc.GeneAlignmentFinder(genes, transcripts)
    res = finder(alignments)
    res = list(res)

    expected = [
        gc.GeneAlignment(genes[0], alignments[0]),
        gc.GeneAlignment(genes[1], alignments[1]),
        gc.GeneAlignment(genes[2], alignments[1]),
        gc.GeneAlignment(genes[0], alignments[2]),
        gc.GeneAlignment(genes[3], alignments[3]),
    ]

    eq_(expected, res)


def test_counts_genes():

    class Alignment(object):
        def __init__(self, read_ID):
            self.read_ID = read_ID

    class GeneAlignment(object):
        def __init__(self, gene, read_ID):
            self.gene = gene
            self.alignment = Alignment(read_ID)


    gene_alignments = [
        GeneAlignment('gene_1', 'read_1'),

        # won't be counted because this read aligns to multiple genes
        GeneAlignment('gene_1', 'read_2'),
        GeneAlignment('gene_2', 'read_2'),

        GeneAlignment('gene_2', 'read_3'),
    ]

    res = gc.count_genes(gene_alignments)
    expected = Counter({'gene_1': 1, 'gene_2': 1})

    eq_(expected, res)
