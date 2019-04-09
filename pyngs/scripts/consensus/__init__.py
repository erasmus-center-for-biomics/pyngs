import sys
import gzip
import itertools
import multiprocessing
from pyngs import sam

from .header import consensus_header
from .task import PairedConsensus
from .worker import Worker
from .writer import Writer


def group_per_umi(reader, tag="um"):
    """Get alignments per UMI from a UMI sorted samparser."""
    batch = []
    umi = None
    for alignment in reader:

        # get the UMI for the alignment
        curumi = alignment.get_tag(tag)
        if curumi is None:
            yield "", [alignment]
            continue

        # if the UMI switches more
        if umi != curumi[2]:
            if batch:
                yield umi, batch
            umi = curumi[2]
            batch = []

        # add the alignment to the buffer
        batch.append(alignment)

    # yield the last umi with its alignments
    if batch:
        yield umi, batch


def make_consensus(inpath, outpath, tag="um", max_distance=20, discard=False, nworkers=8):
    """Run make consensus alignments."""
    # prepare the queues
    to_workers = multiprocessing.JoinableQueue()
    to_writer = multiprocessing.JoinableQueue()

    # prepare reading the input file
    instream = sys.stdin
    if inpath != "stdin":
        if inpath.endswith(".gz"):
            instream = gzip.open(inpath, "rt")
        else:
            instream = open(inpath, "rt")

    reader = sam.Reader(instream)
    header, rgid = consensus_header(reader.header, reader.readgroups())

    # prepare the read processor
    conobj = PairedConsensus()
    conobj.discard = discard
    conobj.consensus_id = rgid
    conobj.max_distance = max_distance

    # initialize the writer
    writer = Writer(to_writer, outpath, header)
    writer.start()

    # initialize the workers
    workers = []
    for _ in range(nworkers):
        workers.append(Worker(to_workers, to_writer, conobj))
    for w in workers:
        w.start()

    # feed the data to the workers
    for umi, alignments in group_per_umi(reader, tag):
        to_workers.put((umi, alignments))

    # add the poison pills at the end of the stack and join the workers
    for _ in range(nworkers):
        to_workers.put(None)
    to_workers.join()
    print("Joined worker queue")
    # add the poison pill at the end of the writer and join it
    to_writer.put(None)

    print("Joining writer")
    to_writer.join()
    writer.join()
    print("Joined writer")

    # joining workers
    for w in workers:
        w.join()
    print("joined workers")
    # close the input file
    if not instream.closed and instream != sys.stdin:
        instream.close()
    print("returning")
    return


def run_consensus(args):
    """Run the consensus calling."""
    # open the files
    make_consensus(
        args.sam,
        args.out,
        args.tag,
        args.distance,
        args.discard,
        nworkers=args.workers)
