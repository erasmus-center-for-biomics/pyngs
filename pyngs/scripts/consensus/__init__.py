import sys
import gzip
import itertools
import argparse
from multiprocessing import Queue 
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


def make_consensus(
        inpath: str, outpath: str, 
        tag: str="um", max_distance: int=20, discard: bool=False, nworkers: int=8, queue_size: int=1000):
    """Run make consensus alignments."""
    # prepare the queues
    to_workers = Queue(maxsize=queue_size)
    to_writer = Queue(maxsize=queue_size)

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

    # print("Joining workers")
    for w in workers:
        w.join()

    # add the poison pill at the end of the writer and join it
    # print("Joining writer")
    to_writer.put(None)
    writer.join()
    # print("Joined writer and workers")

    # close the input file
    if not instream.closed and instream != sys.stdin:
        instream.close()
    # print("returning")


def run_consensus(args: argparse.Namespace):
    """Run the consensus calling."""
    # open the files
    make_consensus(
        args.sam,
        args.out,
        args.tag,
        args.distance,
        args.discard,
        nworkers=args.workers,
        queue_size=args.queue_size)

if __name__ == '__main__':
    sparser = argparse.ArgumentParser(
        prog="extract_subread_from_fastq", 
        description="""Extract a subread from a FastQ file.""")
    sparser.add_argument(
        "-s", "--sam", dest="sam",
        type=str, default="stdin",
        help="""The input SAM file sorted on the tag
        with samtools sort -t {tag}.""")
    sparser.add_argument(
        "-t", "--tag", dest="tag",
        type=str, default="um",
        help="The tag-name for the UMI tag.")
    sparser.add_argument(
        "-r", "--discard", dest="discard",
        type=bool, default=False,
        help="Will non consensus alignments be discarded.")
    sparser.add_argument(
        "-d", "--distance", dest="distance",
        type=int,  default=20,
        help="""The allowed distance between the start postion of
        alignments with the same UMI.""")
    sparser.add_argument(
        "-w", "--workers", dest="workers",
        type=int, default=8,
        help="""The number of workers for the consensus alignments.""")
    sparser.add_argument(
        "-q", "--queue", dest="queue_size",
        type=int, default=100,
        help="""The number of objects to hold in the Queues.""")
    sparser.add_argument(
        "-o", "--output", dest="out",
        type=str, default="stdout",
        help="The output SAM file with the consensus sequences.")
    sparser.set_defaults(func=run_consensus)
