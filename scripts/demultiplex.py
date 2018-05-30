
import sys
import argparse
import gzip
from operator import itemgetter
from itertools import groupby
import pyngs.sequence


class AssignBarcode(object):
    """."""

    def __init__(self, read=0, start=0, end=0, remove=True):
        """Initialize a new assigner."""
        self.read = read
        self.start = start
        self.end = end
        self.remove = remove
        self.buffer = []
        self.buffersize = 1000000

    def __call__(self, files, prefixes, suffixes):
        """Process FastQ files."""
        assert len(files) > self.read
        assert len(files) == len(prefixes)
        assert len(files) == len(suffixes)

        # open the files
        streams = []
        for path in files:
            if path.endswith(".gz"):
                streams.append(gzip.open(path, "rt"))
            else:
                streams.append(open(path, "r"))

        # make generators from the streams
        fastqs = [pyngs.sequence.FastQ(s) for s in streams]
        record = 0
        while True:
            record += 1
            stop = False
            # get the next set of reads
            data = []
            for fgen in fastqs:
                try:
                    data.append(next(fgen))
                except StopIteration:
                    stop = True
                    break

            # if there are no more reads to process break
            if stop:
                break

            # get the barcode
            bcread = data[self.read]
            bcseq = bcread[1][self.start:self.end]

            # remove the barcode
            if self.remove:
                read = list(data[self.read])
                read[1] = read[1][self.end:]
                read[2] = read[2][self.end:]
                data[self.read] = tuple(read)

            # add the reads to the output
            for idx in range(len(data)):
                outfile = prefixes[idx] + bcseq + suffixes[idx]
                self.buffer.append((outfile, record, data[idx]))

            # add data to the buffer
            if len(self.buffer) > self.buffersize:
                self.process()

        # process the remainder
        self.process()

        # close the input streams
        for stream in streams:
            stream.close()

    def process(self):
        """Process the reads in the buffer."""
        self.buffer.sort(key=itemgetter(0, 1))

        # write the entries to the output files
        for outfile, data in groupby(self.buffer, key=itemgetter(0)):
            handle = gzip.open(outfile, "at") if outfile.endswith(".gz") else open(outfile, "at")
            for value in data:
                read = value[2]
                pyngs.sequence.write_fastq(
                    handle,
                    read[0],
                    read[1],
                    read[2])
            handle.close()

        # clear the buffer
        self.buffer.clear()


if __name__ == "__main__":

    def main():
        """Run the main program loop."""
        parser = argparse.ArgumentParser(
            prog=sys.argv[0],
            description="""
            Split datasets based on an inline barcode.
            """
        )
        parser.add_argument(
            "-f", "--file", dest="file",
            type=str, nargs="+",
            help="The FastQ files with the reads.")
        parser.add_argument(
            "--suffix", dest="suffix",
            type=str, nargs="*",
            help="The suffix for the output FastQ files.")
        parser.add_argument(
            "-p", "--prefix", dest="prefix",
            type=str, nargs="*", default="",
            help="The prefix of the output files.")
        parser.add_argument(
            "-r", "--read", dest="read",
            type=int, default=0,
            help="The read with the barcodes.")
        parser.add_argument(
            "-s", "--start", dest="start",
            type=int, default=0,
            help="The start base of the barcode in the read.")
        parser.add_argument(
            "-e", "--end", dest="end",
            type=int, default=12,
            help="The end base of the barcode in the read.")
        parser.add_argument(
            "--remove", dest="remove_code",
            type=bool,
            help="Remove the barcode.")
        args = parser.parse_args()

        # parse the suffixes and prefixes
        suffixes = []
        if not args.suffix:
            suffixes = ["" for _ in args.file]
        elif len(args.suffix) == 1:
            suffixes = [args.suffix[0] for _ in args.file]
        elif len(args.suffix) == len(args.file):
            suffixes = args.suffix

        prefixes = []
        if not args.prefix:
            prefixes = ["" for _ in args.file]
        elif len(args.prefix) == 1:
            prefixes = [args.prefix[0] for _ in args.file]
        elif len(args.prefix) == len(args.file):
            prefixes = args.prefix
        # assign the reads
        assigner = AssignBarcode(
            read=args.read,
            start=args.start,
            end=args.end,
            remove=args.remove_code)
        assigner(args.file, prefixes, suffixes)

    # run the main program
    main()
