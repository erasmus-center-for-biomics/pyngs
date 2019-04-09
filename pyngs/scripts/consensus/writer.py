import gzip
import sys
import multiprocessing
from pyngs import sam


class Writer(multiprocessing.Process):
    """A worker."""

    def __init__(self, inqueue, filename, header):
        multiprocessing.Process.__init__(self)
        self.q_in = inqueue
        self.filename = filename
        self.header = header

    def run(self):
        # prepare the output stream
        print("Starting %s" % self.name)
        outstream = sys.stdout
        if self.filename != "stdout":
            if self.filename.endswith(".gz"):
                outstream = gzip.open(self.filename, "wt")
            else:
                outstream = open(self.filename, "wt")

        # prepare the sam writer
        writer = sam.Writer(outstream, self.header)

        while True:
            task = self.q_in.get()
            if task is None:
                print("Ending %s" % self.name)
                break
            for aln in task:
                writer.write(aln)


        #
        if outstream is not sys.stdout:
            outstream.close()
