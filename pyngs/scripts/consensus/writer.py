import gzip
import sys
import multiprocessing
from pyngs import sam


class SAMWriter(multiprocessing.Process):
    """A worker."""

    def __init__(self, inqueue, filename, header):
        multiprocessing.Process.__init__(self)
        self.q_in = inqueue
        self.filename = filename
        self.header = header

    def run(self):
        # prepare the output stream
        outstream = sys.stdout
        if self.filename != "stdout":
            if self.filename.endswith(".gz"):
                outstream = gzip.open(self.filename, "wt")
            else:
                outstream = open(self.filename, "wt")
        writer = sam.Writer(outstream, self.header)

        while True:
            task = self.q_in.get()
            if task is None:
                self.q_in.task_done()
                # print("Ending %s" % self.name)
                break
            for aln in task:
                writer.write(aln)
            self.q_in.task_done()

        if outstream is not sys.stdout:
            outstream.close()
