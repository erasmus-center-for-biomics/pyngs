import gzip
import sys
from multiprocessing import Process, Queue
from typing import List
from pyngs import sam


class Writer(Process):
    """A worker."""

    def __init__(self, inqueue: Queue, filename: str, header: List[str]):
        Process.__init__(self)
        self.q_in = inqueue
        self.filename = filename
        self.header = header

    def run(self):
        # prepare the output stream
        # print("Starting %s" % self.name)
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
                # print("Ending %s" % self.name)
                # self.q_in.task_done()
                break
            for aln in task:
                writer.write(aln)
            # self.q_in.task_done()
        #
        if outstream is not sys.stdout:
            outstream.close()
        # print("%s Finished" % self.name)
