import multiprocessing


class Worker(multiprocessing.Process):
    """A worker."""

    def __init__(self, inqueue, outqueue, to_run):
        multiprocessing.Process.__init__(self)
        self.q_in = inqueue
        self.q_out = outqueue
        self.to_run = to_run

    def run(self):
        """."""
        # print("Starting %s" % self.name)
        while True:
            task = self.q_in.get()
            if task is None:
                # print("Ending %s" % self.name)
                # self.q_in.task_done()
                break

            for out in self.to_run(*task):
                self.q_out.put(out)
            # signal that we did our job
            # self.q_in.task_done()
        return