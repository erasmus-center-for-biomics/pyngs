

class Distance:

    def __init__(self, width=0, height=0):
        """Initialize a new Levenshstein distance object."""

        def default(cha, chb):
            """The default matcher."""
            if cha == "N" or chb == "N":
                return True
            elif cha == chb:
                return True
            else:
                return False

        #
        self.__set_matrix__(width, height)

        #
        self.match = 0
        self.substitute = 1
        self.insertion = 1
        self.deletion = 1
        self.equals = default

    def __set_matrix__(self, height, width):
        """Initialize the matrix."""
        assert width > 0
        assert height > 0
        self.__matrix__ = []

        for idx in range(height):
            m = [0] * width
            m[0] = idx
            self.__matrix__.append(m)
        for idx in range(width):
            self.__matrix__[0][idx] = idx

    def __call__(self, source, target):
        """Score the query versus the target."""

        # check the width and height of the
        assert len(source) >= len(self.__matrix__)
        assert len(target) >= len(self.__matrix__[0])

        # foreach column and row
        for col in range(1, len(self.__matrix__[0])):
            for row in range(1, len(self.__matrix__)):

                # cost of a match or a substitution
                subcost = self.__matrix__[row - 1][col - 1]
                if self.equals(target[col - 1], source[row - 1]):
                    subcost += self.match
                else:
                    subcost += self.substitute

                # cost of a deletion
                delcost = self.__matrix__[row - 1][col] + self.deletion

                # cost of an insertion
                inscost = self.__matrix__[row][col - 1] + self.insertion

                # set the next position in the matrix
                self.__matrix__[row][col] = min(subcost, delcost, inscost)

    @property
    def score(self):
        """Get the currently calculated distance score."""
        return self.__matrix__[-1][-1]
