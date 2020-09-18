def forward_match(seq: str, recog: list) -> (bool):
    """
    Match the sequence with a recognition sequence.

    :param seq: the sequence to search in
    :param recogn: a list with bases that should be
                   matched subsequently
    :return: True if the sequence matches the
             recognition sequence, False if not
    """
    if len(recog) > len(seq):
        return False
    for idx, bases in enumerate(recog):
        if bases and not seq[idx] in bases:
            return False
    return True


def reverse_match(seq: str, recog: list) -> (bool):
    """
    Match the sequence with a recognition sequence in reverse.

    :param seq: the sequence to search in
    :param recogn: a list with bases that should be
                   matched subsequently
    :return: True if the sequence matches the
             recognition sequence, False if not
    """
    if len(recog) > len(seq):
        return False
    offset = len(seq) - len(recog)
    for idx, bases in enumerate(recog):
        if bases and not seq[offset + idx] in bases:
            return False
    return True


def match(seq: str, recog: list, from_begin: bool=True) -> (bool):
    """
    Match the sequence with a recognition sequence.

    :param seq: the sequence to search in
    :param recog: a list with bases that should be
                        matched subsequently
    :param from_begin: should we match the sequence from the
                       begin or the end
    :return: True if the sequence matches the
             recognition sequence, False if not
    """
    if from_begin:
        return forward_match(seq, recog)
    else:
        return reverse_match(seq, recog)


def sites(sequence: str, recog: list):
    """
    Find the occurences of a recognition sequence in a larger sequence.

    :param seq: the sequence to search in
    :param recog: a list with bases that should be matched subsequently
    :yield: the start of the recognized sequence
    """
    def submatch(offset):
        """match from start onward."""
        for idx, bases in enumerate(recog):
            if not bases:
                continue
            crd = offset + idx
            if crd >= len(sequence):
                break
            if not sequence[crd] in bases:
                return False
        return True

    for start in range(len(sequence)):
        if submatch(start):
            yield start
