
def rle(noncompressed=""):
    """Run-lenght encode a list."""
    retval = []
    pchar = None
    cnt = 0
    for char in noncompressed:
        if pchar is not None and char != pchar:
            retval.append((cnt, pchar))
            cnt = 0
        cnt += 1
        pchar = char
    retval.append((cnt, pchar))
    return retval
