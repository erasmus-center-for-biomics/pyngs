

def quote_tokenizer(sstring: str="", sep: str=",", quote='"'):
    """Separate the fields while ignoring quotes."""
    tokens = sstring.split(sep)
    keep = False
    current = []
    for token in tokens:
        if keep:
            current.append(token)
            if token[-1] == quote:
                keep = False
                yield sep.join(current)
        elif quote in token:
            # if the quote is ended in the same token,
            # just yield the token
            if token[-1] == quote:
                yield token
            # otherwise accumulate until the ending quote.
            else:
                keep = True
                current = [token]
        # not quotes, no issue
        else:
            yield token
