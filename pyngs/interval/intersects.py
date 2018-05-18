"""
Objects to represent Intersects.

This module holds objects to parse overlaps
determined using bedtools intersect.
"""


class Intersects(object):
    """Parse overlaps determined by bedtools intersect."""

    def __init__(self):
        """Initialize the object."""
        pass

    def __next__(self):
        """Get the next overlap."""
        pass

    def __iter__(self):
        """Mark this object as an iterator."""
        return self
