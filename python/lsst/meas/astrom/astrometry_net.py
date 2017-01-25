from ._astrometry_net import MultiIndex, Solver, healpixDistance, finalize

__all__ = ["MultiIndex", "Solver", "healpixDistance", "finalize"]


def __iter__(self):
    """Get an iterator over the indices in a MultiIndex

    Do not modify the number or location of the indices while using the iterator.
    """
    for i in range(len(self)):
        yield self[i]
MultiIndex.__iter__ = __iter__
