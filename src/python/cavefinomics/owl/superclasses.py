from itertools import chain

def get_superclasses(cls, max_depth=None):
    import owlready2
    for s in cls.is_a:
        if isinstance(s, owlready2.entity.ThingClass):
            yield s
    for v in chain.from_iterable(get_superclasses(s) for s in cls.is_a):
        yield v
