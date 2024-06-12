import lsst.gauss2d.fit as g2f


def test_bindings():
    errors = []
    for name_class in dir(g2f):
        klass = getattr(g2f, name_class)

        if callable(klass):
            try:
                obj = klass()
            except ValueError as e:
                # Likely non-trivial __init__
                continue
            except TypeError as e:
                msg = str(e)
                if not (msg.endswith("No constructor defined!")
                        or msg.startswith("__init__(): incompatible constructor arguments.")):
                    errors.append(e)
            # str and repr should not raise
            try:
                repr(obj)
            except Exception as e:
                errors.append(e)
            try:
                str(obj)
            except Exception as e:
                errors.append(e)

    if errors:
        print('\n'.join((str(e) for e in errors)))
        assert not 'errors not empty; see stdout'
