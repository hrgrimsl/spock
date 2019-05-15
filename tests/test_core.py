import sys


def test_import():
    sys.path.append('../..')
    import spock
    python spock.core.py
