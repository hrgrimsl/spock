import sys
import os

def test_import():
    sys.path.append('..')
    os.system('pwd')
    from spock import core

