import sys
import os

def test_import():
    sys.path.append('.')
    os.system('ls')
    from spock import core
    core()
    

