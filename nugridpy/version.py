"""NuGridPy package version"""

import subprocess
__version__ = '0.7.6'
try:
    commit = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()
    __version__ += f"+{commit}"
except Exception:
    pass