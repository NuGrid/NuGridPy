"""NuGridPy package version"""

import subprocess
__version__ = '0.8.0'
try:
    commit = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()
    __version__ += f"+{commit}"
except Exception:
    pass