"""
Data
====

This module configures the application.
"""

import logging


# Logger
level = logging.INFO
logger = logging.getLogger('nugridpy')
logger.setLevel(level)

sh = logging.StreamHandler()
sh.setLevel(level)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)

# Cache maximum size
CACHE_MAX_SIZE = 64
