from __future__ import absolute_import

try:
    from tempfile import TemporaryDirectory
except:
    from ..py2 import TemporaryDirectory

from .abu_chart import load_chart_files
from .compare_image_entropy import compare_images

def test():
    with TemporaryDirectory() as tdir:
        load_chart_files(tdir)
        compare_images(tdir)

if __name__ == "__main__":
    test()
