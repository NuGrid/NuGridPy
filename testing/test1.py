from __future__ import absolute_import

import tempfile
from .abu_chart import load_chart_files
from .compare_image_entropy import compare_images

def test():
    with tempfile.TemporaryDirectory() as tdir:
        load_chart_files(tdir)
        compare_images(tdir)


if __name__ == "__main__":
    test()
