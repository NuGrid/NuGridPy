from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
import matplotlib.image as mpimg
import matplotlib.pylab as plb
import numpy
import sys
from scipy import stats
import glob
import os.path
import warnings
import time

def compare_entropy(name_img1,name_img2,method="rmq"):
     '''Compare two images by the Kullback-Leibler divergence

     Parameters
     ----------
     name_img1 : string
       filename of image 1 (png format)

     name_img2 : string
       filename of image 2 (png format)

     Returns
     -------
     S : float
        Kullback-Leibler divergence S = sum(pk * log(pk / qk), axis=0)

     Note
     ----
     See http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html
     '''
     img1 = mpimg.imread(name_img1)
     img2 = mpimg.imread(name_img2)
     fimg1 = img1.flatten()
     fimg2 = img2.flatten()
     if method == "KL-div":
          eps = 0.0001
          S = stats.entropy(fimg2+eps,fimg1+eps)
          S = numpy.log10(S)
     elif method == "rmq":
          fdiff=fimg1-fimg2
          fdiff_sqr = fdiff**4
          S = (fdiff_sqr.sum())**(old_div(1.,4))

     return S,fimg1, fimg2


def compare_images(path = '.'):
     S_limit = 10.
     file_list = glob.glob(os.path.join(path, 'Abu*'))
     file_list_master = glob.glob(os.path.join(path, 'MasterAbu*'))
     file_list.sort()
     file_list_master.sort()
     S=[]
     print("Identifying images with rmq > "+'%3.1f'%S_limit)
     ierr_count = 0
     for i in range(len(file_list)):
         this_S,fimg1,fimg2 = compare_entropy(file_list[i],file_list_master[i])
         if this_S > S_limit:
              warnings.warn(file_list[i]+" and "+file_list_master[i]+" differ by "+'%6.3f'%this_S)
              ierr_count += 1
              S.append(this_S)
     if ierr_count > 0:
          print("Error: at least one image differs by more than S_limit")
          sys.exit(1)
     #print ("S: ",S)
     #plb.plot(S,'o')
     #plb.xlabel("image number")
     #plb.ylabel("modified log KL-divergence to previous image")
     #plb.show()


if __name__ == "__main__":
     compare_images()
