import matplotlib.image as mpimg
import matplotlib.pylab as plb
import numpy
import sys
from scipy import stats
import glob

def compare_entropy(name_img1,name_img2,method="KL-div"):
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
     elif method == "rms":
          fdiff=fimg1-fimg2
          fdiff_sqr = fdiff**4
          S = (fdiff_sqr.sum())**(1./4)

     return S,fimg1, fimg2

S_limit = -3.
file_list = glob.glob('Abu*')
S=[]
print("Identifying images with modified log KL-divergence > "+'%3.1f'%S_limit)
for i in range(len(file_list)):
    this_S,fimg1,fimg2 = compare_entropy(file_list[i-1],file_list[i])
    if this_S > S_limit:
         print(file_list[i-1]+" and "+file_list[i]+" differ by "+'%6.3f'%this_S)
    S.append(this_S)
    
plb.plot(S,'o')
plb.xlabel("image number")
plb.ylabel("modified log KL-divergence to previous image")
plb.show()
