import subprocess as sp
import os

if __name__ == '__main__':
  fwritelog = open('run.log', 'w')
  os.chdir('/home/icl/genray-c/test_ech')
  sp.call(['/home/icl/genray-c/genray-c_160826.1/xgenray'], shell = False, stdout = fwritelog)
