This is pysplib - a Python library for spectral methods.
=========================================================


--------------
Fortran code on which it is based is written by Daniele Funaro in 1994.
You may find the original source code at:  
cdm.unimo.it/home/matematica/funaro.daniele/splib.txt


I have removed the calls to FFT subroutines from the NAG library  
and replaced them by calls to FFTW3 (http://www.fftw.org/).   
    
This was all done in order to create a Python extension module.   
It was achieved in a following way:  
  
  
To create a signature file: 
``` 
f2py -m pysplib -h splib.fpy splib.f fftw3.f
```  

And second step, to create extension module:  
```
f2py -c -m --fcompiler=gnu95 pysplib splib.f -lfftw3  
```
  
I'm using gfortran, therefore I use --fcompiler=gnu95 switch.  


After doing this you get pysplib.so file. Copy it to ```examples``` subdirectory:  
```
cp pysplib.so examples/
```  
  
And finally run examples:  
```
python ex1.py
```   
  
You may use pysplib by importing it in Python:  

```python
  import pysplib
  pysplib.gammaf(15)
  87178291200.0
...
```
  
Have fun!  
  

# Licence:    
Available under GNU General Public License version 3  
  
  
--------------
For question write to:  
Nikola Mirkov   
email:  
nmirkov@vinca.rs  
largeddysimulation@gmail.com  
  
[![githalytics.com
alpha](https://cruel-carlota.pagodabox.com/3523d5840db5d0347b41d8003e5ceabf
"githalytics.com")](http://githalytics.com/nikola-m/pysplib)  

