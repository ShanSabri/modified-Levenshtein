# modified-Levenshtein
A direct and efficient translation of the dynamic programming algorithm of the Sequence-Levenshtein distance into Cython. 


This repository contains a Cython source and its associated Setup file, both of which are written in Python 3. Once compiled, the distance function can be imported as a typical module function.

### Installation
modified-Levenshtein requires [Cython](https://cython.org) V.0.28.2 to run.

Building from source:
```sh
$ cython -V
Cython version 0.28.2
$ python3 setup.py build_ext --inplace
```

Importing within Python:
```python
>>> from cLev import dist
>>> dist("CAGG", "CGTC")
2
>>> dist("TAGG", "TCCATGCATA")
3
```

### Development

It's been a while since I wrote this code so if you find any mistakes or optimizations please let me know! Feel free to submit a pull request! 

### References
I found the following refrences particularly helpful while hacking this up:
 - https://github.com/gfairchild/pyxDamerauLevenshtein
 - https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance
 - http://hackmap.blogspot.com/2008/04/levenshtein-in-cython.html

License
----
MIT
