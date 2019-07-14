To compile and locally install PyESP, run:
> python3 setup.py build_ext --inplace

The module can then be used in Python 3 using:
> import PyESP

Description of python code use with comments is given in all_examples.py

Parameters are the same as used in the original ESP package (cf. http://lem.ch.unito.it/chemistry/esp_manual.html)

Working examples:
- EC.py : EC mechanism
- ECcat.py : E C_cat mechanism with increasing amounts of substrate for the reaction A + e <> B, then B + substrate -> A + product
- square_scheme.py : square scheme mechanism

Developed and tested on Ubuntu using Python 3.7, but should run on any modern platform (Python 3 only)