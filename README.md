

## flu-evidence

This repository contains the code used for flu data analysis first published in

Baguelin, Marc, Stefan Flasche, Anton Camacho, Nikolaos Demiris, Elizabeth Miller, and W. John Edmunds. ‘Assessing Optimal Target Populations for Influenza Vaccination Programmes: An Evidence Synthesis and Modelling Study’. PLoS Med 10, no. 10 (8 October 2013): e1001527. doi:10.1371/journal.pmed.1001527.

The code is distributed under a GPLv3 license. Please remember to cite the above manuscript when using this code.


## INSTALL

### Dependencies

The code depends on
- GSL
- Boost
- Mongoclient-cxx (for JSON parsing)

The first two should be readily available with any package manager. Mongoclient can be installed as follows:
```
git clone https://github.com/mongodb/mongo-cxx-driver.wiki.git
cd mongo-cxx-driver
scons build --prefix=/usr/local/
sudo scons install --prefix=/usr/local/
```

### Compile

To successfully install this you need cmake, gsl and boost installed. Then run:

```
cmake .
make
```

and it will build the needed executables in the bin directory:

### Documentation

If you have doxygen installed you can run:

```
make doxygen
```

Which will create documentation under doc/html/index.html
