

## flu-evidence

This repository contains the code used for flu data analysis first published in

Baguelin, Marc, Stefan Flasche, Anton Camacho, Nikolaos Demiris, Elizabeth Miller, and W. John Edmunds. ‘Assessing Optimal Target Populations for Influenza Vaccination Programmes: An Evidence Synthesis and Modelling Study’. PLoS Med 10, no. 10 (8 October 2013): e1001527. doi:10.1371/journal.pmed.1001527.

The code is distributed under a GPLv3 license. Please remember to cite the above manuscript when using this code.

## Install

### Dependencies

The code depends on
- cmake
- GSL
- Boost
- Mongoclient-cxx (for JSON parsing)
- Eigen3 (for matrix mathematics/eigenvectors)

The first three should be readily available with any package manager. Mongoclient can be installed as follows (depends on scons):

```
git clone https://github.com/mongodb/mongo-cxx-driver.git
cd mongo-cxx-driver
scons build --prefix=/usr/local/
sudo scons install --prefix=/usr/local/
```

For eigen3, download the latest stable release from: http://eigen.tuxfamily.org/index.php?title=Main_Page#Download
Then unpack it (tar xvf), cd into the just created directory (eigen-eigen-randomno) and do:

```
mkdir build
cd build
cmake ..
make
sudo make install
```

### Compile

To successfully install this you need the dependencies installed. Then run:

```
cmake .
make
```

and it will build the needed executables in the bin directory. Functionality is divided into multiple binaries. Currently bin/flu-evidence-synthesis runs the MCMC code. And bin/inference can be used to run the different vaccination scenarios using the results from the MCMC. 

At the moment the data we used is not yet included in this repo. Feel free to contact us for further help/details.

## Documentation

If you have doxygen installed you can run:

```
make doxygen
```

Which will create documentation under doc/html/index.html
