Setting up LiLF
=================================
LiLF is generally run as a pipeline script - the majority of the processing is done by external scripts, such as `wsclean <https://gitlab.com/aroffringa/wsclean>`_ and `DP3 <https://git.astron.nl/RD/DP3>`_ . To prevent any problems caused by incorrect installations, we generally recommend to use an Apptainer container.


1a. Containers (recommended)
-----------------------------
LiLF can be run fully within the LOFAR containers either supplied by LiLF itself (see the docker scripts in the `container <https://github.com/revoltek/LiLF/tree/master/container>`_ folder on github or via the default LOFAR containers (`FLoCS <https://github.com/tikk3r/flocs/>`_)
Note that the FLoCS containers are generally built for a specific architecture, and are optimized to be compiled yourself.
                                                                                  

1b. Requirements (not recommended)
-----------------------------
Alternatively, if you enjoy debugging libraries yourself (you have been warned), you can install the requirements yourself. Note that deviating from these requirements may cause unpredictable results. Below we specify the software required for LiLF, along with the latest tested, compatible version:

* DP3 v6.5.1
* WSclean v3.6
* casacore commit a67610f21
* LSMtool v1.6.3
* pyBDSF v1.12.0
* losoto v2.5.1


2. Setting up LiLF
-------------------------
LiLF is relatively easy to set up. First of all, run

.. code-block:: bash

   ulimit -n 4000

This is important when concatenating downloaded measurement sets.

Next, change to the LiLF root directory, and add the LiLF root directory to your PYTHONPATH variable, e.g. by running

.. code-block:: bash

   cd /path/to/LilF/
   export PYTHONPATH=$PYTHONPATH:$PWD

Also, add the scripts folder to the PATH variable:

.. code-block:: bash

   cd /path/to/LiLF/scripts/
   export PATH=$PATH:$PWD

This is all that is required to set up LiLF.
