FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
LABEL maintainer="Francesco de Gasperin"
LABEL version="20230417"

USER root

# Support large mlocks
RUN echo "*        -   memlock     unlimited" > /etc/security/limits.conf

RUN apt-get update;
RUN apt-get install -y g++ gfortran cmake git wget flex bison sudo vim hdf5-tools rename python3-magic \
    python3-dev ipython3 python3-pip python3-setuptools python3-astropy python3-pyregion python3-regions python3-h5py python3-sshtunnel python3-pymysql python3-requests python3-numexpr python3-numpy python3-astroquery python3-aplpy python3-cytoolz python3-shapely python3-tqdm python3-numba \
    libaio1 libcfitsio-dev libpng-dev libxml2-dev libarmadillo-dev liblua5.3-dev libfftw3-dev python3-pybind11 wcslib-dev libgsl-dev libblas-dev \
    libboost-dev libboost-date-time-dev libboost-filesystem-dev libboost-system-dev libboost-test-dev libboost-python-dev libboost-numpy-dev libboost-program-options-dev

# Some python packages
#RUN pip3 install pyregion
#RUN pip3 install astropy --upgrade

# Add group. We chose GID 65527 to try avoiding conflicts.
RUN groupadd -g 65527 lofar

# Add user. We chose UID 65527 to try avoiding conflicts.
RUN useradd lofar -d /home/lofar -u 65527 -g 65527 -m -s /bin/bash

# Setup home directories
RUN chown lofar:lofar /home/lofar

####################################################################
## CASACORE v3.5.0
###################################################################
RUN cd /opt && git clone https://github.com/casacore/casacore.git \
    && cd /opt/casacore && git checkout tags/v3.5.0
RUN cd /opt/casacore && mkdir data && cd data && wget --retry-connrefused ftp://anonymous@ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar && tar xf WSRT_Measures.ztar
RUN cd /opt/casacore && mkdir build && cd build && cmake -DPORTABLE=True -DCMAKE_BUILD_TYPE=Release -DDATA_DIR=/opt/casacore/data -DBUILD_PYTHON=True -DUSE_OPENMP=True -DUSE_FFTW3=TRUE -DUSE_HDF5=True -DBUILD_PYTHON:BOOL=False -DBUILD_PYTHON3:BOOL=True .. && make -j 64 && make install

#####################################################################
## CASACORE-python v3.5.1
#####################################################################
RUN cd /opt && git clone https://github.com/casacore/python-casacore.git \
    && cd /opt/python-casacore && git checkout tags/v3.5.1
RUN cd /opt/python-casacore && pip3 install .

#####################################################################
## LOFARBeam (master 6/2020)
#####################################################################
#RUN cd /opt && git clone https://github.com/lofar-astron/LOFARBeam.git \
#    && cd /opt/LOFARBeam && git checkout 36e56402bd701087948d2c65046a8d72e9668bbe
#RUN cd /opt/LOFARBeam && mkdir build && cd build \
#    && cmake -DPORTABLE=True .. \
#    && make -j 64 && make install
#
#ENV PYTHONPATH /usr/local/lib/python3.8/site-packages/:$PYTHONPATH
#ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH

#####################################################################
## EveryBeam v0.4.0
#####################################################################
RUN cd /opt && git clone https://git.astron.nl/RD/EveryBeam.git \
    && cd /opt/EveryBeam && git checkout tags/v0.4.0
RUN cd /opt/EveryBeam && mkdir build && cd build \
    && cmake -DDOWNLOAD_LOBES=Off .. \
    && make -j 64 && make install
    
#####################################################################
## idg v1.0.0
#####################################################################
RUN cd /opt && git clone https://gitlab.com/astron-idg/idg.git \
    && cd /opt/idg && git checkout tags/1.0.0
RUN cd /opt/idg && mkdir build && cd build \
    && cmake -DPORTABLE=True -DCMAKE_INSTALL_PREFIX:PATH=/ .. \
    && make -j 64 && make install

#####################################################################
## AOFlagger v3.2.0
#####################################################################
RUN cd /opt && git clone https://gitlab.com/aroffringa/aoflagger.git \
    && cd /opt/aoflagger && git checkout tags/v3.2.0
RUN cd /opt/aoflagger && mkdir build && cd build && mkdir install \
    && cmake -DPORTABLE=True .. \
    && make -j 64 && make install

#####################################################################
## Dysco v1.2
#####################################################################
#RUN cd /opt && git clone https://github.com/aroffringa/dysco.git \
#    && cd /opt/dysco && git checkout tags/v1.2
#RUN cd /opt/dysco && mkdir build && cd build \
#    && cmake -DPORTABLE=True .. \
#    && make -j 64 && make install

#####################################################################
## LofarStMan
#####################################################################
RUN cd /opt && git clone https://github.com/lofar-astron/LofarStMan
#   && cd /opt/LofarStMan
RUN cd /opt/LofarStMan && mkdir build && cd build \
    && cmake -DPORTABLE=True .. \
    && make -j 64 && make install

#####################################################################
## Wsclean master 19/4/23
#####################################################################
RUN cd /opt && git clone https://gitlab.com/aroffringa/wsclean.git \
    && cd /opt/wsclean && git checkout 69cba6c20dc6b6596481e4bb11e4d5ad2d0673d6
#    && cd /opt/wsclean && git checkout tags/v2.10.1
RUN cd /opt/wsclean && mkdir build && cd build \
    && cmake -DPORTABLE=True .. \
    && make -j 64 && make install

#####################################################################
## Dp3 master 21/4/23
#####################################################################
RUN cd /opt && git clone https://github.com/lofar-astron/DP3.git \
    && cd /opt/DP3 && git checkout 4fbb6120a702397d81429c3b13fabcfa5962110b
RUN ls /usr/local/share/everybeam
RUN ls /usr/local/lib
RUN cd /opt/DP3 && mkdir build && cd build \
    && cmake -DPORTABLE=True .. \
    && make -j 64 && make install

#####################################################################
## LTA
#####################################################################
RUN cd /opt && wget https://lta.lofar.eu/software/lofar_lta-2.7.1.tar.gz && tar xf lofar_lta-2.7.1.tar.gz && rm lofar_lta-2.7.1.tar.gz
RUN cd /opt/lofar_lta-2.7.1 && python3 setup.py install_oracle
RUN cd /opt/lofar_lta-2.7.1 && python3 setup.py install
ENV LD_LIBRARY_PATH /usr/lib/instantclient_11_2:$LD_LIBRARY_PATH

#####################################################################
## Pybdsf (master 1/2021)
#####################################################################
RUN cd /opt && git clone https://github.com/lofar-astron/PyBDSF.git \
    && cd /opt/PyBDSF; git checkout tags/v1.10.1
RUN cd /opt/PyBDSF && pip3 install .

#####################################################################
## LSMtool (1.4.5)
#####################################################################
RUN cd /opt && git clone https://github.com/darafferty/LSMTool.git \
    && cd /opt/LSMTool; git checkout tags/v1.4.5
RUN cd /opt/LSMTool && pip3 install .

####################################################################
## LoSoTo master 19/4/23
###################################################################
RUN cd /opt && git clone https://github.com/revoltek/losoto.git \
    && cd /opt/losoto && git checkout 74becaae02d44889b27af739b79bf4bbbe31a367
RUN cd /opt/losoto && pip3 install .
RUN pip3 install --upgrade scipy

#####################################################################
# msoverview
#####################################################################
RUN apt-get install -y subversion
RUN cd /opt && mkdir -p msoverview/src && cd msoverview/src && svn export https://svn.astron.nl/LOFAR/tags/LOFAR-Release-4_0_17/CEP/MS/src/msoverview.cc \
    && gcc -I/usr/local/include/casacore -L/usr/local/lib msoverview.cc -o /usr/local/bin/msoverview -lcasa_casa -lcasa_ms -lcasa_tables -lstdc++

#####################################################################
## BUILD DDF FROM SOURCE
#####################################################################
#RUN cd /tmp && git clone https://github.com/cyriltasse/DDFacet.git
#RUN cd /tmp/DDFacet && git checkout c0bd5f9ccaa9b7c28376e3e2674428902fa5bf85

COPY DDFacet /tmp/DDFacet
RUN mkdir /opt/DDFacet

# Copy DDFacet and SkyModel into the image
RUN cp -a /tmp/DDFacet/DDFacet /opt/DDFacet/DDFacet && \
    cp -a /tmp/DDFacet/SkyModel /opt/DDFacet/SkyModel && \
    cp -a /tmp/DDFacet/MANIFEST.in /opt/DDFacet/MANIFEST.in && \
    cp -a /tmp/DDFacet/setup.py /opt/DDFacet/setup.py && \
    cp -a /tmp/DDFacet/setup.cfg /opt/DDFacet/setup.cfg && \
    cp -a /tmp/DDFacet/README.rst /opt/DDFacet/README.rst && \
    cp -a /tmp/DDFacet/pyproject.toml /opt/DDFacet/pyproject.toml && \
    cp -a /tmp/DDFacet/.git /opt/DDFacet/.git && \
    cp -a /tmp/DDFacet/.gitignore /opt/DDFacet/.gitignore && \
    cp -a /tmp/DDFacet/.gitmodules /opt/DDFacet/.gitmodules

# Finally install DDFacet
RUN rm -rf /opt/DDFacet/DDFacet/cbuild
#RUN pip3 install -U pip setuptools wheel
#RUN python3 -m pip install tensorflow-cpu

# Force numba <0.55 otherwise it breaks for numpy
#RUN pip3 install numba==0.54.1
# Forse ipython not to be upgraded too much otherwise we lose compatibility with come astropy
#RUN pip3 install ipython==7.17.0

RUN python3 -m pip install -U "/opt/DDFacet/[moresane-support,testing-requirements,fits-beam-support]"
RUN cd /opt/DDFacet/ && python3 setup.py build
RUN cd /opt/DDFacet/ && python3 setup.py install

# Set MeqTrees Cattery path to installation directory
#ENV MEQTREES_CATTERY_PATH /usr/local/lib/python3.6/dist-packages/Cattery/
#ENV PYTHONPATH $MEQTREES_CATTERY_PATH:$PYTHONPATH

#####################################################################
## LiLF
#####################################################################
#RUN cd /opt && git clone https://github.com/revoltek/LiLF.git
#RUN cd /opt/LiLF && git checkout LBAdevel
#ENV PYTHONPATH /opt/LiLF:$PYTHONPATH
#ENV PATH /opt/LiLF/scripts:$PATH
   
#####################################################################

RUN ln -s /usr/bin/python3 /usr/bin/python
RUN ln -s /usr/bin/ipython3 /usr/bin/ipython
ENV HDF5_USE_FILE_LOCKING FALSE
ENV OMP_NUM_THREADS 1
ENV OPENBLAS_NUM_THREADS 1
RUN ulimit -n 4000
