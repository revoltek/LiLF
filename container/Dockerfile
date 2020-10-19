FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
LABEL maintainer="Francesco de Gasperin"
LABEL version="0.1"

USER root

# Support large mlocks
RUN echo "*        -   memlock     unlimited" > /etc/security/limits.conf

RUN apt-get update
RUN apt-get install -y g++ gfortran cmake git wget flex bison vim \
    python3-dev python3-ipython python3-pip python3-setuptools \
    libcfitsio-dev libpng-dev libxml2-dev libarmadillo-dev liblua5.3-dev libfftw3-dev python3-pybind11 wcslib-dev libgsl-dev \
    libboost-dev libboost-date-time-dev libboost-filesystem-dev libboost-system-dev libboost-test-dev libboost-python-dev libboost-numpy-dev libboost-program-options-dev

RUN pip3 install pyregion
RUN pip3 install astropy --upgrade

# Add group. We chose GID 65527 to try avoiding conflicts.
RUN groupadd -g 65527 lofar

# Add user. We chose UID 65527 to try avoiding conflicts.
RUN useradd lofar -d /home/lofar -u 65527 -g 65527 -m -s /bin/bash

# Add metuaser user to sudoers
RUN adduser lofar sudo

# No pass sudo (for everyone, actually)
COPY sudoers /etc/sudoers

# Setup home directories
RUN chown lofar:lofar /home/lofar

####################################################################
## CASACORE
###################################################################
RUN cd /tmp && git clone https://github.com/casacore/casacore.git
RUN cd /tmp/casacore && git checkout tags/v3.3.0
RUN cd /tmp/casacore && mkdir data && cd data && wget --retry-connrefused ftp://anonymous@ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar && tar xf WSRT_Measures.ztar
RUN cd /tmp/casacore && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DDATA_DIR=${INSTALLDIR}/casacore/data -DBUILD_PYTHON=True -DUSE_OPENMP=True -DUSE_FFTW3=TRUE -DUSE_HDF5=True -DBUILD_PYTHON:BOOL=False -DBUILD_PYTHON3:BOOL=True .. && make -j 64 && make install

#####################################################################
## CASACORE-python
#####################################################################
RUN cd /tmp && git clone https://github.com/casacore/python-casacore.git
RUN cd /tmp/python-casacore && git checkout tags/v3.3.1
RUN cd /tmp/python-casacore && pip3 install .

#####################################################################
## LOFARBeam
#####################################################################
RUN cd /tmp && git clone https://github.com/lofar-astron/LOFARBeam.git
RUN cd /tmp/LOFARBeam && git checkout 5b5d80ff308784adbf00985293257b16afff86b2
RUN cd /tmp/LOFARBeam && mkdir build && cd build \
    && cmake .. \
    && make -j 64 && make install
ENV PYTHONPATH /usr/local/lib/python3.8/site-packages/:$PYTHONPATH
ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH

#####################################################################
## idg
#####################################################################
RUN cd /tmp &&  git clone https://gitlab.com/astron-idg/idg.git
RUN cd /tmp/idg && git checkout f4020d47a9e20f85d9e79356ca98d3bdc9ef6ade
RUN cd /tmp/idg && mkdir build && cd build \
    && cmake -DCMAKE_INSTALL_PREFIX:PATH=/ .. \
    && make -j 64 && make install

#####################################################################
## AOFlagger v2.15.0
#####################################################################
#RUN cd /tmp && git clone https://gitlab.com/aroffringa/aoflagger.git
RUN cd /tmp && git clone git://git.code.sf.net/p/aoflagger/code aoflagger
RUN cd /tmp/aoflagger && git checkout 58c25d
RUN cd /tmp/aoflagger && mkdir build && cd build && mkdir install \
    && cmake .. \
    && make -j 64 && make install

#####################################################################
## Dysco
#####################################################################
RUN cd /tmp && git clone https://github.com/aroffringa/dysco.git
RUN cd /tmp/dysco && git checkout 7ef0f4716b6772b78b39b09d72821a84f6471a21
RUN cd /tmp/dysco && mkdir build && cd build \
    && cmake .. \
    && make -j 64 && make install

#####################################################################
## Wsclean
#####################################################################
RUN cd /tmp && git clone https://gitlab.com/aroffringa/wsclean.git
RUN cd /tmp/wsclean && git checkout 6567e489
RUN cd /tmp/wsclean && mkdir build && cd build \
    && cmake .. \
    && make -j 64 && make install

#####################################################################
## Dp3
#####################################################################
RUN cd /tmp && git clone https://github.com/lofar-astron/DP3.git
RUN cd /tmp/DP3 && git checkout tags/v4.2
RUN cd /tmp/DP3 && mkdir build && cd build \
    && cmake .. \
    && make -j 64 && make install

#####################################################################
## Pybdsf
#####################################################################
RUN cd /tmp && git clone https://github.com/lofar-astron/PyBDSF.git
RUN cd /tmp/PyBDSF; git checkout 99573ee3799c8ce657d959c3bb2fdf6a61623b2c
RUN cd /tmp/PyBDSF; pip3 install .

#####################################################################
## LiLF
#####################################################################
RUN cd /opt && git clone https://github.com/revoltek/LiLF.git
RUN cd /opt/LiLF && git checkout LBAdevel
ENV PYTHONPATH /opt/LiLF:$PYTHONPATH
ENV PATH /opt/LiLF/scripts:$PATH

#####################################################################
## LSMtool
#####################################################################
RUN cd /opt && git clone https://github.com/darafferty/LSMTool.git
RUN cd /opt/LSMTool && git checkout 8730aaffd455cdd3071c68a1d56299ac1bece90c
RUN cd /opt/LSMTool; pip3 install .

####################################################################
## LoSoTo
###################################################################
RUN cd /opt && git clone https://github.com/revoltek/losoto.git
RUN cd /opt/losoto && git checkout 921f807fa5220bf5c2daeb5c68c3e1e0e66f9c66
RUN cd /opt/losoto; pip3 install .
#ENV PYTHONPATH /opt/losoto:$PYTHONPATH

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
RUN python3 -m pip install tensorflow-cpu

RUN python3 -m pip install -U "/opt/DDFacet/[moresane-support,testing-requirements,fits-beam-support]"
RUN cd /opt/DDFacet/ && python3 setup.py build

# Set MeqTrees Cattery path to installation directory
ENV MEQTREES_CATTERY_PATH /usr/local/lib/python3.6/dist-packages/Cattery/
ENV PYTHONPATH $MEQTREES_CATTERY_PATH:$PYTHONPATH
