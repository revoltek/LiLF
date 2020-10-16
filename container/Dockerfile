FROM lofar/lofar-pipeline:LOFAR-Release-4_0_16

USER root

# Support large mlocks
RUN echo "*        -   memlock     unlimited" > /etc/security/limits.conf

RUN apt-get update -y
RUN apt-get install -y libcfitsio-dev libpng-dev git cmake python3-dev libxml2-dev python3-pip
RUN apt-get install -y wget libarmadillo-dev liblua5.3-dev libfftw3-dev python3-pybind11
RUN apt-get install -y libboost-dev libboost-tools-dev libboost-atomic-dev libboost-chrono-dev libboost-container-dev libboost-context-dev libboost-coroutine-dev libboost-date-time-dev libboost-exception-dev libboost-filesystem-dev libboost-iostreams-dev libboost-locale-dev libboost-log-dev libboost-math-dev libboost-numpy-dev libboost-program-options-dev libboost-python-dev libboost-random-dev libboost-regex-dev libboost-serialization-dev libboost-stacktrace-dev libboost-system-dev libboost-test-dev libboost-thread-dev libboost-timer-dev
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
RUN cd /tmp && git clone https://github.com/darafferty/LSMTool.git
RUN cd /tmp/LSMTool && git checkout 8730aaffd455cdd3071c68a1d56299ac1bece90c
RUN cd /tmp/LSMTool; pip3 install .

####################################################################
## LoSoTo
###################################################################
RUN cd /tmp && git clone https://github.com/revoltek/losoto.git
RUN cd /tmp/losoto; pip3 install .
#ENV PYTHONPATH /opt/losoto:$PYTHONPATH

#####################################################################
## CMake
#####################################################################
RUN cd /root && wget https://github.com/Kitware/CMake/releases/download/v3.18.3/cmake-3.18.3.tar.gz \
    && gunzip cmake-3.18.3.tar.gz && tar xvf cmake-3.18.3.tar
RUN cd /root/cmake-3.18.3 && cmake . -DCMAKE_USE_OPENSSL=OFF && make && make install
RUN export PATH=/usr/local/bin/:$PATH

#####################################################################
## idg
#####################################################################
RUN cd /opt &&  git clone https://gitlab.com/astron-idg/idg.git \
    && cd idg && mkdir build && cd build \
	&& cmake .. -DCMAKE_INSTALL_PREFIX=/opt/idg && make && make install

#####################################################################
## AOFlagger
#####################################################################
RUN cd /root && git clone https://gitlab.com/aroffringa/aoflagger.git
RUN cd /root/aoflagger && mkdir build && cd build && mkdir install \
    && cmake .. -DCASACORE_ROOT_DIR=/opt/casacore -DCMAKE_INSTALL_PREFIX=/root/aoflagger/build/install \
    && make && make install

#####################################################################
## DP3
#####################################################################
RUN cd /tmp/ && git clone https://github.com/lofar-astron/DP3.git \
	&& cd DP3 && git checkout d0fee2a \
	&& mkdir build

RUN cd /tmp/DP3/build && cmake -DCASACORE_ROOT_DIR=/opt/casacore -DAOFLAGGER_ROOT_DIR=/opt/aoflagger/ -DIDGAPI_LIBRARIES=/opt/idg/lib/libidg-api.so -DIDGAPI_INCLUDE_DIRS=/opt/idg/include -DCMAKE_INSTALL_PREFIX=/opt/DP3 .. 
RUN cd /tmp/DP3/build && make && make install

#####################################################################
# WSCLEAN
#####################################################################
RUN apt-get install -y libgsl-dev libhdf5-dev libfftw3-dev ipython3 python3-ipython
RUN mkdir /opt/wsclean
RUN cd /tmp && git clone https://gitlab.com/aroffringa/wsclean.git
RUN cd /tmp/wsclean && git checkout 6567e489
RUN mkdir /tmp/wsclean/build
RUN cd /tmp/wsclean/build && cmake -DCMAKE_PREFIX_PATH=/opt/casacore  -DBoost_FILESYSTEM_LIBRARY_DEBUG=/usr/lib/x86_64-linux-gnu/libboost_filesystem.so.1.65.1 -DBoost_FILESYSTEM_LIBRARY_RELEASE=/usr/lib/x86_64-linux-gnu/libboost_filesystem.so.1.65.1 -DBoost_DATE_TIME_LIBRARY_DEBUG=/usr/lib/x86_64-linux-gnu/libboost_date_time.so.1.65.1 -DBoost_DATE_TIME_LIBRARY_RELEASE=/usr/lib/x86_64-linux-gnu/libboost_date_time.so.1.65.1 -DBoost_PROGRAM_OPTIONS_LIBRARY_DEBUG=/usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.65.1 -DBoost_PROGRAM_OPTIONS_LIBRARY_RELEASE=/usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.65.1 -DBoost_SYSTEM_LIBRARY_DEBUG=/usr/lib/x86_64-linux-gnu/libboost_system.so.1.65.1 -DBoost_SYSTEM_LIBRARY_RELEASE=/usr/lib/x86_64-linux-gnu/libboost_system.so.1.65.1  -DCMAKE_INSTALL_PREFIX=/opt/wsclean/ ../ && make -j8 && make install
#-DLOFAR_STATION_RESPONSE_INCLUDE_DIR="/opt/lofar/include/StationResponse/" -DLOFAR_STATION_RESPONSE_LIB="/opt/lofar/lib/libstationresponse.so"
ENV PATH /opt/wsclean/bin:$PATH

###################################################################
## Dysco
##################################################################
#RUN cd /tmp; git clone https://github.com/aroffringa/dysco.git
#RUN cd /tmp/dysco && git checkout 7ef0f47
#RUN mkdir /tmp/dysco/build
#RUN cd /tmp/dysco/build && cmake -DCMAKE_INSTALL_PREFIX=/opt/dysco -DCASACORE_ROOT_DIR=/opt/casacore .. && make && make install
#
#ENV LD_LIBRARY_PATH /opt/dysco/lib:$LD_LIBRARY_PATH
#ENV LD_LIBRARY_PATH /opt/casacore/lib:$LD_LIBRARY_PATH

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
RUN pip3 install -U pip setuptools wheel
RUN python3 -m pip install tensorflow==1.8.0
RUN python3 -m pip install -U "/opt/DDFacet/[dft-support,moresane-support,testing-requirements,fits-beam-support]"
RUN cd /opt/DDFacet/ && python3 setup.py build && cd /

# Set MeqTrees Cattery path to installation directory
ENV MEQTREES_CATTERY_PATH /usr/local/lib/python3.6/dist-packages/Cattery/
ENV PYTHONPATH $MEQTREES_CATTERY_PATH:$PYTHONPATH

RUN python3 -m pip install numpy==1.17.0

RUN python3 -c "import Siamese"
RUN python3 -c "import bdsf"

# Perform some basic tests
#RUN DDF.py --help
#RUN ClusterCat.py --help
#RUN MakeMask.py --help
#RUN MakeCatalog.py --help
#RUN MakeModel.py --help
#RUN MaskDicoModel.py --help

#----------------------
# Entrypoint
#----------------------

# Rename lofar home folder as a "vanilla" home folder
#RUN mv /home/lofar /lofar_home_vanilla

# Give write access to anyone to the home folder so the entrypoint will be able
# to copy over the /home/matauser_vanilla into /home/lofar (for Singularity)
RUN chmod 777 /home

# Copy entrypoint
#COPY entrypoint.sh /

# Give right permissions
#RUN chmod 755 /entrypoint.sh

# Set entrypoint
#ENTRYPOINT ["/entrypoint.sh"]

RUN rm /usr/bin/python
RUN ln -s /usr/bin/python3.6 /usr/bin/python
#USER lofar
ENV OMP_NUM_THREADS 1
RUN ulimit -n 4000
ENV CONTAINER_NAME='LOFAR-LBA-pipeline'
