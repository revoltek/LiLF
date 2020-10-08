FROM lofar/lofar-pipeline:LOFAR-Release-4_0_16

USER root

# Support large mlocks
RUN echo "*        -   memlock     unlimited" > /etc/security/limits.conf

RUN apt-get update -y
RUN apt-get install -y libcfitsio-dev libpng-dev git cmake python3-dev libxml2-dev python3-pip

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
## LSMtool
#####################################################################

RUN cd /tmp && git clone https://github.com/darafferty/LSMTool.git
RUN cd /tmp/LSMTool && git checkout 8730aaffd455cdd3071c68a1d56299ac1bece90c

RUN cd /tmp/LSMTool; pip3 install .

#####################################################################
## LiLF
#####################################################################

RUN cd /opt && git clone https://github.com/revoltek/LiLF.git
RUN cd /opt/LiLF && git checkout d4d74188b819e7932bfcdb6a27148d6a5c17ae30



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
RUN python3 -m pip install pybind11
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
RUN DDF.py --help
RUN ClusterCat.py --help
RUN MakeMask.py --help
RUN MakeCatalog.py --help
RUN MakeModel.py --help
RUN MaskDicoModel.py --help



RUN pip3 install pyregion 
RUN pip3 install astropy --upgrade


#####################################################################
# WSCLEAN
#####################################################################
RUN apt-get install -y libgsl-dev libhdf5-dev libfftw3-dev ipython3 python3-ipython
RUN mkdir /opt/wsclean
RUN cd /tmp; git clone https://gitlab.com/aroffringa/wsclean.git
RUN cd /tmp/wsclean && git checkout 6567e489
RUN mkdir /tmp/wsclean/build
RUN cd /tmp/wsclean/build && cmake -DCMAKE_PREFIX_PATH=/opt/casacore  -DBoost_FILESYSTEM_LIBRARY_DEBUG=/usr/lib/x86_64-linux-gnu/libboost_filesystem.so.1.65.1 -DBoost_FILESYSTEM_LIBRARY_RELEASE=/usr/lib/x86_64-linux-gnu/libboost_filesystem.so.1.65.1 -DBoost_DATE_TIME_LIBRARY_DEBUG=/usr/lib/x86_64-linux-gnu/libboost_date_time.so.1.65.1 -DBoost_DATE_TIME_LIBRARY_RELEASE=/usr/lib/x86_64-linux-gnu/libboost_date_time.so.1.65.1 -DBoost_PROGRAM_OPTIONS_LIBRARY_DEBUG=/usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.65.1 -DBoost_PROGRAM_OPTIONS_LIBRARY_RELEASE=/usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.65.1 -DBoost_SYSTEM_LIBRARY_DEBUG=/usr/lib/x86_64-linux-gnu/libboost_system.so.1.65.1 -DBoost_SYSTEM_LIBRARY_RELEASE=/usr/lib/x86_64-linux-gnu/libboost_system.so.1.65.1  -DCMAKE_INSTALL_PREFIX=/opt/wsclean/ ../ && make -j8 && make install
#-DLOFAR_STATION_RESPONSE_INCLUDE_DIR="/opt/lofar/include/StationResponse/" -DLOFAR_STATION_RESPONSE_LIB="/opt/lofar/lib/libstationresponse.so"
ENV PATH /opt/wsclean/bin:$PATH

#----------------------
# Entrypoint
#----------------------

# Copy entrypoint
COPY entrypoint.sh /

# Give right permissions
RUN chmod 755 /entrypoint.sh

# Set entrypoint
ENTRYPOINT ["/entrypoint.sh"]

RUN rm /usr/bin/python; ln -s /usr/bin/python3.6 /usr/bin/python
USER lofar
ENV CONTAINER_NAME='ddfacet-gain_screens_premerge_rebased-lofar-release-4_0_16'
