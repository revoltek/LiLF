FROM ubuntu:24.04
SHELL ["/bin/bash", "-c"]
ARG DEBIAN_FRONTEND=noninteractive
LABEL maintainer="Francesco de Gasperin"
LABEL version="20250228"

ENV MARCH "x86-64"
# AMD v3
#ENV MARCH "znver3"
# Intel
#ENV MARCH "broadwell"

USER root

# Support large mlocks
RUN echo "*        -   memlock     unlimited" > /etc/security/limits.conf

RUN apt-get update
RUN apt-get install -y cmake texinfo git wget flex bison sudo vim hdf5-tools rename python3-magic htop curl \
    python3-dev ipython3 python3-pip python3-setuptools python3-astropy python3-pyregion python3-regions python3-h5py python3-sshtunnel python3-pymysql python3-requests python3-numexpr python3-numpy python3-astroquery python3-aplpy python3-cytoolz python3-shapely python3-tqdm \
    libcfitsio-dev libpng-dev libxml2-dev libarmadillo-dev liblua5.3-dev libfftw3-dev python3-pybind11 wcslib-dev libgsl-dev libblas-dev libaio1t64 libaio-dev \
    libboost-dev libboost-date-time-dev libboost-filesystem-dev libboost-system-dev libboost-test-dev libboost-python-dev libboost-numpy-dev libboost-program-options-dev

# GCC 14. Should be removed after upgrading the container base to Ubuntu version with native GCC >= 14.
RUN apt-get install -y software-properties-common
RUN add-apt-repository -y universe && apt-get update
RUN apt-get install -y gcc-14 g++-14 gfortran-14
RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 14 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 14 && \
    update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-14 14

RUN curl -O http://launchpadlibrarian.net/646633572/libaio1_0.3.113-4_amd64.deb
RUN dpkg -i libaio1_0.3.113-4_amd64.deb && rm libaio1_0.3.113-4_amd64.deb

# Add group. We chose GID 65527 to try avoiding conflicts.
RUN groupadd -g 65527 lofar

# Add user. We chose UID 65527 to try avoiding conflicts.
RUN useradd lofar -d /home/lofar -u 65527 -g 65527 -m -s /bin/bash

# Setup home directories
RUN chown lofar:lofar /home/lofar

# To avoid installing Python virtualenv
ENV PIP_BREAK_SYSTEM_PACKAGES 1

##################################################################
## AMD specific

#RUN if [ $MARCH == "znver3" ]; then \
#    cd /opt && wget https://download.amd.com/developer/eula/aocl/aocl-4-2/aocl-linux-aocc-4.2.0_1_amd64.deb && dpkg -i aocl-linux-aocc-4.2.0_1_amd64.deb \
#    && cd /opt/AMD/aocl/aocl-linux-aocc-4.2.0/aocc && ./set_aocl_interface_symlink.sh lp64; fi
#ENV CPLUS_INCLUDE_PATH /opt/AMD/aocl/aocl-linux-aocc-4.2.0/aocc/include:$CPLUS_INCLUDE_PATH
#ENV CPATH /opt/AMD/aocl/aocl-linux-aocc-4.2.0/aocc/include:$CPATH
#ENV LD_LIBRARY_PATH /opt/AMD/aocl/aocl-linux-aocc-4.2.0/aocc/lib:$LD_LIBRARY_PATH
#ENV AOCL_ROOT /opt/AMD/aocl/aocl-linux-aocc-4.2.0/aocc

#RUN mkdir /opt/aocl && cd /opt/aocl && wget https://download.amd.com/developer/eula/aocl/aocl-4-2/aocl-linux-aocc-4.2.0.tar.gz && tar xf aocl-linux-aocc-4.2.0.tar.gz
#RUN cd /opt/aocl/aocl-linux-aocc-4.2.0 && ./install.sh -t /opt/aocl -i lp64
#RUN source /opt/aocl/4.2.0/amd-libs.cfg
#ENV CPLUS_INCLUDE_PATH /opt/aocl/4.2.0/include:$CPLUS_INCLUDE_PATH
#ENV CPATH /opt/aocl/4.2.0/include:$CPATH
#ENV LD_LIBRARY_PATH /opt/aocl/4.2.0/lib:$LD_LIBRARY_PATH
#ENV AOCL_ROOT /opt/aocl/4.2.0

# Get libflang.so
#RUN mkdir /opt/flang && cd /opt/flang && wget https://github.com/flang-compiler/flang/releases/download/flang_20190329/flang-20190329-x86-70.tgz && tar xf flang-20190329-x86-70.tgz
#ENV LD_LIBRARY_PATH /opt/flang/lib:$LD_LIBRARY_PATH
#ENV PATH /opt/flang/bin:$PATH
#ENV CPATH /opt/flang/include:$CPATH
#ENV CPLUS_INCLUDE_PATH /opt/flang/include:$CPLUS_INCLUDE_PATH

# Explicitely reinstall FFTW since it wasn't built with threads.
#RUN rm -rf $AOCL_ROOT/amd-fftw $AOCL_ROOT/lib/libfftw*.so $AOCL_ROOT/include/fftw3*
#RUN mkdir /opt/fftw/ && cd /opt/fftw/ && git clone https://github.com/amd/amd-fftw.git
# Install single precision version
#RUN cd /opt/fftw/amd-fftw && ./configure --enable-float --enable-threads --enable-sse2 --enable-avx --enable-avx2 --enable-avx512 --enable-openmp --enable-shared --enable-amd-opt --enable-amd-mpifft --prefix=/opt/aocl/4.0 && make -j `nproc --all` && make install && make clean
# Install double precision version
#RUN cd /opt/fftw/amd-fftw && ./configure --enable-threads --enable-sse2 --enable-avx --enable-avx2 --enable-avx512 --enable-openmp --enable-shared --enable-amd-opt --enable-amd-mpifft --prefix=/opt/aocl/4.0 && make -j `nproc --all` && make install

## end AMD specific

##################################################################
## Intel specific
RUN if [ $MARCH == "broadwell" ]; then \
    apt-get install -y gpg-agent && \
    wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null && \
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list && \
    apt-get update && \
    apt-get install -y intel-oneapi-mkl intel-oneapi-mkl-devel; fi

# Outside if condition, but does no harm if empty.
ENV PKG_CONFIG_PATH /opt/intel/oneapi/mkl/latest/lib/pkgconfig:$PKG_CONFIG_PATH
ENV CMAKE_PREFIX_PATH /opt/intel/oneapi/mkl/latest/lib/cmake:$CMAKE_PREFIX_PATH
ENV LIBRARY_PATH /opt/intel/oneapi/mkl/latest/lib:$LIBRARY_PATH
ENV LD_LIBRARY_PATH /opt/intel/oneapi/mkl/latest/lib:$LD_LIBRARY_PATH
ENV MKLROOT /opt/intel/oneapi/mkl/latest
ENV PATH /opt/intel/oneapi/mkl/latest/bin:$PATH
ENV CPATH /opt/intel/oneapi/mkl/latest/include:$CPATH
## end Intel specific

RUN apt-get clean

##################################################################
## CASACORE v3.6.1
##################################################################
RUN cd /opt && git clone --single-branch --branch master https://github.com/casacore/casacore.git \
    && cd /opt/casacore && git checkout tags/v3.6.1
RUN cd /opt/casacore && mkdir data && cd data && wget --retry-connrefused ftp://anonymous@ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar && tar xf WSRT_Measures.ztar && rm WSRT_Measures.ztar
RUN cd /opt/casacore && mkdir build && cd build && cmake -DPORTABLE=False -DCMAKE_BUILD_TYPE=Release -DDATA_DIR=/opt/casacore/data -DBUILD_PYTHON=False -DBUILD_PYTHON3=True -DUSE_OPENMP=True -DUSE_HDF5=True .. && make -j `nproc --all` && make install
RUN cd /opt/casacore && rm -r $(ls -A | grep -v data)

#####################################################################
## CASACORE-python v3.6.1
#####################################################################
RUN cd /opt && git clone https://github.com/casacore/python-casacore.git \
    && cd /opt/python-casacore && git checkout tags/v3.6.1
RUN cd /opt/python-casacore && pip install .

#####################################################################
## EveryBeam v0.7.0
#####################################################################
RUN cd /opt && git clone https://git.astron.nl/RD/EveryBeam.git \
    && cd /opt/EveryBeam && git checkout tags/v0.7.0
RUN cd /opt/EveryBeam && mkdir build && cd build \
    && cmake -DDOWNLOAD_LOBES=Off -DBUILD_WITH_PYTHON=On .. \
    && make -j `nproc --all` && make install
RUN rm -rf /opt/EveryBeam/
ENV EVERYBEAM_DATADIR /usr/local/share/everybeam
    
#####################################################################
## idg v1.2.0
#####################################################################
RUN cd /opt && git clone https://gitlab.com/astron-idg/idg.git \
    && cd /opt/idg && git checkout tags/1.2.0
RUN if [ $MARCH == "x86-64" ]; then cd /opt/idg && mkdir build && cd build && cmake -DPORTABLE=True -DCMAKE_INSTALL_PREFIX:PATH=/ .. && make -j `nproc --all` && make install; fi
RUN if [ $MARCH != "x86-64" ]; then cd /opt/idg && mkdir build && cd build && cmake -DPORTABLE=False -DCMAKE_INSTALL_PREFIX:PATH=/ .. && make -j `nproc --all` && make install; fi

#####################################################################
## AOFlagger v3.4.0
#####################################################################
RUN cd /opt && git clone https://gitlab.com/aroffringa/aoflagger.git \
    && cd /opt/aoflagger && git checkout tags/v3.4.0
RUN if [ $MARCH == "x86-64" ]; then cd /opt/aoflagger && mkdir build && cd build && mkdir install && cmake -DPORTABLE=True .. && make -j `nproc --all` && make install; fi
RUN if [ $MARCH != "x86-64" ]; then cd /opt/aoflagger && mkdir build && cd build && mkdir install && cmake -DPORTABLE=False .. && make -j `nproc --all` && make install; fi

#####################################################################
## Dysco v1.3
## not needed, it comes with casacore
#####################################################################
#RUN cd /opt && git clone https://github.com/aroffringa/dysco.git \
#    && cd /opt/dysco && git checkout tags/v1.3
#RUN cd /opt/dysco && mkdir build && cd build \
#    && cmake -DPORTABLE=True .. \
#    && make -j `nproc --all` && make install

#####################################################################
## LofarStMan
## not needed, only for raw (from the correlator) MS
#####################################################################
#RUN cd /opt && git clone https://github.com/lofar-astron/LofarStMan
#RUN cd /opt/LofarStMan && mkdir build && cd build \
#    && cmake -DPORTABLE=True .. \
#    && make -j `nproc --all` && make install

####################################################################
## Dp3 master (03/02/25)
#####################################################################
RUN ls /usr/local/share/everybeam
RUN ls /usr/local/lib
RUN cd /opt && git clone https://git.astron.nl/RD/DP3.git \
    && cd /opt/DP3 && git checkout f475bc30

RUN if [ $MARCH == "x86-64" ]; then \
    cd /opt/DP3 && mkdir build && cd build \
    && cmake .. \
    && make -j `nproc --all` && make install; fi

RUN if [ $MARCH == "broadwell" ]; then \
    cd /opt/DP3 && mkdir build && cd build \
    && cmake -DTARGET_CPU=${MARCH} .. \
    && make -j `nproc --all` && make install; fi

RUN if [ $MARCH == "znver3" ]; then \
    cd /opt/DP3 && mkdir build && cd build \
    && cmake -DTARGET_CPU=${MARCH} .. \
    && make -j `nproc --all` && make install; fi

#        -DBLAS_flexiblas_LIBRARY=${AOCL_ROOT}/lib/libblis.so \
#        -DFFTW3F_LIB="${AOCL_ROOT}/lib/libfftw3f.so" \
#        -DFFTW3F_THREADS_LIBRARY="${AOCL_ROOT}/lib/libfftw3f_threads.so" \
#        -DFFTW3_INCLUDE_DIR="${AOCL_ROOT}/include/" \
#        -DLAPACK_flexiblas_LIBRARY=${AOCL_ROOT}/lib/libflame.so \

RUN rm -rf /opt/DP3/
# add "import dp3"
ENV PYTHONPATH /usr/local/lib/python3.12/site-packages/:$PYTHONPATH

#####################################################################
## Wsclean master (11/2/25)
#####################################################################
RUN cd /opt && git clone https://gitlab.com/aroffringa/wsclean.git \
    && cd /opt/wsclean && git checkout af847af6b39986b18dd06c8e0846fb01f87f36c3

RUN if [ $MARCH == "x86-64" ]; then \
    cd /opt/wsclean && mkdir build && cd build \
    && cmake .. \
    && make -j `nproc --all` && make install; fi

RUN if [ $MARCH == "broadwell" ]; then \
    cd /opt/wsclean && mkdir build && cd build \
    && cmake -DTARGET_CPU=${MARCH} .. \
    && make -j `nproc --all` && make install; fi

RUN if [ $MARCH == "znver3" ]; then \
    cd /opt/wsclean && mkdir build && cd build \
    && cmake -DTARGET_CPU=${MARCH} .. \
    && make -j `nproc --all` && make install; fi

#    -DFFTW3_LIB="${AOCL_ROOT}/lib/libfftw3.so" \
#    -DFFTW3_THREADS_LIB="${AOCL_ROOT}/lib/libfftw3_threads.so" \
#    -DFFTW3F_LIB="${AOCL_ROOT}/lib/libfftw3f.so" \
#    -DFFTW3F_THREADS_LIB="${AOCL_ROOT}/lib/libfftw3f_threads.so" \
#    -DFFTW3_INCLUDE_DIR="${AOCL_ROOT}/include" \
#    -DBLAS_openblas_LIBRARY="${AOCL_ROOT}/lib/libblis.so" \
#    -DLAPACK_LIBRARIES="${AOCL_ROOT}/lib/libflame.so" \
#    -DCMAKE_CXX_FLAGS="-march=${MARCH} -L${AOCL_ROOT}/lib -lamdlibm -lm -O3 -Ofast" \
#    -DCMAKE_C_FLAGS="-march=${MARCH} -L${AOCL_ROOT}/lib -lamdlibm -lm -O3 -Ofast" \

#RUN rm -rf /opt/wsclean/
  
#####################################################################
## LTA
#####################################################################
RUN apt-get install -y python3-cx-oracle
RUN cd /opt && wget https://lta.lofar.eu/software/lofar_lta-2.8.0.tar.gz && tar xf lofar_lta-2.8.0.tar.gz && rm lofar_lta-2.8.0.tar.gz
# local patched version with correct url for download
#COPY lofar_lta-2.8.0 /opt/lofar_lta-2.8.0
RUN cd /opt/lofar_lta-2.8.0 && python3 setup.py install_oracle
RUN cd /opt/lofar_lta-2.8.0 && python3 setup.py install
RUN ln -s /usr/lib/x86_64-linux-gnu/libaio.so.1t64 /usr/lib/instantclient_11_2/libaio.so.1
ENV LD_LIBRARY_PATH /usr/lib/instantclient_11_2:$LD_LIBRARY_PATH

#####################################################################
## Pybdsf (v1.11.0)
#####################################################################
RUN cd /opt && git clone https://github.com/lofar-astron/PyBDSF.git \
    && cd /opt/PyBDSF; git checkout tags/v1.11.0
RUN cd /opt/PyBDSF && pip install .

#####################################################################
## LSMtool (fix_query_tempdir)
#####################################################################
RUN cd /opt && git clone https://git.astron.nl/RD/LSMTool.git \
    && cd /opt/LSMTool; git checkout fix_query_tempdir
RUN cd /opt/LSMTool && pip install . --upgrade

####################################################################
## LoSoTo (master 28/2/25)
###################################################################
RUN cd /opt && git clone https://github.com/revoltek/losoto.git \
    && cd /opt/losoto && git checkout da8d88a3d42a72a0e91e6f36ced072c73a34e6ee
RUN cd /opt/losoto && pip install .

#####################################################################
# msoverview
#####################################################################
#RUN cd /opt && mkdir -p msoverview/src && cd msoverview/src && wget https://git.astron.nl/ro/lofar/-/raw/LOFAR-Release-4_0_17/CEP/MS/src/msoverview.cc \
#    && gcc -I/usr/local/include/casacore -L/usr/local/lib msoverview.cc -o /usr/local/bin/msoverview -lcasa_casa -lcasa_ms -lcasa_tables -lstdc++

#####################################################################
## BUILD DDF FROM SOURCE
#####################################################################
#RUN cd /tmp && git clone https://github.com/cyriltasse/DDFacet.git
#RUN cd /tmp/DDFacet && git checkout c0bd5f9ccaa9b7c28376e3e2674428902fa5bf85

#COPY DDFacet /tmp/DDFacet
#RUN mkdir /opt/DDFacet
#
## Copy DDFacet and SkyModel into the image
#RUN cp -a /tmp/DDFacet/DDFacet /opt/DDFacet/DDFacet && \
#    cp -a /tmp/DDFacet/SkyModel /opt/DDFacet/SkyModel && \
#    cp -a /tmp/DDFacet/MANIFEST.in /opt/DDFacet/MANIFEST.in && \
#    cp -a /tmp/DDFacet/setup.py /opt/DDFacet/setup.py && \
#    cp -a /tmp/DDFacet/setup.cfg /opt/DDFacet/setup.cfg && \
#    cp -a /tmp/DDFacet/README.rst /opt/DDFacet/README.rst && \
#    cp -a /tmp/DDFacet/pyproject.toml /opt/DDFacet/pyproject.toml && \
#    cp -a /tmp/DDFacet/.git /opt/DDFacet/.git && \
#    cp -a /tmp/DDFacet/.gitignore /opt/DDFacet/.gitignore && \
#    cp -a /tmp/DDFacet/.gitmodules /opt/DDFacet/.gitmodules
#
## Finally install DDFacet
#RUN rm -rf /opt/DDFacet/DDFacet/cbuild
##RUN pip3 install -U pip setuptools wheel
##RUN python3 -m pip install tensorflow-cpu
#
## Force numba <0.55 otherwise it breaks for numpy
##RUN pip3 install numba==0.54.1
## Forse ipython not to be upgraded too much otherwise we lose compatibility with come astropy
##RUN pip3 install ipython==7.17.0
#
#RUN python3 -m pip install -U "/opt/DDFacet/[moresane-support,testing-requirements,fits-beam-support]"
#RUN cd /opt/DDFacet/ && python3 setup.py build
#RUN cd /opt/DDFacet/ && python3 setup.py install
#
## Set MeqTrees Cattery path to installation directory
##ENV MEQTREES_CATTERY_PATH /usr/local/lib/python3.6/dist-packages/Cattery/
##ENV PYTHONPATH $MEQTREES_CATTERY_PATH:$PYTHONPATH

#####################################################################
## LiLF
#####################################################################
#RUN cd /opt && git clone https://github.com/revoltek/LiLF.git
#RUN cd /opt/LiLF && git checkout LBAdevel
#ENV PYTHONPATH /opt/LiLF:$PYTHONPATH
#ENV PATH /opt/LiLF/scripts:$PATH
   
#####################################################################

RUN mkdir /usr/share/casacore/
RUN ln -s /usr/share/casacore/data/ /opt/casacore/data

ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN ln -s /usr/bin/ipython3 /usr/bin/ipython
ENV HDF5_USE_FILE_LOCKING FALSE
ENV OMP_NUM_THREADS 1
ENV OPENBLAS_NUM_THREADS 1
#ENTRYPOINT [ "ulimit -n 8192" ]
