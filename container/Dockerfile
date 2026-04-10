FROM ubuntu:25.10
SHELL ["/bin/bash", "-c"]
ARG DEBIAN_FRONTEND=noninteractive
LABEL maintainer="Francesco de Gasperin"
LABEL version="20260215"

ENV MARCH="x86-64-v2"
# AMD v3
#ENV MARCH="znver3"
# Intel
#ENV MARCH="broadwell"

USER root

# Support large mlocks
RUN echo "*        -   memlock     unlimited" > /etc/security/limits.conf

RUN apt-get update && \
    apt-get install -y --no-install-recommends bison build-essential cmake curl flex gfortran git hdf5-tools htop ipython3 libaio1t64 libaio-dev libarmadillo-dev \
    libblas-dev \
    libboost-date-time-dev libboost-dev libboost-filesystem-dev libboost-numpy-dev libboost-program-options-dev libboost-python-dev libboost-system-dev \
    libcfitsio-dev libdeflate-dev libfftw3-dev libgsl-dev liblua5.3-dev libpng-dev libxml2-dev \
    libopenmpi-dev openmpi-bin \
    python3-cytoolz python3-dev python3-magic python3-numpy python3-pip python3-pybind11 python3-pymysql python3-requests python3-setuptools python3-sshtunnel \
    rename rust-all sudo texinfo vim wcslib-dev wget

RUN curl -O https://launchpadlibrarian.net/646633572/libaio1_0.3.113-4_amd64.deb
RUN dpkg -i libaio1_0.3.113-4_amd64.deb && rm libaio1_0.3.113-4_amd64.deb

# Add group. We chose GID 65527 to try avoiding conflicts.
RUN groupadd -g 65527 lofar

# Add user. We chose UID 65527 to try avoiding conflicts.
RUN useradd lofar -d /home/lofar -u 65527 -g 65527 -m -s /bin/bash

# Setup home directories
RUN chown lofar:lofar /home/lofar

# To avoid installing Python virtualenv
ENV PIP_BREAK_SYSTEM_PACKAGES=1

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
    apt-get install -y gpg-agent gnupg && \
    wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null && \
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list && \
    apt-get update && \
    apt-get install -y intel-oneapi-mkl intel-oneapi-mkl-devel \
    && apt-get clean \
    && mv /opt/intel/oneapi/compiler/latest/lib/libiomp5.so /opt/intel/oneapi/mkl/latest/lib/ \
    && rm -rf /opt/intel/oneapi/compiler/ \
    && find /opt/intel/oneapi -type f \( -name "*.a" -o -name "libmkl_sycl*" \) -delete; fi

# Silence UndefinedVar warning from BuildKit Dockerfile linter
ENV PKG_CONFIG_PATH="" CMAKE_PREFIX_PATH="" LIBRARY_PATH="" LD_LIBRARY_PATH="" CPATH="" PYTHONPATH=""
# Outside MKLs if condition, but does no harm if empty
ENV PKG_CONFIG_PATH=/opt/intel/oneapi/mkl/latest/lib/pkgconfig:$PKG_CONFIG_PATH
ENV CMAKE_PREFIX_PATH=/opt/intel/oneapi/mkl/latest/lib/cmake:$CMAKE_PREFIX_PATH
ENV LIBRARY_PATH=/opt/intel/oneapi/mkl/latest/lib:$LIBRARY_PATH
ENV LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/latest/lib:$LD_LIBRARY_PATH
ENV MKLROOT=/opt/intel/oneapi/mkl/latest
ENV PATH=/opt/intel/oneapi/mkl/latest/bin:$PATH
ENV CPATH=/opt/intel/oneapi/mkl/latest/include:$CPATH
## end Intel specific

##################################################################
## CASACore master (2/9/25)
##################################################################
RUN cd /opt && \
    git clone --single-branch --filter=blob:none --branch master https://github.com/casacore/casacore.git && \
    cd casacore && \
    git checkout 450a568 && \
    mkdir data && cd data && \
    wget --retry-connrefused --no-verbose https://www.astron.nl/iers/WSRT_Measures.ztar && \
    tar xf WSRT_Measures.ztar && \
    rm WSRT_Measures.ztar && \
    cd .. && \
    mkdir build && cd build && \
    if [ "$MARCH" = "broadwell" ]; then \
          CMAKE_EXTRA="-DBLA_VENDOR=Intel10_64lp"; else CMAKE_EXTRA=""; fi && \
    cmake -Wno-dev -DPORTABLE=False -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release \
          -DDATA_DIR=/opt/casacore/data -DUSE_OPENMP=True -DUSE_HDF5=True \
          -DBUILD_SISCO=ON $CMAKE_EXTRA .. && \
    make -j $(nproc) && \
    make install && \
    cd /opt/casacore && \
    find . -mindepth 1 -maxdepth 1 ! -name 'data' -exec rm -rf {} +

#####################################################################
## CASACore-python v3.7.1
#####################################################################
#RUN cd /opt && git clone https://github.com/casacore/python-casacore.git \
#    && cd /opt/python-casacore && git checkout tags/v3.7.1 && rm -rf .git/objects/pack/
#RUN cd /opt/python-casacore && pip install .

#####################################################################
## EveryBeam v0.8.2
#####################################################################
RUN cd /opt && git clone https://git.astron.nl/RD/EveryBeam.git \
    && cd /opt/EveryBeam && git checkout tags/v0.8.2 && rm -rf .git/objects/pack/
RUN cd /opt/EveryBeam && mkdir build && cd build \
    && cmake -Wno-dev -DDOWNLOAD_LOBES=Off -DBUILD_WITH_PYTHON=On .. \
    && make -j `nproc --all` && make install \
    && rm -rf /opt/EveryBeam/

ENV EVERYBEAM_DATADIR=/usr/local/share/everybeam

#####################################################################
## idg v1.2.0
#####################################################################
RUN cd /opt && git clone https://gitlab.com/astron-idg/idg.git \
    && cd /opt/idg && git checkout tags/1.2.0 && rm -rf .git/objects/pack/
RUN if [ $MARCH == "x86-64-v2" ]; then cd /opt/idg && mkdir build && cd build && cmake -DPORTABLE=True -DCMAKE_INSTALL_PREFIX:PATH=/ .. && make -j `nproc --all` && make install; fi
RUN if [ $MARCH != "x86-64-v2" ]; then cd /opt/idg && mkdir build && cd build && cmake -DPORTABLE=False -DCMAKE_INSTALL_PREFIX:PATH=/ .. && make -j `nproc --all` && make install; fi

#####################################################################
## AOFlagger (master 10/2/26)
#####################################################################
RUN cd /opt && git clone https://gitlab.com/aroffringa/aoflagger.git \
    && cd /opt/aoflagger && git checkout 398fcd5b && rm -rf .git/objects/pack/
RUN if [ $MARCH == "x86-64-v2" ]; then cd /opt/aoflagger && mkdir build && cd build && mkdir install && cmake -Wno-dev -DPORTABLE=True -DBUILD_TESTING=OFF -DENABLE_GUI=OFF .. && make -j `nproc --all` && make install; fi
RUN if [ $MARCH != "x86-64-v2" ]; then cd /opt/aoflagger && mkdir build && cd build && mkdir install && cmake -Wno-dev -DPORTABLE=False -DBUILD_TESTING=OFF -DENABLE_GUI=OFF .. && make -j `nproc --all` && make install; fi

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
## DP3 (master 12/2/26)
#####################################################################
RUN ls /usr/local/share/everybeam
RUN ls /usr/local/lib
RUN cd /opt && git clone https://git.astron.nl/RD/DP3.git \
    && cd /opt/DP3 && git checkout 56995db9 && rm -rf .git/objects/pack/

RUN if [ $MARCH == "x86-64-v2" ]; then \
    cd /opt/DP3 && mkdir build && cd build \
    && cmake -Wno-dev .. \
    && make -j `nproc --all` && make install \
    && rm -rf /opt/DP3/; fi

RUN if [ $MARCH == "broadwell" ]; then \
    cd /opt/DP3 && mkdir build && cd build \
    && cmake -Wno-dev -DTARGET_CPU=${MARCH} .. \
    && make -j `nproc --all` && make install \
    && rm -rf /opt/DP3/; fi

RUN if [ $MARCH == "znver3" ]; then \
    cd /opt/DP3 && mkdir build && cd build \
    && cmake -Wno-dev -DTARGET_CPU=${MARCH} .. \
    && make -j `nproc --all` && make install \
    && rm -rf /opt/DP3/; fi

#        -DBLAS_flexiblas_LIBRARY=${AOCL_ROOT}/lib/libblis.so \
#        -DFFTW3F_LIB="${AOCL_ROOT}/lib/libfftw3f.so" \
#        -DFFTW3F_THREADS_LIBRARY="${AOCL_ROOT}/lib/libfftw3f_threads.so" \
#        -DFFTW3_INCLUDE_DIR="${AOCL_ROOT}/include/" \
#        -DLAPACK_flexiblas_LIBRARY=${AOCL_ROOT}/lib/libflame.so \

# add "import dp3"
ENV PYTHONPATH=/usr/local/lib/python3.13/site-packages/:$PYTHONPATH

#####################################################################
## Wsclean (master 14/2/26)
#####################################################################
RUN cd /opt && git clone https://gitlab.com/aroffringa/wsclean.git \
    && cd /opt/wsclean && git checkout e4c63ab8 && rm -rf .git/objects/pack/

RUN if [ $MARCH == "x86-64-v2" ]; then \
    cd /opt/wsclean && mkdir build && cd build \
    && cmake -Wno-dev .. \
    && make -j `nproc --all` && make install \
    && rm -rf /opt/wsclean; fi

RUN if [ $MARCH == "broadwell" ]; then \
    cd /opt/wsclean && mkdir build && cd build \
    && cmake -Wno-dev -DTARGET_CPU=${MARCH} .. \
    && make -j `nproc --all` && make install \
    && rm -rf /opt/wsclean; fi

RUN if [ $MARCH == "znver3" ]; then \
    cd /opt/wsclean && mkdir build && cd build \
    && cmake -Wno-dev -DTARGET_CPU=${MARCH} .. \
    && make -j `nproc --all` && make install \
    && rm -rf /opt/wsclean; fi

#    -DFFTW3_LIB="${AOCL_ROOT}/lib/libfftw3.so" \
#    -DFFTW3_THREADS_LIB="${AOCL_ROOT}/lib/libfftw3_threads.so" \
#    -DFFTW3F_LIB="${AOCL_ROOT}/lib/libfftw3f.so" \
#    -DFFTW3F_THREADS_LIB="${AOCL_ROOT}/lib/libfftw3f_threads.so" \
#    -DFFTW3_INCLUDE_DIR="${AOCL_ROOT}/include" \
#    -DBLAS_openblas_LIBRARY="${AOCL_ROOT}/lib/libblis.so" \
#    -DLAPACK_LIBRARIES="${AOCL_ROOT}/lib/libflame.so" \
#    -DCMAKE_CXX_FLAGS="-march=${MARCH} -L${AOCL_ROOT}/lib -lamdlibm -lm -O3 -Ofast" \
#    -DCMAKE_C_FLAGS="-march=${MARCH} -L${AOCL_ROOT}/lib -lamdlibm -lm -O3 -Ofast" \

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
ENV LD_LIBRARY_PATH=/usr/lib/instantclient_11_2:$LD_LIBRARY_PATH

#####################################################################
## PyBDSF master 3/11/25 | LSMtool master 27/3/25 | LoSoTo master 12/2/26 | RMExtract 0.5.1 | spinifex master 14/3/25 | mocpy master 21/11/25 | spherical-geometry v1.3.3
#####################################################################
RUN pip install --root-user-action=ignore \
    git+https://github.com/lofar-astron/PyBDSF.git@8b33037 \
    git+https://git.astron.nl/RD/LSMTool.git@367dfd3e \
    git+https://github.com/revoltek/losoto.git@86358eb \
    git+https://github.com/lofar-astron/RMextract@v0.5.1 \
    git+https://git.astron.nl/RD/spinifex@1616fe1e \
    git+https://github.com/cds-astro/mocpy@945b877 \
    git+https://github.com/spacetelescope/spherical_geometry@1.3.3 \
    && pip cache purge \
    && rm -rf /root/.cache/pip \
    && rm -rf /root/.cargo/{git,registry}

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

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN ln -s /usr/bin/ipython3 /usr/bin/ipython
ENV HDF5_USE_FILE_LOCKING=FALSE
ENV OMP_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1
#ENTRYPOINT [ "ulimit -n 8192" ]
