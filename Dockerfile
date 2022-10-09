FROM python:3.9-slim
RUN apt-get update
RUN apt-get -y install git gfortran
RUN pip install --no-cache --upgrade pip
RUN pip install --no-cache notebook jupyterlab
RUN git clone https://github.com/mgaitan/fortran_magic && \
    cd fortran_magic && \
    python setup.py install
ENV HOME=/tmp
WORKDIR ${HOME}
