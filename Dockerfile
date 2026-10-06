FROM firedrakeproject/firedrake-vanilla:2025-01@sha256:83395337ecbc40ea5d00211c54d77b5114328b2def700581d98d07cf42864cc1
USER root
ENV PATH="/home/firedrake/firedrake/bin:${PATH}"
ENV PYOP2_CACHE_DIR=/tmp/pyop2-cache
ENV FIREDRAKE_TSFC_KERNEL_CACHE_DIR=/tmp/tsfc-cache
RUN python -m pip install \
    git+https://github.com/icepack/icepack.git@c9a29780cd0f7d068d206cb5a170fa367a7655b0 \
    nbclient nbconvert pytest siphash24
WORKDIR /work
