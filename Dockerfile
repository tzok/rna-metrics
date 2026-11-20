# syntax=docker/dockerfile:1

# ---- Build stage ---------------------------------------------------------
# Everything that is only needed to compile things lives here, so the
# runtime image stays free of toolchains.
FROM python:3.15-slim AS builder

# g++/wget build USalign; cmake, make, flex and bison build the mmcif
# extension, which has no wheels for Python 3.15 yet and is compiled from
# source.
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ \
    make \
    wget \
    cmake \
    flex \
    bison \
    && rm -rf /var/lib/apt/lists/*

# Download and compile USalign. It is linked statically so the resulting
# binary needs nothing from this stage at runtime.
WORKDIR /tmp
RUN wget https://zhanggroup.org/US-align/bin/module/USalign.cpp && \
    g++ -static -O3 -ffast-math -o USalign USalign.cpp && \
    mv USalign /usr/local/bin/ && \
    rm USalign.cpp

# Install Python dependencies into /install. They are copied into the
# runtime image's /usr/local, so the interpreter picks them up without any
# PATH tweaks and runtime `pip install` keeps working as before.
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir --root=/install -r /tmp/requirements.txt


# ---- Runtime stage -------------------------------------------------------
FROM python:3.15-slim

# USalign binary and the Python packages built above
COPY --from=builder /usr/local/bin/USalign /usr/local/bin/USalign
COPY --from=builder /install/usr/local /usr/local

# Set the working directory in the container
WORKDIR /app

# Copy all Python scripts and other files
COPY requirements.txt /app/
COPY clashscore.py /app/
COPY inf.py /app/
COPY lddt.py /app/
COPY mcq.py /app/
COPY rmsd.py /app/
COPY tm_score.py /app/
COPY torsion.py /app/

# Make all scripts executable
RUN chmod +x /app/*.py

# Set environment variable for scripts location
ENV PATH="/app:${PATH}"

# Default command (can be overridden)
CMD ["python", "-c", "import sys; print('Available scripts: clashscore.py, inf.py, lddt.py, mcq.py, rmsd.py, tm_score.py, torsion.py')"]

COPY pytest.ini /app
COPY test_requirements.txt /app
COPY tests /app/tests
