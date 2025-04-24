# Use official Python slim image based on Debian
FROM python:3.11

# Set DEBIAN_FRONTEND to noninteractive to prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini

# Install tzdata first and set timezone
RUN apt-get update && \
    apt-get install -y --no-install-recommends tzdata && \
    ln -fs /usr/share/zoneinfo/Etc/UTC /etc/localtime && \
    dpkg-reconfigure --frontend noninteractive tzdata && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install build tools and required C libraries for Python packages
# python3-dev, python3-pip, python3-venv are not needed as they come with the base image
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    libhdf5-dev \
    libsodium-dev \
    libzmq3-dev \
    gcc \
    g++ \
    # Install sudo for the adduser command below if needed
    sudo \
    # Keep essential utilities
    curl \
    wget \
    git \
    vim \
    nano \
    unzip \
    zip \
    tini \
    # R-base and openssh-server removed as likely not needed
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# --- Create non-root user ---
# Using different UID/GID than default root/staff/users
ARG NB_USER="sandboxuser"
ARG NB_UID=1001
ARG NB_GID=1001
ENV USER=${NB_USER}
ENV HOME=/home/${NB_USER}
ENV PATH=${HOME}/.local/bin:${PATH} 

# Create the group and user, add to sudo (optional)
RUN groupadd -g ${NB_GID} ${NB_USER} && \
    useradd -m -s /bin/bash -u ${NB_UID} -g ${NB_GID} ${NB_USER} && \
    adduser ${NB_USER} sudo && \
    echo "${NB_USER} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers && \
    chown -R ${NB_UID}:${NB_GID} ${HOME}

# Create directories needed for jupyter and set ownership
# Run this as root before switching user
RUN mkdir -p ${HOME}/.local/share/jupyter \
             ${HOME}/.ipython/profile_default/startup \
             ${HOME}/.ipython/profile_default/static \
    && chown -R ${NB_UID}:${NB_GID} ${HOME}/.local \
                                   ${HOME}/.ipython

# Copy requirements *before* installing them, into user's home
COPY ./requirements.txt ${HOME}/requirements.txt
RUN chown ${NB_UID}:${NB_GID} ${HOME}/requirements.txt

# Switch to the non-root user BEFORE installing pip packages
USER ${NB_USER}
WORKDIR ${HOME}

# Install python packages as the non-root user
# Upgrade pip first
# No need for --break-system-packages with official python images
RUN python -m pip install --no-cache-dir --upgrade pip --user
# Install specific versions of core Jupyter components + requirements
RUN python -m pip install --no-cache-dir --user \
    ipython==8.12.0 \
    traitlets==5.9.0 \
    jupyter_client==8.3.0 \
    jupyter_core==5.3.1 \
    pyzmq==25.1.0 \
    tornado==6.3.2 \
    ipykernel==6.25.1 
    # -r ${HOME}/requirements.txt

# Expose the ports Jupyter kernel might use (optional, for info)
EXPOSE 9000-9010

# Define the command to run the kernel and save connection info
# Use 'python' instead of 'python3' as it's the standard in official images
# --- Add right above the CMD (anywhere after USER/WORKDIR) -------------
ENV IPY_BASE_PORT=4000  
# ----------------------------------------------------------------------

# --- REPLACE the old single-line CMD with everything below -------------
# Start ipykernel on fixed ports 4000-4004 so the host can map them.
COPY start_kernel.py /usr/local/bin/start_kernel.py
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["python", "/usr/local/bin/start_kernel.py"]

