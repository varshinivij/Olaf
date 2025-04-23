# Use a standard Ubuntu base image compatible with ARM
FROM ubuntu:latest

# Set DEBIAN_FRONTEND to noninteractive to prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Update and install tzdata non-interactively first
RUN apt-get update && \
    apt-get install -y tzdata && \
    # Set timezone to UTC non-interactively
    ln -fs /usr/share/zoneinfo/Etc/UTC /etc/localtime && \
    dpkg-reconfigure --frontend noninteractive tzdata && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install gcc/g++ separately (good practice)
RUN apt-get update && apt-get install -y --no-install-recommends \
    pkg-config \
    libhdf5-dev \
    gcc \
    g++ \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install other basic dependencies
# Use --no-install-recommends to potentially reduce image size
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    wget \
    git \
    vim \
    nano \
    unzip \
    zip \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    build-essential \
    openssh-server \
    r-base \
    sudo \
    # Add any other essential packages here
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy requirements *before* installing them
COPY ./requirements.txt /home/requirements.txt

# Create the target directory for data copying
RUN mkdir -p /home/user

RUN python3 -m pip install --no-cache-dir --break-system-packages -r /home/requirements.txt
COPY ./startup.sh /startup.sh
RUN chmod +x /startup.sh

RUN mkdir -p /home/user

EXPOSE 8888
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", "--NotebookApp.token=''", "--NotebookApp.password=''"]
# CMD ["/startup.sh"]

