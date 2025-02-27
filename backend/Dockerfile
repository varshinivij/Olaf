# Use a lightweight Python base image
FROM python:3.11-slim

# Fix GPG Errors Before Running apt-get update
RUN apt-get update --allow-releaseinfo-change && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    gpg \
    gnupg \
    curl \
    && rm -rf /var/lib/apt/lists/*  # Force fresh package lists

# Install required system dependencies: Java (for Firebase Emulators)
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre \
    && apt-get clean

# Install Node.js 18 (for Firebase)
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - && \
    apt-get install -y --no-install-recommends nodejs

# Verify installations
RUN node -v && npm -v && python --version && pip --version && java -version

# Install firebase-tools globally
RUN npm install -g firebase-tools

# Set working directory
WORKDIR /usr/src/app

# Copy over Python dependencies
COPY functions/requirements.txt .

# Install Python dependencies (including Functions Framework)
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the project files
COPY . .

# Expose Firebase Emulator ports
EXPOSE 5001 8080 9099 5000

# Start Firebase Emulators
CMD ["firebase", "emulators:start", "--project", "twocube-web"]