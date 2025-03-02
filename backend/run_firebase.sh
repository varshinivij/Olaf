#!/bin/bash

set -e  # Exit script on first error

IMAGE_NAME="firebase_emulator"
CONTAINER_NAME="firebase_emulator"

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "❌ Docker is not installed. Please install Docker and try again."
    exit 1
fi

# Check if Docker is running
if ! docker info &> /dev/null; then
    echo "❌ Docker is not running. Please start Docker and try again."
    exit 1
fi

echo "🔍 Checking if the Docker image is up to date..."

# Check if the image exists
if docker images | grep -q "$IMAGE_NAME"; then
    echo "✅ Image $IMAGE_NAME already exists."
    
    # Check if the Dockerfile or relevant files have changed since last build
    if [ -f .docker_timestamp ] && [ "$(find Dockerfile functions/requirements.txt -newer .docker_timestamp)" = "" ]; then
        echo "🚀 No changes detected in Docker dependencies. Skipping rebuild."
    else
        echo "🔄 Changes detected. Rebuilding the Docker image..."
        docker compose build
        touch .docker_timestamp  # Update timestamp file
    fi
else
    echo "🚀 No existing image found. Building the Docker image..."
    docker compose build
    touch .docker_timestamp  # Create timestamp file after build
fi

# Check if the container is already running
if docker ps | grep -q "$CONTAINER_NAME"; then
    echo "✅ Firebase Emulator is already running."
else
    echo "🚀 Starting Firebase Emulator..."
    docker compose up -d
fi

echo "✅ Firebase Emulator is now running! Use 'docker compose logs -f' to see logs."