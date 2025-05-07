#!/usr/bin/env bash

set -e  # Exit on error

# Process IDs for cleanup
FRONTEND_PID=""
BACKEND_PID=""

# Cleanup function to kill both frontend and backend
cleanup() {
    echo "Cleaning up..."
    if [ -n "$FRONTEND_PID" ] && ps -p $FRONTEND_PID > /dev/null; then
        echo "Killing frontend (PID $FRONTEND_PID)..."
        kill $FRONTEND_PID
    fi
    if [ -n "$BACKEND_PID" ] && ps -p $BACKEND_PID > /dev/null; then
        echo "Killing backend (PID $BACKEND_PID)..."
        kill $BACKEND_PID
    fi
    wait
    echo "Shutdown complete."
}

# Trap errors and signals to run cleanup
trap cleanup EXIT SIGINT SIGTERM ERR

# Create dummy environment.ts if missing
ENV_FILE="/Users/riffled/Desktop/Olaf/frontend/src/environments/environment.ts"
if [ ! -f "$ENV_FILE" ]; then
    echo "environment.ts not found. Creating dummy version..."
    mkdir -p "$(dirname "$ENV_FILE")"
    cat <<EOT > "$ENV_FILE"
export const environment = {
  dummy: true,
};
EOT
fi

# Start backend
echo "=== Starting Olaf Backend (Firebase) ==="
cd backend || exit
chmod +x run_firebase.sh
./run_firebase.sh
BACKEND_PID=$!
echo "Backend started (PID $BACKEND_PID)"

# Start frontend
echo "=== Starting Olaf Frontend in Development Mode ==="
cd ../frontend || exit
npm install
npm run dev &
FRONTEND_PID=$!
echo "Frontend started (PID $FRONTEND_PID)"

# Wait for backend to finish (or be killed)
wait $BACKEND_PID