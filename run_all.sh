#!/usr/bin/env bash
#
# run_all.sh
#
# Runs both the frontend (Angular dev mode) and backend (Firebase).
# Update commands/path if your project or run commands differ.
echo "=== Starting Olaf Backend (Firebase) ==="
cd backend || exit
chmod +x run_firebase.sh
./run_firebase.sh &

echo "=== Starting Olaf Frontend in Development Mode ==="
cd ../frontend || exit
npm install
npm run dev &  # run in the background

echo "Frontend started in background."
echo
# Once Firebase stops or is killed, the script ends.
# If you want to also kill the frontend automatically upon exit,
# you can track its PID and kill it here.