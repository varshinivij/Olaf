#!/usr/bin/env bash
#
# run_all_with_google_cloud.sh
#
# Runs the frontend in production mode (possibly connecting to real GCP),
# then starts Firebase in the background.

echo "=== Starting Olaf Frontend in Production Mode ==="
cd /frontend || exit
npm install
npm run prod &  # run in the background

echo "Frontend started in background."
echo

echo "=== Starting Olaf Backend (Firebase) ==="
cd /backend || exit
chmod +x run_firebase
./run_firebase