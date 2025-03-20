#!/usr/bin/env bash
#
# stop_all.sh
#
# Stops running processes for the backend

echo "=== Stopping Olaf Backend (Firebase) ==="
cd /backend || exit
chmod +x stop_firebase
./stop_firebase
kill $(lsof -t -i :4200) 2>/dev/null || echo "Frontend not running on port 4200"