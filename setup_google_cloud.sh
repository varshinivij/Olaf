#!/usr/bin/env bash
#
# setup_frontend_google_cloud.sh
#
# Interactive script to configure your Angular frontend environment for
# Google Cloud / Firebase. It creates environment.ts in the specified path
# with user-provided Firebase and API endpoint settings.

# Flavor text
echo "============================================================"
echo "  Welcome to the Olaf Frontend Setup for Google Cloud!"
echo "  This script will configure your environment.ts file"
echo "  so your Angular app connects to a remote Firebase project."
echo "============================================================"
echo

# Prompt user for each Firebase key
read -rp "Enter your Firebase projectId (e.g., twocube-web): " FIREBASE_PROJECT_ID
read -rp "Enter your Firebase appId (e.g., 1:123456:web:abcd1234): " FIREBASE_APP_ID
read -rp "Enter your Firebase storageBucket (e.g., myapp.appspot.com): " FIREBASE_STORAGE_BUCKET
read -rp "Enter your Firebase apiKey: " FIREBASE_API_KEY
read -rp "Enter your Firebase authDomain (e.g., myapp.firebaseapp.com): " FIREBASE_AUTH_DOMAIN
read -rp "Enter your Firebase messagingSenderId (e.g., 123456): " FIREBASE_MESSAGING_SENDER_ID
read -rp "Enter your Firebase measurementId (e.g., G-ABC123XYZ): " FIREBASE_MEASUREMENT_ID

echo
echo "Now let's configure your Production API endpoints..."
read -rp "Enter your chatAPIEndpoint: " CHAT_API_ENDPOINT
read -rp "Enter your l3chatAPIEndpoint: " L3CHAT_API_ENDPOINT
read -rp "Enter your nameMakerAPIEndpoint: " NAME_MAKER_API_ENDPOINT
read -rp "Enter your boxRequestApi: " BOX_REQUEST_API
read -rp "Enter your boxExecuteApi: " BOX_EXECUTE_API
read -rp "Enter your boxCloseApi: " BOX_CLOSE_API
read -rp "Enter your boxStatusApi: " BOX_STATUS_API
read -rp "Enter your boxUploadApi: " BOX_UPLOAD_API
read -rp "Enter your addFirebaseFileApi: " ADD_FIREBASE_FILE_API
read -rp "Enter your boxTerminalApi: " BOX_TERMINAL_API

# Directory where environment.ts should be placed
FRONTEND_ENV_PATH="/frontend/src/environments"
mkdir -p "$FRONTEND_ENV_PATH"

# Write out the environment.ts file
cat <<EOF > "$FRONTEND_ENV_PATH/environment.ts"
export const environment = {
  production: true,
  firebase: {
    projectId: '$FIREBASE_PROJECT_ID',
    appId: '$FIREBASE_APP_ID',
    storageBucket: '$FIREBASE_STORAGE_BUCKET',
    apiKey: '$FIREBASE_API_KEY',
    authDomain: '$FIREBASE_AUTH_DOMAIN',
    messagingSenderId: '$FIREBASE_MESSAGING_SENDER_ID',
    measurementId: '$FIREBASE_MEASUREMENT_ID',
  },

  // Production API Endpoints
  chatAPIEndpoint: '$CHAT_API_ENDPOINT',
  l3chatAPIEndpoint: '$L3CHAT_API_ENDPOINT',
  nameMakerAPIEndpoint: '$NAME_MAKER_API_ENDPOINT',

  // Box-related APIs
  boxRequestApi: '$BOX_REQUEST_API',
  boxExecuteApi: '$BOX_EXECUTE_API',
  boxCloseApi: '$BOX_CLOSE_API',
  boxStatusApi: '$BOX_STATUS_API',
  boxUploadApi: '$BOX_UPLOAD_API',
  addFirebaseFileApi: '$ADD_FIREBASE_FILE_API',
  boxTerminalApi: '$BOX_TERMINAL_API',
};
EOF

echo
echo "============================================================"
echo "Your environment.ts file has been created at:"
echo "  $FRONTEND_ENV_PATH/environment.ts"
echo "============================================================"
echo "Setup complete! You can now build or run your Angular app."