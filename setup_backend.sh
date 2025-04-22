#!/usr/bin/env bash
#
# setup_backend.sh
#
# Interactive script to configure your backend environment variables.
# Creates a .env file in the specified path with user-provided keys.

echo "============================================================"
echo "  Welcome to the Olaf Backend Setup!"
echo "  This script will create a .env file in your backend folder"
echo "  with the API keys you provide."
echo "============================================================"
echo

# Prompt user for keys
read -rp "Enter your E2B_API_KEY: " E2B_API_KEY
read -rp "Enter your E2B_TEMPLATE ID (default: vh7kehbtf0t4xbx9ec9u, OpenTechBio's original template): " E2B_TEMPLATE
E2B_TEMPLATE=${E2B_TEMPLATE:-vh7kehbtf0t4xbx9ec9u}
read -rp "Enter your OPENAI_API_KEY (e.g., sk-...): " OPENAI_API_KEY

BACKEND_ENV_PATH="backend/functions"
mkdir -p "$BACKEND_ENV_PATH"

# Write out the .env file
cat <<EOF > "$BACKEND_ENV_PATH/.env"
# This is a template environment file for your backend (Firebase Functions).
# Replace placeholders with your real keys.

E2B_API_KEY=${E2B_API_KEY}
E2B_TEMPLATE=${E2B_TEMPLATE}
OPENAI_API_KEY=${OPENAI_API_KEY}
EOF

echo
echo "============================================================"
echo "Your .env file has been created at:"
echo "  $BACKEND_ENV_PATH/.env"
echo "============================================================"
echo "Setup complete! You can now deploy or run your Firebase functions."
echo