import datetime
import os
from firebase_admin import initialize_app, credentials
from google.auth.credentials import Credentials

# Create a dummy credentials object that fulfills the google.auth.credentials.Credentials interface.
class DummyGoogleCredentials(Credentials):
    def __init__(self):
        super().__init__()
        self.token = "dummy-token"
        self.expiry = datetime.datetime.max

    def refresh(self, request):
        # No need to actually refresh; just keep a dummy token.
        self.token = "dummy-token"
        self.expiry = datetime.datetime.max

    def apply(self, headers, token=None):
        headers["Authorization"] = "Bearer dummy-token"

# Create a dummy credential for Firebase Admin by subclassing its Base.
class DummyCredential(credentials.Base):
    def get_credential(self):
        return DummyGoogleCredentials()

# Determine if we're running in emulator mode
IS_EMULATOR = os.getenv("FUNCTIONS_EMULATOR") == "true"
print(f"Running in emulator: {IS_EMULATOR}")

if IS_EMULATOR:
    print("Running in Firebase Emulator. Connecting to local emulators...")
    # Set environment variables for local emulators
    os.environ["FIRESTORE_EMULATOR_HOST"] = "localhost:8080"
    os.environ["FIREBASE_AUTH_EMULATOR_HOST"] = "localhost:9099"
    os.environ["STORAGE_EMULATOR_HOST"] = "http://localhost:9199"
    os.environ["GOOGLE_CLOUD_PROJECT"] = "twocube-web"

    # Use the dummy credential and explicitly set the project ID.
    dummy_cred = DummyCredential()
    initialize_app(dummy_cred, {"projectId": "twocube-web"})
else:
    initialize_app()

# Expose endpoints to Firebase Functions
from endpoints.agent_endpoints import *
from endpoints.e2b_endpoints import *
from endpoints.file_storage_endpoints import *
from endpoints.session_endpoints import *