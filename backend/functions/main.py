from firebase_admin import initialize_app


initialize_app()


# expose endpoints to Firebase Functions
from functions.endpoints.agent_endpoints import *
from functions.endpoints.e2b_endpoints import *
from functions.endpoints.file_storage_endpoints import *
from functions.endpoints.session_endpoints import *
