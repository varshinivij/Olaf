from firebase_admin import initialize_app


initialize_app()


# expose endpoints to Firebase Functions
from endpoints.agent_endpoints import *
from endpoints.e2b_endpoints import *
from endpoints.file_storage_endpoints import *
from endpoints.session_endpoints import *
