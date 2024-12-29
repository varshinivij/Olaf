from firebase_admin import initialize_app
from firebase_functions.https_fn import Request, Response, on_request
from firebase_functions.options import CorsOptions

from routes.master_agent_interaction import master_agent_interaction_handler
from routes.name_maker import name_maker_handler


initialize_app()

# import other modules' Cloud Functions
from functions.endpoints.e2b_endpoints import *
from functions.endpoints.file_storage_endpoints import *
from functions.endpoints.session_endpoints import *

"""
TODO ASK IN MEETING:

1. How should the services be structured, move session util methods there?
2. e2b_functions or executor, which one are we using? executor should just be e2bservice?
   + e2b_functions have syntax errors like sandbox_status
3. Separate folder for endpoints?
4. Main idea:
    - endpoints folder for storing actual endpoints
    - any util methods used to implement endpoints should be in services folder
    - maybe rename some files (ie. executor.py) to be services
        - agent-utils could be agent-service? including storing the api keys.
        - later on it would be easy to put in logic to load the API keys from
          some secret store.
    - maybe wrap services in a class? in case we want to stick in API keys for
      example, lets us do work upon loading a service and increasing modularity
      basically like Angular.
    - should a "type" value be put in the messages of the history class?
        - new role called "system"?
"""


@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def master_agent_interaction(req: Request) -> Response:
    return master_agent_interaction_handler(req)


# TODO move this into /routes and treat it as a separate route
@on_request(cors=CorsOptions(cors_origins="*", cors_methods=["post"]))
def name_maker(req: Request) -> Response:
    return name_maker_handler(req)
