import json
from firebase_functions.https_fn import Request, Response
from utils.validation import validate_name_maker_request
from functions.services.agent_service import chat_completion
from functions.models.history import History

def name_maker_handler(req: Request) -> Response:
    try:
        history = validate_name_maker_request(req)

        system = """
        You are given a conversation history.
        Please construct a super short title for the conversation.
        """
        hist = History(history)
        hist.log("system", system, "text")

        response = chat_completion(hist)
        response_str = ""
        for chunk in response:
            response_str += chunk.get("choices", [{}])[0].get("delta", {}).get("content", "")

        response_data = {"message": response_str}
        return Response(json.dumps(response_data), status=200)
    except Exception as e:
        return Response(json.dumps({"error": str(e)}), status=500)
