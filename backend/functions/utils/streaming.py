import json
from typing import Any, Dict, Generator, Callable, Union

def stream_llm_response(
    api_response: Generator, 
    handle_function: Callable[[str, Dict[str, Any]], Union[str, Dict[str, Any], Generator]],
    store_interaction: Callable[[str, str], None],
    role: str = "assistant"
) -> Generator[Dict[str, str], None, None]:
    """
    A reusable generator that yields partial text chunks or handles function calls 
    from an LLM's streaming response.

    Args:
        api_response: The generator returned by chat_completion_api (stream of tokens).
        handle_function: A function that takes (function_name, arguments) 
                         and returns either a string, dict, or a generator of dict chunks.
        store_interaction: A callback to store the final response in conversation history.
        role: The role to store at the end, e.g., "assistant".

    Yields:
        A dictionary with at least:
            {"type": "text"|"code"|"function_router"|..., "content": "..."}
    """
    content_accumulated = ""
    current_chunk_type = "text"

    for chunk in api_response:
        try:
            # In typical usage: chunk["choices"][0]["delta"]
            choice = chunk["choices"][0]
            delta = choice.get("delta", {})

            # Check for function call
            if "function_call" in delta:
                func_name = delta["function_call"]["name"]
                args_str = delta["function_call"].get("arguments", "{}")
                try:
                    func_args = json.loads(args_str)
                except json.JSONDecodeError:
                    func_args = {}

                # The handler might return a generator or a string/dict
                response = handle_function(func_name, func_args)

                if hasattr(response, "__iter__") and not isinstance(response, (str, dict)):
                    # It's a generator - yield each sub-chunk
                    for sub_chunk in response:
                        yield sub_chunk
                        # Accumulate
                        if "content" in sub_chunk:
                            content_accumulated += sub_chunk["content"]
                elif isinstance(response, str):
                    # It's a single string - yield once
                    yield {"type": "text", "content": response}
                    content_accumulated += response
                elif isinstance(response, dict):
                    # Possibly routing instructions or other data
                    yield {"type": "function_router", "content": response} # type: ignore

            else:
                # Normal text chunk
                text = delta.get("content", "")
                if text:
                    # Identify code vs. text
                    if text.startswith("```"):
                        current_chunk_type = "code"
                    elif text.endswith("```"):
                        current_chunk_type = "text"

                    content_accumulated += text
                    yield {"type": current_chunk_type, "content": text}

        except KeyError:
            # If the chunk doesn't have the expected structure, skip
            continue

    # Once done, store the entire accumulated text
    store_interaction(role, content_accumulated)