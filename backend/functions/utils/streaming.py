import json
from typing import Any, Dict, Generator, Callable, Union

def stream_llm_response(
    api_response: Generator, 
    handle_function: Callable[[str, Dict[str, Any,], str], Union[str, Dict[str, Any], Generator]],
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
    function_argument = ""
    current_chunk_type = "text"
    function_called = False
    func_name = ""
    func_args = {}
    for chunk in api_response:
        try:
            # In typical usage: chunk["choices"][0]["delta"]
            choice = chunk["choices"][0]
            delta = choice.get("delta", {})

            # Check for function call
            if "tool_calls" in delta and (not function_called):
                print("Function call detected")
                print(delta)
                func_name = delta["tool_calls"][0]["function"]["name"] #TODO handle all tool calls and don't assume first
                args_str = delta["tool_calls"][0]["function"].get("arguments", "{}")
                try:
                    func_args = json.loads(args_str)
                except json.JSONDecodeError:
                    func_args = {}
                function_called = True
            if function_called:
                argument = delta["tool_calls"][0]["function"].get("argument", "")
                if argument:
                    if argument.startswith("```"):
                        current_chunk_type = "code"
                    elif argument.endswith("```"):
                        current_chunk_type = "text"
                    function_argument += argument
                    content_accumulated += argument
                    yield {"type": current_chunk_type, "content": argument}
            else:
                print("no function call detected")
                # Normal text chunk
                text = delta.get("content", "")
                print(delta)
                print(text)
                if text:
                    # Identify code vs. text
                    if text.startswith("```"):
                        current_chunk_type = "code"
                    elif text.endswith("```"):
                        current_chunk_type = "text"

                    content_accumulated += text
                    yield {"type": current_chunk_type, "content": text}

        except KeyError as e:
            print(f"Error processing chunk: {e}")
            # If the chunk doesn't have the expected structure, skip
            continue
        # The handler might return a generator or a string/dict
    if(function_called):
        print(f"Calling function: {func_name} with args: {func_args}")
        response = handle_function(func_name, func_args, function_argument)
        print(response)
        if hasattr(response, "__iter__") and not isinstance(response, (str, dict)):
            # It's a generator - yield each sub-chunk
            print("response is a generator")
            for sub_chunk in response:
                # Accumulate
                if "content" in sub_chunk:
                    content_accumulated += sub_chunk["content"]
                print(sub_chunk)
                yield sub_chunk
        elif isinstance(response, str):
            # It's a single string - yield once
            yield {"type": "text", "content": response}
            content_accumulated += response
        elif isinstance(response, dict):
            # Possibly routing instructions or other data
            yield {"type": "function_router", "content": response} # type: ignore
        function_called = False
    # Once done, store the entire accumulated text
    store_interaction(role, content_accumulated)