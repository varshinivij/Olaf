from openai import OpenAI
import requests
import re
import json

api_key = "REMOVED"
url = "https://api.openai.com/v1/chat/completions"
client = OpenAI(api_key=api_key)

def chat_completion(history, tools=None):

    headers = {
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json"
        }

    payload = {
        "model": "gpt-4o",
        "messages": history.get_history(),
        "temperature": 0.1,
        "stream": True
    }

    response = requests.post(url, headers=headers, json=payload, stream=True)

    if response.status_code == 200:
        for chunk in response.iter_lines(decode_unicode=True):
            if chunk:
                try: 
                    chunk = chunk.replace("data: ", "").strip()
                    json_chunk = json.loads(chunk)
                    yield json_chunk
                except:
                    continue

    else:
        print(f"Error: {response.status_code}, {response.text}")
        return None

    
def chat_completion_function(history, tools=None):
    stream = client.chat.completions.create(
        model="gpt-4o",
        messages=[
            {"role": "system", "content": "Call the function best suited to find the response"},
            *history.get_history(),
        ],
        temperature=0.1,
        tools= tools,
    )

    return stream


def chat_completion_api(history, system_prompt, tools=None):

    headers = {
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json"
        }

    payload = {
        "model": "gpt-4o",
        "messages": [
            {"role": "system", "content": system_prompt},
            *history.get_history(),
        ],
        "temperature": 0.1,
        "tools": tools,
        "stream": True
    }

    response = requests.post(url, headers=headers, json=payload, stream=True)

    if response.status_code == 200:
        for chunk in response.iter_lines(decode_unicode=True):
            if chunk:
                try: 
                    chunk = chunk.replace("data: ", "").strip()
                    json_chunk = json.loads(chunk)
                    yield json_chunk
                except:
                    continue

    else:
        print(f"Error: {response.status_code}, {response.text}")
        return None
    

def chat_completion_plan(history, system_prompt, tools=None):

    headers = {
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json"
        }

    payload = {
        "model": "gpt-4o-2024-08-06",
        "messages": [
            {"role": "system", "content": system_prompt},
            *history.get_history(),
        ],
        "temperature": 0.1,
        "tools": tools,
        "response_format": {
            "type": "json_schema",
            "json_schema": {
                "name": "plan_response",
                "schema": {
                    "type": "object",
                    "properties": {
                        "overview": {
                            "type": "string",
                            "description": "A brief overview of the entire plan."
                        },
                        "steps": {
                            "type": "array",
                            "items": {
                                "type": "object",
                                "properties": {
                                    "step_number": {"type": "integer"},
                                    "title": {"type": "string", "description": "A short title or name for the step."},
                                    "description": {"type": "string", "description": "A detailed explanation of what this step involves."},
                                    "dependencies": {
                                        "type": "array",
                                        "items": {"type": "string"},
                                        "description": "Any prerequisites or dependencies required for this step."
                                    },
                                    "expected_output": {"type": "string", "description": "The expected result or output of this step."}
                                },
                                "required": ["step_number", "title", "description", "dependencies", "expected_output"],
                                "additionalProperties": False
                            }
                        },
                        "final_check": {
                            "type": "string",
                            "description": "A final review or checklist after all steps are completed."
                        }
                    },
                    "required": ["overview", "steps", "final_check"],
                    "additionalProperties": False
                },
                "strict": True
            }
        },
        "stream": True
    }

    response = requests.post(url, headers=headers, json=payload, stream=True)

    if response.status_code == 200:
        for chunk in response.iter_lines(decode_unicode=True):
            if chunk:
                try: 
                    chunk = chunk.replace("data: ", "").strip()
                    json_chunk = json.loads(chunk)
                    yield json_chunk
                except:
                    continue

    else:
        print(f"Error: {response.status_code}, {response.text}")
        return None


def extract_python_code(text):
    """
    Extract the Python code from the provided text string.

    Args:
    text (str): The input text containing the Python code

    Returns:
    str: The extracted Python code
    """
    code_pattern = re.compile(r'```python(.*?)```', re.DOTALL)
    match = code_pattern.search(text)

    if match:
        return match.group(1).strip()
    else:
        return "No Python code found in the input text."
    
def extract_code_and_text(content: str) -> tuple[str, str]:
    code_blocks = re.findall(r'```(.*?)```', content, re.DOTALL)
    text_parts = re.split(r'```.*?```', content, flags=re.DOTALL)
    text = ' '.join([part.strip() for part in text_parts if part.strip()])
    code = '\n'.join([block.strip() for block in code_blocks if block.strip()])
    return text, code
    

def stream(agent):
    for chunk in agent.generate():
        print(chunk)
        try:
            content = chunk['choices'][0]['delta']['content']
            yield content.encode('utf-8')
        except:
            if (chunk.startswith("Response: ")):
                yield chunk.replace("Response: ", "").encode('utf-8')
            continue