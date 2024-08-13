import requests
import re

def chat_completion(history):
        api_key = "REMOVED"
        url = "https://api.openai.com/v1/chat/completions"

        headers = {
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json"
        }

        payload = {
            "model": "gpt-4o",
            "messages": history.get_history(),
            "temperature": 0.1
        } 

        response = requests.post(url, headers=headers, json=payload)

        if response.status_code == 200:
            result = response.json()
            response_message = result['choices'][0]['message']['content']
            return response_message

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