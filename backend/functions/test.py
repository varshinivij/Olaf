import requests
from history import History
def call_generate_code(history, language="Python"):
    url = "http://127.0.0.1:5001/twocube-web/us-central1/generate_plan"
    headers = {"Content-Type": "application/json"}
    payload = {
        "history": history,
        "language": language
    }
    print("Hello")
    print(payload)

    response = requests.post(url, headers=headers, json=payload)
    
    print("Hello1")
    
    
    if response.status_code == 200:
        return response.json().get("message")
    else:
        return f"Error: {response.status_code} - {response.text}"

# Example usage
history_data = [
    {"role": "system", "content": "You are a helpful assistant."},
    {"role": "user", "content": "Can you write a python function calculate the average"},
    
]
history = History(history_data)
generated_code = call_generate_code(history_data)
print(generated_code)