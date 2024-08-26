import requests
from history import History

def call_generate_code(history, language="Python"):
    url = "http://127.0.0.1:5001/twocube-web/us-central1/generate_plan"
    headers = {"Content-Type": "application/json"}
    payload = {
        "history": history,
        "language": language
    }
    
    try:
        print("Sending request...")
        response = requests.post(url, headers=headers, json=payload)
        print("Request sent.")

        if response.status_code == 200:
            return response.json().get("message")
        else:
            return f"Error: {response.status_code} - {response.text}"

    except requests.RequestException as e:
        return f"Request failed: {str(e)}"

# Example usage
history_data = [
    {"role": "system", "content": "You are a helpful assistant."},
    {"role": "user", "content": "I have uploaded 16 files, 4 files per cell-type. For each cell type there are three negative sequence files and one positive sequence file. Build a convolution neural network based model to classify positive and negative DNA sequences. For evaluation results, plot the area under precision recall curve and area under the receiver operator characteristic curve"},
]

history = History(history_data)
generated_code = call_generate_code(history_data)
print(generated_code)