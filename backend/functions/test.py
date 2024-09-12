import requests
from history import History
import json

def call_generate_code(history, language="Python"):
    url = "http://127.0.0.1:5001/twocube-web/us-central1/generate_plan"
    headers = {"Content-Type": "application/json"}
    print(history)
    payload = {
        "history": history,
        "language": language
    }

    with requests.post(url, headers=headers, json=payload, stream=True) as response:
        result = ""
        response.raise_for_status()
        for line in response:
            line = line.decode('utf-8')
            if line:
                result += line

        print(result)

# Example usage
history_data = [
    {"role": "user", "content": "I have uploaded 16 files, 4 files per cell-type. For each cell type there are three negative sequence files and one positive sequence file. Build a convolution neural network based model to classify positive and negative DNA sequences. For evaluation results, plot the area under precision recall curve and area under the receiver operator characteristic curve"}
]

generated_code = call_generate_code(history_data)


# I have uploaded 16 files, 4 files per cell-type. For each cell type there are three negative sequence files and one positive sequence file. Build a convolution neural network based model to classify positive and negative DNA sequences. For evaluation results, plot the area under precision recall curve and area under the receiver operator characteristic curve


# # Define the API key and the endpoint URL
# api_key = 'REMOVED'
# url = 'https://api.openai.com/v1/chat/completions'

# # Set up the headers with your API key for authentication
# headers = {
#     'Authorization': f'Bearer {api_key}',
#     'Content-Type': 'application/json'
# }

# # Define the data for the POST request including the model and prompt
# data = {
#     "model": "gpt-4o",
#     "messages": [{"role": "user", "content": "Say this is a test"}],
#     "max_tokens": 7,
#     "temperature": 0,
#     "stream": True  # Note: `requests` does not support real-time streaming as `curl` does
# }

# # Make the POST request to the OpenAI API
# response = requests.post(url, headers=headers, json=data)

# # Check if the request was successful and print the response
# if response.status_code == 200:
#     try:
#         print(response.json())
#     except ValueError as e:
#         # Handle JSON decoding error (empty response or malformed JSON)
#         print("JSON decoding failed:", e)
#         print("Response content:", response.text)
# else:
#     print("Failed to fetch data:", response.status_code, response.text)