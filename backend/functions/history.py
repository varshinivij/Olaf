import json
class History:
    def __init__(self, system):
        if isinstance(system, str):
            self.history = json.loads(system)[0]
        else:
            self.history = system
    
    def log(self, role, content):
        entry = {
            "role": role,
            "content": content
        }
        self.history.append(entry)
    
    def get_history(self):
        return self.history
    
    def most_recent_entry(self):
        return self.history[-1]
    
    

# # Example usage:
# history = History()

# # Logging text content
# history.log("user", {"type": "text", "text": "What's in this image?"})

# # Logging an image
# history.log("user", {
#     "type": "image_url",
#     "image_url": {
#         "url": "https://upload.wikimedia.org/wikipedia/commons/thumb/d/dd/Gfp-wisconsin-madison-the-nature-boardwalk.jpg/2560px-Gfp-wisconsin-madison-the-nature-boardwalk.jpg"
#     }
# })

# # Logging a function call
# history.log("user", {
#     "type": "function",
#     "function": {
#         "name": "get_current_weather",
#         "description": "Get the current weather in a given location",
#         "parameters": {
#             "type": "object",
#             "properties": {
#                 "location": {
#                     "type": "string",
#                     "description": "The city and state, e.g. San Francisco, CA"
#                 },
#                 "unit": {
#                     "type": "string",
#                     "enum": ["celsius", "fahrenheit"]
#                 }
#             },
#             "required": ["location"]
#         }
#     }
# })