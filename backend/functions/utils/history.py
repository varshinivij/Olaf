import json
class History:
    def __init__(self, system):
        if isinstance(system, str):
            self.history = json.loads(system)[0]
        else:
            self.history = system

    def log(self, role, content):
        entry = {"role": role, "content": content}
        self.history.append(entry)

    def get_history(self):
        return self.history

    def most_recent_entry(self):
        return self.history[-1]

    def remove_system_messages(self):
        self.history = [entry for entry in self.history if entry["role"] != "system"]

    def upsert(self, prompt):
        self.history.insert(0, {"role": "system", "content": prompt})

