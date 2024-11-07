class Pipe:
    def __init__(self, name: str, operations: list = None):
        self.name = name
        self.operations = operations if operations else []
        
    async def process(self, session: dict) -> dict:
        processed_session = session.copy()
        for operation in self.operations:
            try:
                processed_session = await operation(processed_session)
            except Exception as e:
                print(f"Error in pipe {self.name}: {str(e)}")
                raise
        return processed_session