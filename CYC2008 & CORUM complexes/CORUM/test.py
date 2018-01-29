import json

with open("allComplexes.json", "r") as f:
    data = f.read()

data = json.loads(data)

print(json.dumps(data, indent=4))
