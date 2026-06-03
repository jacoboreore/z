import os, json
from z_lambda_monitor import handler

# ----- craft an event just like API Gateway would -----
event = {"body": json.dumps({
    "jobId": "46ddde63-a716-4492-9d3e-e59b43e73452"
})}

# ----- invoke -----
resp = handler(event, None)
print(resp)
