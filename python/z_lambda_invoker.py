import os, json
from z_lambda_batch import handler

# ----- craft an event just like API Gateway would -----
event = {"body": json.dumps({
    "image": "ami-03a152f600ede95ee", # your AMI with /home/ubuntu/jacob/z
    "t": 8888888888         # your test t
})}

# ----- invoke -----
resp = handler(event, None)
print(resp)
