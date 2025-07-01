# -------------------------------
# File: src/summarizer.py
# -------------------------------
from openai import OpenAI
import os
from dotenv import load_dotenv
load_dotenv()

api_key = os.environ.get("OPENAI_API_KEY")
org_id = os.environ.get("OPENAI_ORG_ID")

# Automatically include org ID only if it's set and valid
if org_id:
    client = OpenAI(api_key=api_key, organization=org_id)
else:
    client = OpenAI(api_key=api_key)

def summarize(docs):
    joined = "\n".join(docs)
    response = client.chat.completions.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "system", "content": "Summarize the following research excerpts:"},
            {"role": "user", "content": joined}
        ]
    )
    return response.choices[0].message.content