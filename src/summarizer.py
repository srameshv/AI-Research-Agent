# -------------------------------
# File: src/summarizer.py
# -------------------------------
from openai import OpenAI
import os
from dotenv import load_dotenv
from tenacity import retry, retry_if_exception_type
import time, json
from datetime import datetime

load_dotenv()

api_key = os.environ.get("OPENAI_API_KEY")
org_id = os.environ.get("OPENAI_ORG_ID")

# Automatically include org ID only if it's set and valid
if org_id:
    client = OpenAI(api_key=api_key, organization=org_id)
else:
    client = OpenAI(api_key=api_key)


@retry(retry=retry_if_exception_type(ConnectionError))
def summarize(docs, topic):
    joined = "\n".join(docs)
    start = time.time()
    response = client.chat.completions.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "system", "content": "Summarize the following research excerpts:"},
            {"role": "user", "content": joined}
        ]
    )
    end = time.time()

    summary = response.choices[0].message.content
    usage = response.usage
    duration = round(end - start, 2)
    # Log to file
    log_entry = {
        "timestamp": datetime.utcnow().isoformat(),
        "topic": topic,
        "prompt_tokens": usage.prompt_tokens,
        "completion_tokens": usage.completion_tokens,
        "total_tokens": usage.total_tokens,
        "duration_seconds": duration,
        "summary_preview": summary[:200] + ("..." if len(summary) > 200 else "")
    }
    os.makedirs("logs", exist_ok=True)
    with open("logs/summary_log.jsonl", "a") as f:
        f.write(json.dumps(log_entry) + "\n")

    return summary, usage, duration
