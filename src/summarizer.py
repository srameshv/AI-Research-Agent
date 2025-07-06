# -------------------------------
# File: src/summarizer.py
# -------------------------------
import json
import os, time

from datetime import datetime
from dotenv import load_dotenv
from tenacity import retry, retry_if_exception_type
from cache import llm
from langchain_openai import OpenAI

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
    '''
    docs = f"**{title}** {meta}\n\n{summary}"
    '''
    joined = "\n".join(docs)
    prompt = f"Summarize the following research excerpts:\n\n{joined}"

    # Try to fetch from GPTCache
    start = time.time()
    '''
    role: system => Provides high-level instructions or context-setting messages to guide the model's behavior.
    role: user => Represents the messages or queries from the user.
    In future if we want to move to the memory save for context , add this 
    thread_id=response.id

    
    response = client.responses.create(
        model="gpt-3.5-turbo",
        input=[
            {"role": "system", "content": "Summarize the following research excerpts:"},
            {"role": "user", "content": joined}
        ]
    )
    '''
    summary = llm.predict(prompt)  # Auto-cached by LangChain
    end = time.time()

    # summary = response.output[0].content[0].text
    # Save to GPTCache

    usage = {
        "input_tokens": "N/A",
        "output_tokens": "N/A",
        "total_tokens": "N/A"
    }
    duration = round(end - start, 2)

    # Log to file
    log_entry = {
        "timestamp": datetime.utcnow().isoformat(),
        "topic": topic,
        "duration_seconds": duration,
        "summary_preview": summary[:200] + ("..." if len(summary) > 200 else "")
    }
    # Create Log file directory if it doesn't exist.
    os.makedirs("logs", exist_ok=True)
    with open("logs/summary_log.jsonl", "a") as f:
        f.write(json.dumps(log_entry) + "\n")

    return summary, usage, duration
