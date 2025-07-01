# -------------------------------
# File: src/agent.py
# -------------------------------
import os
os.environ["TOKENIZERS_PARALLELISM"] = "false"

from retriever import retrieve_documents
from summarizer import summarize

def research_topic(topic):
    docs = retrieve_documents(topic)
    summary = summarize(docs, topic)
    return summary

if __name__ == "__main__":
    topic = input("Enter a research topic: ")
    print(research_topic(topic))