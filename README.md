# AI-Research-Agent
This project builds a simple research assistant that takes a topic, searches the web/papers, retrieves relevant documents, and summarizes findings using an LLM (OpenAI/Gemini).


## What It Does
- Takes a research topic
- Scans trusted sources (arxiv.org, pubmed, more coming)
- Extracts & filters relevant papers
- Summarizes them automatically, into a human-friendly format
You just describe what you want to learn — it handles the rest.

## Example Usage

from agent import research_topic
topic = "Sinus Infections"
sources = ["arxiv.org", "pubmed"]
summary = research_topic(topic, sources)
print(summary)

## Built With
- tenacity for fault-tolerant retries

- feedparser, biopython for data retrieval

- Summarization layer (plug in OpenAI/GPT, extractive models, etc.)

## Roadmap
- Add support for PapersWithCode and Semantic Scholar
- Integrate sentence-transformer ranking
- Let the agent explain why it picked each paper
- Use LangChain-style memory + context for follow-ups
- Surface metadata: citations, impact factor, etc.
