# -------------------------------
# File: src/agent.py
# -------------------------------
import os
os.environ["TOKENIZERS_PARALLELISM"] = "false"
from collections import defaultdict
from summarizer import summarize
from typing import List, Tuple, Optional
from Bio import Entrez
import time

Entrez.email = "oppaboo456@gmail.com"  # Required by NCBI


Entrez.email = "your_email@example.com"


def fetch_from_pubmed(topic: str, max_results: int = 5) -> List[dict]:
    # Step 1: Search PubMed
    search_handle = Entrez.esearch(db="pubmed", term=topic, retmax=max_results)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    ids = search_results["IdList"]

    if not ids:
        return []

    # Step 2: Fetch summaries
    fetch_handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    records = Entrez.read(fetch_handle)
    fetch_handle.close()

    papers = []
    for article in records["PubmedArticle"]:
        try:
            citation = article["MedlineCitation"]
            article_info = citation["Article"]

            title = article_info["ArticleTitle"]

            # Handle multipart abstracts
            abstract_sections = article_info.get("Abstract", {}).get("AbstractText", [])
            abstract = " ".join(
                str(sec) for sec in abstract_sections) if abstract_sections else "No abstract available."

            journal = article_info.get("Journal", {}).get("Title", "Unknown Journal")
            pub_date_info = article_info.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
            pub_year = pub_date_info.get("Year") or pub_date_info.get("MedlineDate", "Unknown Year")

            papers.append({
                "title": title,
                "summary": abstract,
                "journal": journal,
                "year": pub_year,
                "source": "PubMed"
            })
        except Exception as e:
            print(f"Skipping a record due to parsing error: {e}")
            continue

    return papers


def research_topic_basic(topic: str, sources: List[str]) -> str:
    all_papers = []
    for source in sources:
        if "pubmed" in source:
            all_papers.extend(fetch_from_pubmed(topic))
        else:
            print(f"Source '{source}' not yet supported.")

    if not all_papers:
        raise ValueError("No papers found from provided sources.")

    docs = [f"{paper['title']}\n\n{paper['summary']}" for paper in all_papers]
    summary, usage, duration = summarize(docs, topic)
    return summary


def research_topic(topic: str, sources: List[str]) -> Tuple[str, Optional[dict], float]:
    all_papers: list[dict] = []

    for source in sources:
        if "pubmed" in source:
            all_papers.extend(fetch_from_pubmed(topic))
        else:
            print(f"[WARN] Source '{source}' not yet supported.")

    if not all_papers:
        raise ValueError("No papers found from provided sources.")

    # Group papers by source
    docs_by_source = defaultdict(list)
    for paper in all_papers:
        source = paper.get("source", "Unknown")
        title = paper.get("title", "Untitled")
        summary = paper.get("summary", "No summary available.")
        journal = paper.get("journal", None)
        year = paper.get("year", None)

        meta = f"({journal}, {year})" if journal and year else ""
        formatted = f"**{title}** {meta}\n\n{summary}"
        docs_by_source[source].append(formatted)

    final_sections = []
    total_usage = {"prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0}
    start_time = time.time()

    for source, docs in docs_by_source.items():
        print(f"[INFO] Summarizing papers from: {source} ({len(docs)} papers)")
        try:
            summary, usage, _ = summarize(docs, topic)

            # Aggregate token usage if available
            if usage:
                total_usage["prompt_tokens"] += usage.prompt_tokens
                total_usage["completion_tokens"] += usage.completion_tokens
                total_usage["total_tokens"] += usage.total_tokens

            section = f"### Summary from {source}\n\n{summary}"
            final_sections.append(section)

        except Exception as e:
            print(f"[ERROR] Failed to summarize source {source}: {e}")
            continue

    duration = round(time.time() - start_time, 2)
    final_output = "\n\n---\n\n".join(final_sections)

    return final_output, total_usage, duration


if __name__ == "__main__":
    topic = "Alzheimer's disease gene expression"
    sources = ["pubmed"]
    summary = research_topic(topic, sources)
    print(summary)
