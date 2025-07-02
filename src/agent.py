# -------------------------------
# File: src/agent.py
# -------------------------------
import os
os.environ["TOKENIZERS_PARALLELISM"] = "false"

from retriever import retrieve_documents
from summarizer import summarize
from typing import List


from Bio import Entrez
Entrez.email = "oppaboo456@gmail.com"  # Required by NCBI

def fetch_from_pubmed(topic: str, max_results: int = 5) -> List[dict]:
    # Step 1: Search PubMed for article IDs retmax is a parameter used with the Bio.Entrez.esearch function to
    # specify the maximum number of unique identifiers (UIDs) to be returned in the search results.
    handle = Entrez.esearch(db="pubmed", term=topic, retmax=max_results)
    record = Entrez.read(handle)
    ids = record["IdList"]
    handle.close()

    if not ids:
        return []

    # Step 2: Fetch summaries for those IDs
    summaries = []
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    for article in records["PubmedArticle"]:
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", [""])[0]
        summaries.append({
            "title": title,
            "summary": abstract,
            "source": "PubMed"
        })

    return summaries


def research_topic(topic: str, sources: List[str]) -> str:
    all_papers = []
    for source in sources:
        if "pubmed" in source:
            all_papers.extend(fetch_from_pubmed(topic))
        else:
            print(f"Source '{source}' not yet supported.")

    if not all_papers:
        raise ValueError("No papers found from provided sources.")

    docs = [f"{paper['title']}\n\n{paper['summary']}" for paper in all_papers]
    summary = summarize(docs, topic)
    return summary


if __name__ == "__main__":
    topic = "Alzheimer's disease gene expression"
    sources = ["pubmed"]
    summary = research_topic(topic, sources)
    print(summary)
