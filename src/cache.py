# src/cache.py

from langchain.globals import set_llm_cache
from langchain_community.cache import SQLiteCache
from langchain_openai import OpenAI
import os
from dotenv import load_dotenv
load_dotenv()
# Set up SQLite-based prompt/response cache
set_llm_cache(SQLiteCache(database_path=".llm_cache.db"))
api_key = os.environ.get("OPENAI_API_KEY")

# Create a single shared LLM instance (caching is automatic)
# use text completion model for now
llm = OpenAI(
    model="gpt-3.5-turbo-instruct",
    api_key=os.environ.get("OPENAI_API_KEY"),
)
