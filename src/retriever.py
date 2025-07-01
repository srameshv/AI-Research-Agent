# -------------------------------
# File: src/retriever.py
# -------------------------------
from llama_index.core import SimpleDirectoryReader, VectorStoreIndex
from llama_index.embeddings.huggingface import HuggingFaceEmbedding

def retrieve_documents(topic):
    documents = SimpleDirectoryReader(input_dir="data/sample_docs").load_data()
    embed_model = HuggingFaceEmbedding(model_name="all-MiniLM-L6-v2")
    index = VectorStoreIndex.from_documents(documents, embed_model=embed_model)
    query_engine = index.as_query_engine()
    response = query_engine.query(topic)
    return [response.response]