import streamlit as st
import time
from agent import research_topic

st.set_page_config(page_title="AI Research Agent", layout="centered")
st.title("AI Research Agent")
st.markdown("Ask a research question. I'll retrieve and summarize relevant info from online trusted sources.")
topic = st.text_input("Enter a research topic")
sources = st.multiselect("Select sources", ["arxiv.org", "pubmed"], default=["pubmed"])

if st.button("Summarize"):
    if not topic:
        st.warning("Please enter a topic.")
    else:
        with st.spinner("Generating summary..."):
            start_time = time.time()
            formatted, summary, usage, raw_duration = research_topic(topic, sources)
            duration = round(raw_duration, 2)


        st.markdown("## Summary")
        st.write(summary)

        for block in formatted:
            st.markdown(block)
            st.markdown("---")

        # Copy to clipboard
        st.code(summary, language="markdown")

        # Download as Markdown
        st.download_button(
            label="Download Summary",
            data=summary,
            file_name="summary.md",
            mime="text/markdown"
        )

        # Time taken
        st.markdown(f" Time taken: `{duration}` seconds")

        # Actual token usage (if available)
        # Prompt Tokens: These are the tokens that make up your input to the LLM. 
        # Completion Tokens: These are the tokens that the LLM generates as its response to your prompt.
        # Total Tokens: This is the sum of the prompt tokens and the completion tokens for a single API request.
       

        # Rough token estimate
        token_estimate = len(summary.split()) / 0.75  # avg 0.75 words/token
        st.markdown(f"Estimated tokens: `{int(token_estimate)}`")
