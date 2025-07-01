import streamlit as st
import time
from agent import research_topic

st.set_page_config(page_title="AI Research Agent", layout="centered")
st.title("🔍 AI Research Agent")
st.markdown("Ask a research question. I'll retrieve and summarize relevant info from your local document store.")

topic = st.text_input("Enter a research topic")

if st.button("Summarize"):
    if not topic:
        st.warning("Please enter a topic.")
    else:
        with st.spinner("Generating summary..."):
            start_time = time.time()
            summary = research_topic(topic)
            duration = round(time.time() - start_time, 2)

        st.markdown("### 🧠 Summary")
        st.write(summary)

        # ✅ Copy to clipboard
        st.code(summary, language="markdown")

        # 📥 Download as Markdown
        st.download_button(
            label="📥 Download Summary",
            data=summary,
            file_name="summary.md",
            mime="text/markdown"
        )

        # ⏱️ Response time
        st.markdown(f"⏱️ Time taken: `{duration}` seconds")

        # 📊 Rough token estimate
        token_estimate = len(summary.split()) / 0.75  # avg 0.75 words/token
        st.markdown(f"📊 Estimated tokens: `{int(token_estimate)}`")
