# ==============================================================================
# 1. 导入项目所需的全部库
# ==============================================================================
import os
import streamlit as st
import google.generativeai as genai
from Bio import Entrez
from langchain_google_genai import GoogleGenerativeAIEmbeddings
from langchain.vectorstores import FAISS

# ==============================================================================
# 2. API 配置
# ==============================================================================
try:
    api_key = st.secrets["GEMINI_API_KEY"]
    genai.configure(api_key=api_key)
except Exception as e:
    st.error(f"Gemini API密钥配置失败: {e}")

# ==============================================================================
# 3. 定义核心函数 (双知识源检索)
# ==============================================================================

# --- 函数1：从我们自建的知识库(FAISS)中检索 ---
# 使用Streamlit的缓存功能，避免每次都重新加载模型和索引
@st.cache_resource
def load_vector_store():
    try:
        # 注意：这里的模型名称是用于“嵌入”的，和后面用于“生成”的模型可以不同
        embeddings = GoogleGenerativeAIEmbeddings(model="models/text-embedding-004")
        # allow_dangerous_deserialization=True 是加载FAISS索引所必需的
        vector_store = FAISS.load_local("faiss_index_gut_microbiome", embeddings, allow_dangerous_deserialization=True)
        return vector_store
    except Exception as e:
        st.error(f"加载本地知识库失败！请确保'faiss_index_gut_microbiome'文件夹已上传至GitHub仓库。错误: {e}")
        return None

def search_local_kb(query, k=3):
    vector_store = load_vector_store()
    if vector_store:
        results = vector_store.similarity_search(query, k=k)
        return "\n\n".join([doc.page_content for doc in results])
    return ""

# --- 函数2：从PubMed实时检索 ---
def search_pubmed(query, max_results=3):
    # (这个函数和我们之前V2版本中的完全一样)
    try:
        Entrez.email = "your_email@example.com" # <-- (重要) 请务必替换为您自己的真实邮箱!
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]
        if not id_list: return "", []
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        records_text = handle.read()
        handle.close()
        
        formatted_references = []
        raw_text_for_llm = ""
        papers = records_text.strip().split("\n\n")
        for i, paper_text in enumerate(papers):
            try:
                pmid = [line.split("- ")[1] for line in paper_text.split('\n') if line.startswith("PMID-")][0]
                title = [line.split("- ")[1] for line in paper_text.split('\n') if line.startswith("TI  -")][0]
                ref_string_for_display = f"**{i+1}. {title}** (PMID: {pmid})"
                formatted_references.append(ref_string_for_display)
                raw_text_for_llm += f"--- PubMed文献 {i+1} (PMID: {pmid}) ---\n{paper_text}\n\n"
            except IndexError: continue
        return raw_text_for_llm, formatted_references
    except Exception as e: return f"PubMed搜索时发生错误: {e}", []

# --- 函数3：调用Gemini API进行最终生成 ---
def ask_gemini(prompt_text):
    try:
        model = genai.GenerativeModel('gemini-1.5-pro-latest') # 使用更强大的Pro模型
        response = model.generate_content(prompt_text, request_options={'timeout': 180})
        return response.text
    except Exception as e: return f"调用Gemini API时发生错误: {e}"

# ==============================================================================
# 4. 构建Streamlit用户界面
# ==============================================================================
st.set_page_config(layout="wide", page_title="AI科研思维训练工具 V4.0")
st.title("🔬 AI 科研思维训练工具 (双知识源版)")
st.markdown("本工具的回答基于**您的专属知识库**和**实时的PubMed文献**，由Gemini 1.5 Pro模型进行综合分析。")

with st.sidebar:
    st.header("⚙️ 输入参数")
    user_question = st.text_area("1. 请输入您的核心问题或发现陈述", height=150, placeholder="例如：我的实验发现，在小鼠模型中，补充丁酸盐能够显著改善高脂饮食诱导的肠道屏障功能障碍。")
    submit_button = st.button("🚀 获取AI导师的综合分析", use_container_width=True)

st.header("💬 AI导师的综合反馈")
output_container = st.container()

if submit_button and user_question:
    with output_container:
        # 1. 并行检索
        with st.spinner('正在检索您的专属知识库和PubMed...'):
            local_context = search_local_kb(user_question)
            pubmed_context, pubmed_references = search_pubmed(user_question)

        # 2. 整合信息并生成最终Prompt
        with st.spinner('已获取资料，正在由Gemini 1.5 Pro进行综合分析...'):
            final_prompt = f"""
            你是一位顶级的生物医学科研专家，擅长整合多源信息进行严谨的分析。

            **用户的核心问题/发现是：**
            {user_question}

            ---
            **信息源一：来自用户专属知识库的核心内容：**
            {local_context if local_context else "（未在专属知识库中找到直接相关内容）"}
            ---
            **信息源二：来自PubMed的最新相关文献：**
            {pubmed_context if pubmed_context else "（未在PubMed中找到直接相关内容）"}
            ---

            **你的任务：**
            请综合以上所有信息，为用户提供一个全面、深入、且有条理的分析报告。报告应包含：
            1.  **核心观点总结**：直接回答或分析用户提出的问题/发现。
            2.  **证据支持**：你的关键论点必须有依据。如果依据来自用户的知识库，请说明。如果依据来自PubMed，请务必以 `[PMID: XXXXXX]` 的格式引用。
            3.  **洞见与建议**：基于综合分析，提出可能的机制、潜在的矛盾点或下一步的关键实验建议。
            """
            
            # 3. 获取最终反馈
            final_feedback = ask_gemini(final_prompt)

            # 4. 显示结果
            st.subheader("📝 AI生成的综合分析报告")
            st.markdown(final_feedback)
            
            if pubmed_references:
                st.subheader("📚 报告中引用的实时PubMed文献")
                for ref in pubmed_references:
                    st.markdown(f"- {ref}")
else:
    output_container.info("请在左侧输入您的问题或发现，然后点击按钮开始分析。")