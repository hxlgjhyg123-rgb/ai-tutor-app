# ==============================================================================
# 1. 导入项目所需的全部库
# ==============================================================================
import os
import re
import streamlit as st
import google.generativeai as genai
from Bio import Entrez

# ==============================================================================
# 2. API 配置
# ==============================================================================
try:
    api_key = st.secrets["GEMINI_API_KEY"]
    genai.configure(api_key=api_key)
except Exception:
    st.error("无法配置API密钥！请确保您已在Streamlit Cloud的设置中正确添加了GEMINI_API_KEY。")

# ==============================================================================
# 3. 定义核心函数 (已全面升级)
# ==============================================================================

# --- PubMed搜索函数 (稍作优化) ---
def search_pubmed(queries, max_results=2):
    """根据一个关键词列表，去PubMed搜索并返回格式化的文献信息。"""
    try:
        Entrez.email = "307660349@qq.com" # <-- (重要) 请务必替换为您自己的真实邮箱!
        
        all_id_list = []
        for query in queries:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            all_id_list.extend(record["IdList"])
        
        # 去重
        unique_id_list = list(set(all_id_list))
        if not unique_id_list:
            return "在PubMed中未找到与您输入内容高度相关的文献。", []

        handle = Entrez.efetch(db="pubmed", id=unique_id_list, rettype="medline", retmode="text")
        records_text = handle.read()
        handle.close()
        
        formatted_references = []
        raw_text_for_llm = ""
        papers = records_text.strip().split("\n\n")
        
        for i, paper_text in enumerate(papers):
            try:
                pmid = [line.split("- ")[1] for line in paper_text.split('\n') if line.startswith("PMID-")][0]
                title = [line.split("- ")[1] for line in paper_text.split('\n') if line.startswith("TI  -")][0]
                abstract_lines = [line[6:] for line in paper_text.split('\n') if line.startswith("AB  -")]
                abstract = " ".join(abstract_lines)

                ref_string_for_display = f"**{i+1}. {title}** (PMID: {pmid})"
                formatted_references.append(ref_string_for_display)
                
                raw_text_for_llm += f"--- 参考资料 {i+1} (PMID: {pmid}) ---\n标题: {title}\n摘要: {abstract}\n\n"
            except IndexError:
                continue

        return raw_text_for_llm, formatted_references
    except Exception as e:
        return f"PubMed搜索时发生错误: {e}", []

# --- Gemini API调用函数 (保持不变) ---
def ask_gemini(prompt_text):
    try:
        model = genai.GenerativeModel('gemini-1.5-flash-latest')
        response = model.generate_content(prompt_text, request_options={'timeout': 120}) # 增加超时时间
        return response.text
    except Exception as e:
        return f"调用Gemini API时发生错误: {e}"

# --- 新增：智能生成搜索关键词的函数 ---
def generate_search_queries(role, experiment_type, observation, conclusion, method):
    """第一步AI调用：根据角色和输入，生成PubMed搜索关键词。"""
    
    # 针对不同角色，生成关键词的指令也不同
    query_generation_prompts = {
        "苏格拉底导师 (Socratic Tutor)": f"我是一名学生，正在学习关于“{experiment_type}”的知识，我对“{observation}”这个现象感到困惑。请为我生成2个PubMed搜索关键词，帮助我查找相关的基础原理或背景知识。",
        "魔鬼代言人 (Devil's Advocate)": f"我的结论是“{conclusion}”，方法是“{method}”。请为我生成3个PubMed搜索关键词，用于查找能够挑战我这个结论的文献，比如寻找替代解释、实验方法的缺陷或矛盾的证据。",
        "文献连接者 (Knowledge Connector)": f"我的核心发现是“{conclusion}”。请为我生成3个PubMed搜索关键词，用于查找与这个发现最相关的研究背景和前沿进展。",
        "下一步战略家 (Next-Step Strategist)": f"我的发现是“{conclusion}”。请为我生成3个PubMed搜索关键词，用于查找可以用来验证这个发现的实验方法，或者可以探索的下游机制研究。"
    }
    
    prompt = f"""{query_generation_prompts[role]}
    请只返回关键词列表，用换行符分隔，不要任何多余的解释。例如：
    Keyword 1
    Keyword 2
    """
    
    response_text = ask_gemini(prompt)
    # 从返回的文本中解析出关键词列表
    queries = [line.strip() for line in response_text.split('\n') if line.strip()]
    return queries


# --- 最终的核心引擎函数 (实现两步走逻辑) ---
def get_ai_feedback_with_rag(role, experiment_type, observation, conclusion, method):
    # 1. 第一步：智能生成搜索关键词
    status_text.info("第一步：AI正在思考应该查阅哪些资料...")
    search_queries = generate_search_queries(role, experiment_type, observation, conclusion, method)
    if not search_queries:
        return "无法生成有效的搜索关键词，请尝试更详细地描述您的问题。", []
    
    status_text.info(f"已生成搜索关键词：{', '.join(search_queries)}")

    # 2. 第二步：用生成的关键词去PubMed检索
    status_text.info("第二步：正在PubMed数据库中检索相关文献...")
    retrieved_papers, formatted_references = search_pubmed(search_queries)
    if not formatted_references:
        return "未能在PubMed中找到与AI生成的关键词高度相关的文献。", []
    
    status_text.success(f"已成功检索到 {len(formatted_references)} 篇相关文献！")

    # 3. 第三步：整合所有信息，进行最终的RAG生成
    status_text.info("第三步：AI正在阅读文献并组织最终的反馈...")
    
    final_rag_prompt = f"""
    你是一位严谨的AI科研导师，当前扮演的角色是：**{role}**。

    这是用户的输入信息：
    - 实验领域/类型: {experiment_type}
    - 实验方法: {method}
    - 观察: {observation}
    - 初步结论/核心发现: {conclusion}

    ---
    这是我为你从PubMed检索到的真实文献摘要，请你**严格依据**这些信息来完成你的角色任务。
    {retrieved_papers}
    ---

    请开始你的回答。在你的回答中，当论证关键信息时，必须以 `[PMID: XXXXXX]` 的格式清晰地引用你参考的文献。
    """
    
    feedback = ask_gemini(final_rag_prompt)
    return feedback, formatted_references

# ==============================================================================
# 5. 构建Streamlit用户界面
# ==============================================================================
st.set_page_config(layout="wide", page_title="AI科研思维训练工具 V3.0")
st.title("🔬 AI 科研思维训练工具 (全角色文献引用版)")
st.markdown("本工具的每一个回答都基于实时的PubMed文献检索，确保关键论据有据可查。")

with st.sidebar:
    st.header("⚙️ 输入参数")
    role = st.selectbox(
        "1. 请选择AI导师的角色", 
        ["苏格拉底导师 (Socratic Tutor)", "魔鬼代言人 (Devil's Advocate)", "文献连接者 (Knowledge Connector)", "下一步战略家 (Next-Step Strategist)"]
    )
    experiment_type = st.text_input("2. 实验领域/类型", placeholder="例如：肿瘤免疫治疗, Western Blot")
    method = st.text_area("3. 简述实验方法", height=100)
    observation = st.text_area("4. 描述您的实验观察", height=100)
    conclusion = st.text_area("5. 您的初步结论/核心发现", height=100)
    submit_button = st.button("🚀 获取AI导师的反馈", use_container_width=True)

st.header("💬 AI导师的反馈")

# 使用一个容器来放置状态文本和最终结果
output_container = st.container()

if submit_button:
    # 在主界面显示状态更新
    with output_container:
        status_text = st.empty()
        feedback, references = get_ai_feedback_with_rag(role, experiment_type, observation, conclusion, method)
        
        # 清空状态文本
        status_text.empty()

        # 显示最终结果
        col1, col2 = st.columns([2, 1.5])
        with col1:
            st.subheader("📝 AI生成的分析与洞见")
            st.markdown(feedback)
        if references:
            with col2:
                st.subheader("📚 引用的参考文献 (来自PubMed)")
                for ref in references:
                    st.markdown(f"- {ref}")
else:
    with output_container:
        st.info("请在左侧填写信息并点击按钮，AI的反馈将显示在这里。")