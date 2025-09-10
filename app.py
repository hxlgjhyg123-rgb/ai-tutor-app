# -----------------------------------------------------------------
# 1. 导入必要的库
# -----------------------------------------------------------------
import os
import streamlit as st
import google.generativeai as genai

# -----------------------------------------------------------------
# 2. 配置Google Gemini API
#    程序会从Streamlit的Secrets Manager中读取密钥
# -----------------------------------------------------------------
try:
    api_key = st.secrets["GEMINI_API_KEY"]
    genai.configure(api_key=api_key)
except Exception:
    # 如果在本地测试或未配置密钥，会显示此错误
    st.error("无法配置API密钥！请确保您已在Streamlit Cloud的设置中正确添加了GEMINI_API_KEY。")

# -----------------------------------------------------------------
# 3. 定义核心函数和提示词模板
# -----------------------------------------------------------------

# 提示词模板字典，储存四个角色的指令
PROMPT_TEMPLATES = {
    "苏格拉底导师 (Socratic Tutor)": """你是一位经验丰富的分子生物学导师。我是一名研究生。
我的初步观察是：{observation}
请不要直接告诉我答案，而是通过提问的方式，引导我思考并确认我对实验技术、关键质控点和图表信息的理解是否到位。
只有当我回答不上来时，你再为我解释对应的基础知识。""",
    
    "魔鬼代言人 (Devil's Advocate)": """你是一位以严谨和批判性著称的顶刊审稿人。
我的实验方法是：{method}
我的结果与初步结论是：{conclusion}
你的唯一任务是严格审视我的逻辑链条，并针对我的结论，提出至少三个有力的、可能的其他解释或潜在的实验漏洞。请直接开始提出质疑。""",
    
    "文献连接者 (Knowledge Connector)": """你是一位在 {experiment_type} 领域知识渊博的专家。
我的核心发现是：{conclusion}
请基于你的知识库，帮我分析：
1. 这个发现与当前领域内的哪些主流观点或通路是一致的？
2. 这个发现是否与某些已发表的研究结果相矛盾？如果矛盾，可能的原因是什么？
3. 这个发现可能的新颖性(Novelty)在哪里？""",
    
    "下一步战略家 (Next-Step Strategist)": """你是一位经验丰富的实验室PI，擅长设计简洁而关键的实验。
我们目前的情况是：{observation}
我们的结论和遇到的挑战是：{conclusion}
请为我设计两个关键的后续实验来验证或推进我们的发现。请清晰地说明每个实验的核心目的、核心方法和预期的结果解读。"""
}

# 调用Gemini API的函数
def ask_gemini(prompt_text):
  """发送一个prompt给Gemini模型并返回文本回答。"""
  try:
      model = genai.GenerativeModel('gemini-1.5-flash-latest')
      response = model.generate_content(prompt_text)
      return response.text
  except Exception as e:
      # 返回一个格式化的错误信息，方便在界面上显示
      return f"调用API时发生错误: {e}"

# -----------------------------------------------------------------
# 4. 构建Streamlit用户界面 (UI)
# -----------------------------------------------------------------

# 设置页面标题和布局
st.set_page_config(layout="wide", page_title="AI科研思维训练工具")

# 标题和介绍
st.title("🔬 AI 科研思维训练工具")
st.markdown("本工具旨在通过结构化的AI对话，帮助研究生更深入、更批判性地解读实验结果。")

# 创建一个侧边栏用于放置输入控件
with st.sidebar:
    st.header("⚙️ 输入参数")
    
    # 下拉菜单选择AI角色
    role = st.selectbox(
        label="1. 请选择AI导师的角色",
        options=list(PROMPT_TEMPLATES.keys()),
        help="不同的角色会从不同角度对您的输入进行分析和提问。"
    )
    
    # 文本框输入各项信息
    experiment_type = st.text_input(label="2. 实验领域/类型", placeholder="例如：肿瘤免疫治疗, Western Blot")
    method = st.text_area(label="3. 简述实验方法", height=100)
    observation = st.text_area(label="4. 描述您的实验观察", height=100)
    conclusion = st.text_area(label="5. 您的初步结论/遇到的问题", height=100)
    
    # 提交按钮
    submit_button = st.button("🚀 获取AI导师的反馈", use_container_width=True)

# 主界面用于显示输出
st.header("💬 AI导师的反馈")
output_placeholder = st.empty() # 创建一个占位符，用于显示结果或等待信息
output_placeholder.info("请在左侧填写信息并点击按钮，AI的反馈将显示在这里。")

# 当用户点击提交按钮时，执行下面的逻辑
if submit_button:
    # 显示一个加载提示
    with st.spinner('AI导师正在深度思考中，请稍候...'):
        # 1. 根据角色选择，从字典获取对应的Prompt模板
        prompt_template = PROMPT_TEMPLATES[role]
        
        # 2. 将用户在文本框中输入的内容，填充到Prompt模板中
        full_prompt = prompt_template.format(
            observation=observation,
            conclusion=conclusion,
            method=method,
            experiment_type=experiment_type
        )
        
        # 3. 调用核心函数，获取AI的回答
        feedback = ask_gemini(full_prompt)
        
        # 4. 将AI的回答显示在主界面的占位符中
        output_placeholder.markdown(feedback)