import streamlit as st
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.ExPASy import get_sprot_raw
from Bio.SwissProt import read
import requests
import re

# 设置页面标题
st.set_page_config(page_title="重组蛋白质性质计算器-测试版1", page_icon="🧬", layout="wide")

# 标题和说明
st.title("🧬 蛋白质性质计算器")
st.markdown("""
这个工具可以计算蛋白质的各种物理化学性质，包括：
- **分子量 (kD)** - 蛋白质的分子量
- **等电点 (pI)** - 蛋白质净电荷为零时的pH值
- **消光系数** - 蛋白质在280nm处的摩尔消光系数
- **不稳定指数** - 预测蛋白质在试管中的稳定性
- **GRAVY** - 疏水性平均值
""")

# 定义标签序列
TAG_SEQUENCES = {
    "10his": "HHHHHHHHHH",
    "6his": "HHHHHH",
    "GST": "MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHNMLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALDVVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK",
    "SUMO": "MADLYKQGGKSEVHLTQLHNDLPSLPSPSTVINGLKSKIQTNQKQYSPSVQEAKPEVKPEVKPETHINLKVSDGSSEIFFKIKKTTPLRRLMEAFAKRQGKEMDSLRFLYDGIRIQADQTPEDLDMEDNDIIEAHREQIGG"
}

def search_uniprot_id(protein_name):
    """通过蛋白名称搜索对应的 UniProt ID"""
    try:
        # 清理蛋白名称，移除特殊字符
        clean_name = re.sub(r'[^\w\s-]', '', protein_name).strip()
        
        # 使用 UniProt 搜索 API
        url = f"https://rest.uniprot.org/uniprotkb/search?query={clean_name}+AND+(reviewed:true)&format=tsv&fields=accession,protein_name&size=1"
        response = requests.get(url)
        
        if response.status_code == 200 and response.text.strip():
            lines = response.text.strip().split('\n')
            if len(lines) > 1:  # 有结果
                uniprot_id = lines[1].split('\t')[0]
                return uniprot_id
        return None
    except Exception as e:
        st.error(f"搜索 UniProt ID 时出错: {e}")
        return None

def get_protein_sequence(protein_identifier):
    """根据 UniProt ID 或蛋白名称获取蛋白序列"""
    # 首先检查是否是有效的 UniProt ID 格式
    if re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', protein_identifier):
        # 看起来像 UniProt ID，直接尝试获取
        uniprot_id = protein_identifier
    else:
        # 是通用名称，需要先搜索获取 UniProt ID
        with st.spinner(f"正在搜索蛋白 '{protein_identifier}' 的 UniProt ID..."):
            uniprot_id = search_uniprot_id(protein_identifier)
        if not uniprot_id:
            st.error(f"未找到蛋白 '{protein_identifier}' 的 UniProt ID")
            return None, None
        st.success(f"找到 UniProt ID: {uniprot_id}")
    
    try:
        # 方法1: 通过 ExPASy 获取
        handle = get_sprot_raw(uniprot_id)
        record = read(handle)
        sequence = record.sequence
        return sequence, uniprot_id
    except:
        try:
            # 方法2: 通过 UniProt API 获取
            url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
            response = requests.get(url)
            if response.status_code == 200:
                # 解析FASTA格式，跳过第一行（描述行）
                lines = response.text.strip().split('\n')
                sequence = ''.join(lines[1:])  # 合并所有序列行
                return sequence, uniprot_id
            else:
                return None, None
        except Exception as e:
            st.error(f"获取序列时出错: {e}")
            return None, None

def parse_truncation_range(truncation_text):
    """解析截短范围文本，如 '38-208' 或 '38 208' 返回 (38, 208)"""
    if not truncation_text:
        return None
    
    try:
        # 使用正则表达式提取数字
        numbers = re.findall(r'\d+', str(truncation_text))
        if len(numbers) >= 2:
            start = int(numbers[0])
            end = int(numbers[1])
            # 确保起始位置小于结束位置
            if start < end:
                return (start, end)
        elif len(numbers) == 1:
            # 如果只有一个数字，认为是起始位置，结束位置为序列末尾
            return (int(numbers[0]), None)
        return None
    except Exception as e:
        st.error(f"解析截短范围时出错: {e}")
        return None

def truncate_sequence(sequence, truncation_range):
    """根据截短范围截取序列"""
    if not sequence or not truncation_range:
        return sequence
    
    start, end = truncation_range
    
    # 调整索引（序列从1开始，Python从0开始）
    start_idx = start - 1
    
    # 如果结束位置为None，则截取到序列末尾
    if end is None:
        end_idx = len(sequence)
    else:
        end_idx = end
    
    # 检查范围是否有效
    if start_idx < 0 or end_idx > len(sequence) or start_idx >= end_idx:
        st.error(f"截短范围无效: {start}-{end}，序列长度: {len(sequence)}")
        return sequence
    
    truncated_seq = sequence[start_idx:end_idx]
    st.info(f"序列截短: 从位置 {start} 到 {end if end else '末尾'}，截短后长度: {len(truncated_seq)}")
    return truncated_seq

def add_tag_to_sequence(sequence, tag):
    """给序列添加标签"""
    if not tag or tag == "无标签":
        return sequence
    
    tag_sequence = TAG_SEQUENCES.get(tag)
    if tag_sequence:
        tagged_sequence = tag_sequence + sequence
        st.info(f"添加 {tag} 标签，标签长度: {len(tag_sequence)}，总长度: {len(tagged_sequence)}")
        return tagged_sequence
    else:
        st.warning(f"未知标签: {tag}")
        return sequence

def calculate_protein_properties(sequence):
    """计算蛋白的各种物理化学属性"""
    if not sequence:
        return None, None, None, None, None
    
    try:
        analyzed_seq = ProteinAnalysis(sequence)
        
        # 分子量 (转换为kD)
        molecular_weight = analyzed_seq.molecular_weight() / 1000.0
        
        # 等电点
        isoelectric_point = analyzed_seq.isoelectric_point()
        
        # 消光系数 (选择半胱氨酸形成二硫键的情况)
        extinction_coeff = analyzed_seq.molar_extinction_coefficient()[0]
        
        # 不稳定指数
        instability_index = analyzed_seq.instability_index()
        
        # GRAVY (疏水性)
        gravy = analyzed_seq.gravy()
        
        return molecular_weight, isoelectric_point, extinction_coeff, instability_index, gravy
    
    except Exception as e:
        st.error(f"计算错误: {e}")
        return None, None, None, None, None

# 主界面
def main():
    # 创建两列布局
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.header("输入参数")
        
        # 蛋白质标识符输入
        protein_input = st.text_input(
            "蛋白质名称或UniProt ID",
            placeholder="例如：P01308 或 Insulin",
            help="输入UniProt ID（如P01308）或蛋白质名称（如Insulin）"
        )
        
        # 截短信息输入
        truncation_input = st.text_input(
            "截短范围（可选）",
            placeholder="例如：38-208 或 38 208",
            help="输入截短范围，格式：起始位置-结束位置"
        )
        
        # 标签选择
        tag_selection = st.selectbox(
            "选择标签（可选）",
            options=["无标签", "10his", "6his", "GST", "SUMO"],
            help="选择要添加到蛋白质N端的标签"
        )
    
    with col2:
        st.header("计算结果")
        
        if protein_input:
            with st.spinner("正在获取蛋白质序列并计算性质..."):
                # 获取蛋白质序列
                sequence, uniprot_id = get_protein_sequence(protein_input.strip())
                
                if sequence:
                    # 显示基本信息
                    st.success(f"成功获取序列！UniProt ID: {uniprot_id}")
                    st.info(f"完整序列长度: {len(sequence)} 个氨基酸")
                    
                    # 处理截短
                    truncation_range = parse_truncation_range(truncation_input)
                    processed_sequence = sequence
                    
                    if truncation_range:
                        processed_sequence = truncate_sequence(processed_sequence, truncation_range)
                    
                    # 处理标签
                    if tag_selection != "无标签":
                        processed_sequence = add_tag_to_sequence(processed_sequence, tag_selection)
                    
                    # 计算性质
                    mw, pi, ext_coeff, instab, gravy = calculate_protein_properties(processed_sequence)
                    
                    # 显示结果
                    if mw is not None:
                        # 创建结果表格
                        results_df = pd.DataFrame({
                            '性质': ['分子量 (kD)', '等电点', '消光系数', '不稳定指数', 'GRAVY'],
                            '值': [f"{mw:.2f}", f"{pi:.2f}", f"{ext_coeff:.0f}", f"{instab:.2f}", f"{gravy:.3f}"]
                        })
                        
                        # 隐藏索引并显示表格
                        st.table(results_df.set_index('性质'))
                        
                        # 显示序列信息摘要
                        st.subheader("序列信息摘要")
                        summary_data = {
                            "处理类型": [],
                            "序列长度": []
                        }
                        
                        summary_data["处理类型"].append("原始序列")
                        summary_data["序列长度"].append(len(sequence))
                        
                        if truncation_range:
                            summary_data["处理类型"].append("截短后序列")
                            summary_data["序列长度"].append(len(truncate_sequence(sequence, truncation_range)))
                        
                        if tag_selection != "无标签":
                            summary_data["处理类型"].append("添加标签后")
                            summary_data["序列长度"].append(len(processed_sequence))
                        
                        summary_df = pd.DataFrame(summary_data)
                        st.table(summary_df)
                        
                    else:
                        st.error("无法计算蛋白质性质")
                else:
                    st.error("无法获取蛋白质序列，请检查输入是否正确")
        
        else:
            st.info("请在左侧输入蛋白质名称或UniProt ID")

# 运行主程序
main()

# 页脚信息
st.sidebar.markdown("---")
st.sidebar.info("""
**使用说明：**
- 输入UniProt ID（如P01308）或蛋白质名称（如Insulin）
- 可选：指定截短范围和添加标签
- 系统将自动获取序列并计算各种物理化学性质
""")