# 导入 OpenAI 客户端
from openai import OpenAI
 
# 创建客户端实例
client = OpenAI(
    base_url="https://api.acemap.cn/v1",  # IIoT 提供的 API 基础 URL
    api_key="sk-z50g76lgvZxb5Bpdll1t0X9JAlYJoJuwNGAr2Kpc6Qt8W8Ur",  # 您的 API 密钥
    timeout=30.0
)

try:
    response = client.chat.completions.create(
        model="deepseek-r1:14b",
        messages=[
            {"role": "system", "content": "你是一个有用的助手。"},
            {"role": "user", "content": "解释量子计算的基本原理。"}
        ]
    )
    print(response.choices[0].message.content)
except Exception as e:
    print(f"完整错误信息: {e}")
    print(f"错误类型: {type(e)}")

try:
    models = client.models.list()
    print("可用模型:", [model.id for model in models.data])
except Exception as e:
    print(f"列表模型错误: {e}")