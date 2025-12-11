#!/usr/bin/env python3
"""
简单的HTTP服务器，用于提供前端界面
使用方法：python frontend_server.py
然后在浏览器访问 http://localhost:8090
"""

import http.server
import socketserver
import os
from pathlib import Path

# 配置
PORT = 8090
DIRECTORY = Path(__file__).parent

class MyHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(DIRECTORY), **kwargs)
    
    def do_GET(self):
        # 如果访问根路径，重定向到 frontend.html
        if self.path == '/' or self.path == '':
            self.path = '/frontend.html'
        return super().do_GET()
    
    def end_headers(self):
        # 添加CORS头，允许跨域访问
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        self.send_header('Cache-Control', 'no-store, no-cache, must-revalidate')
        super().end_headers()

def main():
    # 检查frontend.html是否存在
    frontend_file = DIRECTORY / 'frontend.html'
    if not frontend_file.exists():
        print(f"错误: 找不到 {frontend_file}")
        print("请确保 frontend.html 在当前目录下")
        return
    
    with socketserver.TCPServer(("", PORT), MyHTTPRequestHandler) as httpd:
        print("=" * 60)
        print(f"🌊 海洋科学智能问答系统 - 前端服务器")
        print("=" * 60)
        print(f"\n✅ 服务器已启动")
        print(f"📍 访问地址: http://localhost:{PORT}/frontend.html")
        print(f"📁 服务目录: {DIRECTORY}")
        print(f"\n⚠️  请确保后端服务器已在 http://localhost:1030 运行")
        print(f"   (运行命令: python server_v3.py)\n")
        print("按 Ctrl+C 停止服务器\n")
        
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\n\n👋 服务器已停止")

if __name__ == "__main__":
    main()