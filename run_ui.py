"""
ProtFlow Studio Launcher:
Starts the ComfyUI-style node editor web interface on localhost.
"""

import sys
import webbrowser
from protflow.ui.server import start_ui_server

if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8080
    print(f"[*] ProtFlow ComfyUI Studio başlatılıyor...")
    print(f"[*] Tarayıcınızda açın: http://127.0.0.1:{port}")
    start_ui_server(port=port)
