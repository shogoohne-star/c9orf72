# -*- coding: utf-8 -*-
"""テンプレートに style.css と SVG 資産を埋め込んで配布用 HTML を組み立てる。
プレースホルダ:  {{css}}  {{svg:NAME}}  {{html:NAME}}"""
import re, os, sys

def load_css():
    return open("style.css", encoding="utf-8").read()

def load_svg(name):
    s = open(f"assets/{name}.svg", encoding="utf-8").read()
    s = re.sub(r"<\?xml[^>]*\?>\s*", "", s)                 # XML 宣言を除去
    s = re.sub(r"<!--.*?-->", "", s, flags=re.S)            # コメント除去
    if "rdkit" in s or "baseProfile" in s:                  # RDKit 出力の整形
        s = re.sub(r"\s(width|height)='[^']*'", "", s, count=2)
        s = s.replace("#19191E", "currentColor")
        s = s.replace("<svg version='1.1'", "<svg class='molsvg' preserveAspectRatio='xMidYMid meet' version='1.1'")
    return s.strip()

def load_html(name):
    return open(f"assets/{name}.html", encoding="utf-8").read().strip()

def build(tpl, out):
    s = open(tpl, encoding="utf-8").read()
    s = s.replace("{{css}}", load_css())
    s = re.sub(r"\{\{svg:([A-Za-z0-9_]+)\}\}", lambda m: load_svg(m.group(1)), s)
    s = re.sub(r"\{\{html:([A-Za-z0-9_]+)\}\}", lambda m: load_html(m.group(1)), s)
    left = re.findall(r"\{\{[^}]+\}\}", s)
    if left:
        print("  !! 未解決のプレースホルダ:", set(left)); sys.exit(1)
    open(out, "w", encoding="utf-8").write(s)
    print(f"{out}: {len(s):,} bytes")

for t, o in (("summary.tpl.html", "summary.html"), ("detail.tpl.html", "detail.html")):
    if os.path.exists(t):
        build(t, o)
