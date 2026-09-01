# -*- coding: utf-8 -*-
"""docx-js の出力を仕上げる：Normal を既定段落スタイルに、
和文の行頭禁則と英数字まわりの自動空白を有効にする。"""
import re, shutil, sys, zipfile

SETTINGS_EXTRA = (
    '<w:compat><w:compatSetting w:name="compatibilityMode" '
    'w:uri="http://schemas.microsoft.com/office/word" w:val="15"/></w:compat>'
)

def fix(path):
    src = zipfile.ZipFile(path)
    items = [(i, src.read(i.filename)) for i in src.infolist()]
    src.close()
    out = []
    for info, data in items:
        if info.filename == 'word/styles.xml':
            s = data.decode('utf-8')
            s = s.replace('<w:style w:type="paragraph" w:styleId="Normal">',
                          '<w:style w:type="paragraph" w:default="1" w:styleId="Normal">')
            data = s.encode('utf-8')
        elif info.filename == 'word/settings.xml':
            s = data.decode('utf-8')
            if '<w:compat>' not in s:
                s = s.replace('</w:settings>', SETTINGS_EXTRA + '</w:settings>')
            data = s.encode('utf-8')
        out.append((info, data))
    tmp = path + '.tmp'
    with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as z:
        for info, data in out:
            z.writestr(info, data)
    shutil.move(tmp, path)
    print(f'fixed {path}')

for p in sys.argv[1:]:
    fix(p)
