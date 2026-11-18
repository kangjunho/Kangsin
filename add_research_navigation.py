from pathlib import Path
import re
import shutil

root = Path(__file__).resolve().parent
pages = ('index.html','about.html','members.html','publications.html','e-learning.html','bioinformatics-analysis.html','contact.html')
pattern = re.compile(r'(<li\b[^>]*>\s*<a\b[^>]*href=["\']publications\.html["\'][^>]*>\s*Publications\s*</a>\s*</li>)', re.I)
link = '<li class="nav-item"><a class="nav-link" href="research-code.html">Research &amp; Code</a></li>'
for name in pages:
    p = root / name
    if not p.is_file():
        print('Missing:', name); continue
    old = p.read_text(encoding='utf-8')
    if 'href="research-code.html"' in old:
        print('Already linked:', name); continue
    new, n = pattern.subn(lambda m: m.group(0)+'\n'+link, old, count=1)
    if n != 1:
        print('No matching menu; inspect manually:', name); continue
    new = re.sub(r'\bnavbar-expand-(?:sm|md|lg)\b','navbar-expand-xl',new,count=1)
    shutil.copy2(p, p.with_suffix('.html.bak'))
    p.write_text(new, encoding='utf-8')
    print('Updated:', name)
