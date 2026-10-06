"""v8: record the software environment at build time (release_v8/environment.txt)."""
import subprocess, importlib, platform, datetime, os
os.makedirs('release_v8', exist_ok=True)
L = [f'# Environment of the v8 build, captured {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}', ' '.join(platform.uname()[:3] + (platform.machine(),)), 'Python ' + platform.python_version()]
for m in ['numpy', 'pandas', 'scipy', 'matplotlib', 'docx', 'openpyxl', 'PIL']:
    L.append(f'{m} {getattr(importlib.import_module(m), "__version__", "?")}')
run = lambda c: subprocess.run(c, shell=True, capture_output=True, text=True).stdout.strip()
L.append(run('plink1.9 --version | head -1'))
L.append(run('''Rscript -e 'for (p in c("MendelianRandomization","MRPRESSO","coloc","susieR")) cat(p, as.character(packageVersion(p)), "\\n"); cat(R.version.string, "\\n")' '''))
L.append(run('soffice --version | head -1'))
L.append('MVMR reference code: ref_MVMR/ from GitHub WSpiller/MVMR commit 8be0d94 (sourced)')
open('release_v8/environment.txt', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
