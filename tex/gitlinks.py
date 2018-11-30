import subprocess
hash = subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode("utf-8")[:-1]
with open("gitlinks.tex", "w") as f:
    print(r"\newcommand{\pycodelink}[1]{\href{https://github.com/rodluger/limbdark/blob/%s/tex/figures/python/#1.py}{\codeicon}\,\,}" % hash, file=f)
    print(r"\newcommand{\jlcodelink}[1]{\href{https://github.com/rodluger/limbdark/blob/%s/tex/figures/julia/#1.jl}{\codeicon}\,\,}" % hash, file=f)
    print(r"\newcommand{\animlink}[1]{\href{https://github.com/rodluger/limbdark/blob/%s/tex/figures/#1.gif}{\animicon}\,\,}" % hash, file=f)
    print(r"\newcommand{\prooflink}[1]{\href{https://github.com/rodluger/limbdark/blob/%s/proofs/#1.ipynb}{\raisebox{-0.1em}{\prooficon}}}" % hash, file=f)
