rule crw:
    input: "src/data/crw.tar.gz"
    output:
        fig="src/tex/figures/crw.pdf",
        qp0="src/tex/output/qp0.txt",
        qp1="src/tex/output/qp1.txt",
        qp2="src/tex/output/qp2.txt",
        qres0="src/tex/output/qres0.txt",
        qres1="src/tex/output/qres1.txt"
    script: "src/scripts/sas_val_crw.py"