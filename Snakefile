rule ball_fig_gen:
    input: 
        'src/scripts/ball_plot_paper.py',
        'src/data/dataset.tar.xz'
    output:
        'src/tex/figures/Iq_ball.pdf'
    conda:
        'environment.yml'
rule bib_fig_gen:
    input: 
        'src/scripts/bib_ecc_plot_paper.py',
        'src/data/dataset.tar.xz'
    output:
        'src/tex/figures/Iq_bib_ecc.pdf'
    conda:
        'environment.yml'