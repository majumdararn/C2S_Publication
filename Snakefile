rule fig_ball:
    input:
        "src/data/dataset.tar.xz"
    output:
        "src/figures/Iq_ball.pdf"
    conda:
        "environment.yml"
    script:
        "src/scripts/ball_plot_paper.py"