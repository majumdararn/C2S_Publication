rule ball_plot:
    input:
        "src/data/dataset.tar.xz"
    output:
        "src/figures/Iq_ball.pdf"
    conda:
        "environment.yml"
    script:
        "src/scripts/ball_plot_paper.py"

rule box_plot:
    input:
        "src/data/dataset.tar.xz"
    output:
        "src/figures/Iq_box_orien.pdf"
    conda:
        "environment.yml"
    script:
        "src/scripts/box_plot_paper.py"