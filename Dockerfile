FROM ubuntu:21.10
ENV DEBIAN_FRONTEND noninteractive
ENV package_dir /usr/local/lib/python3.9/dist-packages/fretlabel/
RUN apt-get update && \
    apt-get install -y --no-install-recommends pymol pip midori && \
    python3 -m pip install -U pip && \
    rm -rf /var/lib/apt/lists/* && \
    pip install fretlabel && \
    cp "$package_dir"/fretlabel_gui.py /usr/lib/python3/dist-packages/pmg_tk/startup && \
    echo {\"browser\": null, \"local_docs\": null} > "$package_dir"/.fretlabel_settings.json
CMD pymol
