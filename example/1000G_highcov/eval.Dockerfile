# Python + matplotlib for the fast-mode evaluation and its HTML report,
# published by CI as the :eval tag of the repository's image package so
# provisioning is an unprivileged `apptainer pull` on any node. Its own image
# on purpose: the analysis image is a scientific artifact, and plotting
# dependencies change on operational grounds.
FROM python:3.12-slim
RUN pip install --no-cache-dir matplotlib==3.9.2
ENTRYPOINT ["python3"]
