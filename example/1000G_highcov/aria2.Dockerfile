# aria2c for the 1000G download manager, published by CI as the :aria2 tag of
# the repository's image package so provisioning is an unprivileged
# `apptainer pull` on any node - compute nodes included, where fakeroot
# builds of aria2.def have failed. Deliberately its own image: the analysis
# image is a scientific artifact, a downloader is operational plumbing.
FROM alpine:3.20
RUN apk add --no-cache aria2 ca-certificates
ENTRYPOINT ["aria2c"]
