# Minimal image for the rxnlab-pg-backup CronJob.
#
# Why a custom image instead of `alpine:3.18` + apk at runtime?
# OpenShift's restricted SCC runs containers as a non-root arbitrary UID, so
# `apk add` fails with "Permission denied" trying to write the apk database.
# Installing tools at build time (as root in the builder) bakes them into the
# image, and they're then usable at runtime regardless of the runtime UID.
#
# Build (re-)build:
#   oc start-build rxnlab-backup-image --from-dir=rahti/backup-image
#
# First-time build setup (already done — keep for reproducibility):
#   oc new-build --name=rxnlab-backup-image --binary --strategy=docker
#
# Sized for the CronJob — adds ~80 MB on top of base Alpine.

FROM alpine:3.18

RUN apk add --no-cache \
    aws-cli \
    postgresql-client \
    gzip \
    coreutils \
    bash

USER 1001
