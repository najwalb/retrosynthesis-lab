#!/usr/bin/env bash
# One-shot bootstrap for the RxnLab feedback persistence layer on Rahti.
#
#  - Creates the rxnlab-pg Secret (database-user/-password/-name) with a random password.
#  - Creates the rxnlab-db Secret (DATABASE_URL) that the app consumes.
#  - Applies rahti/postgres.yaml.
#  - Applies rahti/postgres-backup.yaml (suspended; enable separately — see README).
#
# Idempotent: existing Secrets are left untouched; existing Deployments/Services
# are reconciled.
#
# Usage:    bash rahti/setup.sh
#           (must be `oc login`-ed and in the right project first)

set -euo pipefail

cd "$(dirname "$0")"

NAMESPACE=$(oc project -q)
echo "Project: $NAMESPACE"
echo

# ── Secrets ───────────────────────────────────────────────────────────────────
if oc get secret rxnlab-pg >/dev/null 2>&1; then
    echo "secret/rxnlab-pg already exists — leaving it alone"
else
    PGPASS=$(openssl rand -base64 24 | tr -d '/+=' | head -c 24)
    oc create secret generic rxnlab-pg \
        --from-literal=database-user=rxnlab \
        --from-literal=database-password="$PGPASS" \
        --from-literal=database-name=rxnlab
    echo "secret/rxnlab-pg created"
fi

if oc get secret rxnlab-db >/dev/null 2>&1; then
    echo "secret/rxnlab-db already exists — leaving it alone"
else
    PGPASS=$(oc get secret rxnlab-pg -o jsonpath='{.data.database-password}' | base64 -d)
    oc create secret generic rxnlab-db \
        --from-literal=DATABASE_URL="postgresql+psycopg2://rxnlab:${PGPASS}@rxnlab-pg:5432/rxnlab"
    echo "secret/rxnlab-db created"
fi

# ── Postgres ──────────────────────────────────────────────────────────────────
oc apply -f postgres.yaml
oc rollout status deploy/rxnlab-pg --timeout=3m

# ── Backup CronJob (suspended) ────────────────────────────────────────────────
oc apply -f postgres-backup.yaml
echo
echo "rxnlab-pg-backup CronJob is SUSPENDED. To enable, see rahti/README.md."

echo
echo "── Next step ─────────────────────────────────────────────────────────────"
echo "Wire the app to the DB:"
echo "  oc set env deployment/<app-name> --from=secret/rxnlab-db"
