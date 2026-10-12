# Run on your NixOS computer

The app runs in a Debian 13 container managed by Podman. Its source comes from
your existing folder, `/home/peter/amir_design`. Runtime packages, the pinned
shRNAI/SeedMatchR tools, and reference data are downloaded during initial setup.
Use a 64-bit Intel/AMD computer and allow roughly 10 GB of free disk space. Larger
screens can need more RAM and disk space for their detailed results.

## Update your existing installation

1. Stop the app by pressing **Ctrl+C** in the terminal where it is running.
2. Download `amir-design-update.zip` from the chat and extract its contents into
   `/home/peter/amir_design`, replacing matching project files. The archive contains
   source, documentation, and tests; it does not contain your `.venv`, databases,
   upstream tools, native runtime, or `.git` directory.
3. Restart the existing container, without `sudo`:

   ```bash
   podman start -ai amir-design-local
   ```

   Setup reuses your current transcript references and adds the matching CDS/UTR
   annotations. Leave the terminal open until Streamlit starts. Existing downloads
   are reused. The update includes the BLAS/LAPACK and Streamlit launcher repairs.
4. Open `http://localhost:8501` in your browser and reload the page.

If you prefer to extract the update from Terminal, and the archive is in Downloads:

```bash
nix-shell -p unzip --run 'unzip -o /home/peter/Downloads/amir-design-update.zip -d /home/peter/amir_design'
```

For subsequent launches, use `podman start -ai amir-design-local` and open the same
browser address. **Ctrl+C** stops the app; the named container preserves its
installed dependencies and reference files.

## First installation

Check Podman without `sudo`:

```bash
podman info --format '{{.Host.Arch}}'
```

It should print `amd64`. Then copy this whole block into Terminal once:

```bash
podman run --name amir-design-local -it \
  -p 127.0.0.1:8501:8501 \
  -v /home/peter/amir_design:/workspace/amir_design \
  docker.io/library/debian:trixie bash -lc '
    set -euo pipefail
    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends \
      ca-certificates curl git python3 ripgrep dpkg-dev libstdc++6
    if [ ! -x /root/.local/bin/uv ]; then
      curl -LsSf https://astral.sh/uv/install.sh | sh
    fi
    export PATH="/root/.local/bin:$PATH"
    cd /workspace/amir_design
    bash scripts/setup.sh --references
    exec bash scripts/start.sh
  '
```

When Streamlit starts, open `http://localhost:8501` in your browser. Keep Terminal
open while using the app. The published port is available only on your own
computer.

## Test the new form

1. Select **Human**.
2. Enter `ENSG00000185962` (LCE3A).
3. Select **CDS only**, leave **Require all isoforms** unchecked, and retain the
   **Uploaded Design.py rules** profile.
4. Click **Design and screen**. With the supplied v47 reference, one candidate is
   expected, with a score near **0.4862** and coverage **1/1**.
5. Change **Sort results by** and download the workbook. The analysis is reused
   when sorting changes. The workbook includes candidate origins and covered/missed
   transcript IDs as well as activity and off-target metrics.

For mouse Fpgs, select **Mouse** and enter `ENSMUSG00000009566`. The supplied M36
reference contains ten annotated transcripts. Unchecked all-isoform coverage
allows partial-coverage candidates; checking it requires coverage of all ten.
Choose a transcript ID to generate candidates from just that transcript while
still reporting coverage across all ten transcripts of its gene.

The reference's annotated transcript count is not an expression measurement.
Scores are uncalibrated activity predictions and seed hits are potential binding
sites. These examples verify operation rather than recommend an experimental
guide.

## Run automated checks

While the app's container is running, open a second terminal and paste:

```bash
podman exec -it amir-design-local bash -lc '
  cd /workspace/amir_design
  export PATH="/workspace/.runtime/bin:$PATH"
  export TF_CPP_MIN_LOG_LEVEL=2
  export TF_NUM_INTRAOP_THREADS=2
  export TF_NUM_INTEROP_THREADS=1
  .venv/bin/python -m pytest -q
'
```

If setup or a check fails, copy its error text before making further changes.
