# Run on your NixOS computer

The app runs in a Debian 13 container managed by Podman. Its source comes from
your existing folder, `/home/peter/amir_design`. Runtime packages, the pinned
shRNAI/SeedMatchR tools, and reference data are downloaded during initial setup.
Use a 64-bit Intel/AMD computer and allow roughly 10 GB of free disk space. Larger
screens can need more RAM and disk space for their detailed results.

## Update your existing installation

1. Stop the app by pressing **Ctrl+C** in the terminal where it is running.
2. Open [the GitHub repository](https://github.com/petikalman003/amir_design), make
   sure the branch is **main**, and click **Code → Download ZIP**.
3. Extract the ZIP. Open its outer folder, usually named `amir_design-main`.
   Copy **the contents of that folder**, including `app.py`, `amir_design`,
   `scripts`, `pyproject.toml`, and `uv.lock`, into `/home/peter/amir_design`.
   Merge folders and replace matching files. Do not place the entire outer folder
   inside your existing project. Show hidden files in your file manager if needed
   so that the updated `.python-version` and `.gitignore` are included.

   Preserve your existing `data`, `tools`, `.venv`, `.runtime`, and `.git`
   directories. The source ZIP does not contain these installed dependencies,
   databases, or repository history; updating the source does not require deleting
   them or downloading the app again on each launch.
4. Restart the existing container, without `sudo`:

   ```bash
   podman start -ai amir-design-local
   ```

   The container runs setup before starting the app. It uses the updated frozen
   Python lockfile to install the new XlsxWriter dependency automatically, and
   reuses your existing reference databases and runtime downloads. Matching CDS/UTR
   annotations are added if missing. Leave the terminal open until Streamlit starts.
   This update also includes the previous BLAS/LAPACK and Streamlit launcher repairs.
5. Open `http://localhost:8501` in your browser and reload the page.

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
   expected, with a score near **0.4862** and coverage **1/1**. Progress shows the
   current stage, how many guides have been screened, and elapsed time.
5. Change **Sort results by**, then click **Download results and sequences**.
   This summary workbook includes every candidate, exact off-target counts,
   candidate origins, covered/missed transcript IDs, activity scores, and synthesis
   sequences. Sorting reuses the analysis and changes the summary download's order.
6. To inspect individual off-target positions, click **Prepare detailed workbook**.
   When it finishes, click **Download detailed workbook**. This optional file adds
   every position to the candidate/sequence report. It can be large and take longer
   to prepare. Its saved sorting order is shown in the app; changing the table sort
   keeps that file intact. Prepare it again if you want a different order in the
   detailed download.

For mouse Fpgs, select **Mouse** and enter `ENSMUSG00000009566`. The supplied M36
reference contains ten annotated transcripts. Unchecked all-isoform coverage
allows partial-coverage candidates; checking it requires coverage of all ten.
Choose a transcript ID to generate candidates from just that transcript while
still reporting coverage across all ten transcripts of its gene.

The default analysis screens all eligible candidates against the full mouse or
human reference and returns exact counts. It does not prepare millions of
off-target position rows before showing your results. An exact-seed prefilter
skips transcripts with no matching seed, then the unchanged SeedMatchR function
finds the positions; no candidates or matches are capped. Detailed positions are
kept on disk when requested, while R still needs enough memory for the full
transcript reference. The first 10,000 positions are shown in the preview; every
position is retained in the detailed workbook.

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
