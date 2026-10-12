#!/usr/bin/env bash
# Rootless Debian 13 R/Biostrings and BLAST+ runtimes for the three-stage pipeline.
# APT checks Debian archive signatures and package hashes; no system files are changed.
set -euo pipefail

runtime_prefix="${AMIR_RUNTIME_DIR:-/workspace/.runtime}"
runtime_prefix="$(mkdir -p "$runtime_prefix" && cd "$runtime_prefix" && pwd)"
r_tree="$runtime_prefix/r"
apt_tree="$runtime_prefix/apt-r"
runtime_user="$(id -un)"

if [[ -x "$runtime_prefix/bin/Rscript" && -x "$runtime_prefix/bin/blastn" && -x "$runtime_prefix/bin/makeblastdb" ]] && "$runtime_prefix/bin/Rscript" --vanilla -e \
  'stopifnot(getRversion() >= "4.3.0", requireNamespace("Biostrings", quietly=TRUE))' >/dev/null 2>&1 \
  && "$runtime_prefix/bin/blastn" -version >/dev/null 2>&1 \
  && "$runtime_prefix/bin/makeblastdb" -version >/dev/null 2>&1; then
  printf 'R runtime ready: %s\n' "$runtime_prefix/bin/Rscript"
  exit 0
fi

if ! rg -q '^VERSION_CODENAME=trixie$' /etc/os-release; then
  printf 'This rootless binary installer requires Debian 13; install R >=4.3 and Biostrings with your platform package manager.\n' >&2
  exit 1
fi
mkdir -p "$apt_tree/lists/partial" "$apt_tree/cache/archives/partial" \
  "$apt_tree/parts" "$apt_tree/sourceparts" "$r_tree" "$runtime_prefix/bin"
cp /var/lib/dpkg/status "$apt_tree/status"
cat > "$apt_tree/sources.list" <<'EOF'
deb [signed-by=/usr/share/keyrings/debian-archive-keyring.gpg] https://deb.debian.org/debian trixie main
EOF
cat > "$apt_tree/apt.conf" <<EOF
Dir::Etc::Parts "$apt_tree/parts";
Dir::Etc::Main "$apt_tree/apt.conf";
Dir::Etc::sourcelist "$apt_tree/sources.list";
Dir::Etc::sourceparts "$apt_tree/sourceparts";
Dir::State "$apt_tree";
Dir::State::status "$apt_tree/status";
Dir::State::lists "$apt_tree/lists";
Dir::Cache "$apt_tree/cache";
APT::Sandbox::User "$runtime_user";
EOF
APT_CONFIG="$apt_tree/apt.conf" /usr/bin/apt-get -o APT::Update::Error-Mode=any update
APT_CONFIG="$apt_tree/apt.conf" /usr/bin/apt-get --download-only --assume-yes \
  --no-install-recommends install r-base-core r-bioc-biostrings ncbi-blast+
for archive in "$apt_tree"/cache/archives/*.deb; do
  dpkg-deb --extract "$archive" "$r_tree"
done

# Debian normally configures these files during privileged package installation.
# Relocate only our extracted copy to make that runtime usable in the workspace.
python3 - "$r_tree" "$runtime_prefix" <<'PY'
from pathlib import Path
import shlex
import sys

tree, prefix = map(Path, sys.argv[1:])
config = tree / 'etc/R'
(config / 'Renviron').write_bytes((tree / 'usr/lib/R/etc/Renviron.ucf').read_bytes())
for path in (tree / 'usr/lib/R/etc').iterdir():
    if path.is_symlink() and str(path.readlink()).startswith('/etc/R/'):
        destination = tree / str(path.readlink()).lstrip('/')
        path.unlink()
        path.symlink_to(destination)
r_launcher = tree / 'usr/lib/R/bin/R'
text = r_launcher.read_text()
for key, value in {
    'R_HOME_DIR': tree / 'usr/lib/R',
    'R_SHARE_DIR': tree / 'usr/share/R/share',
    'R_INCLUDE_DIR': tree / 'usr/share/R/include',
    'R_DOC_DIR': tree / 'usr/share/R/doc',
}.items():
    lines = text.splitlines(keepends=True)
    text = ''.join(f'{key}={shlex.quote(str(value))}\n' if line.startswith(key + '=') else line for line in lines)
r_launcher.write_text(text)
multiarch = __import__('subprocess').check_output(['dpkg-architecture', '-qDEB_HOST_MULTIARCH'], text=True).strip()
# These packages normally expose their libraries through update-alternatives.
# Extracting .deb archives does not run that configuration step, so include the
# actual BLAS/LAPACK package directories in our relocated R library search path.
r_library_dirs = [
    tree / 'usr/lib/R/lib',
    tree / 'usr/lib' / multiarch,
    tree / 'usr/lib' / multiarch / 'blas',
    tree / 'usr/lib' / multiarch / 'lapack',
]
launcher = prefix / 'bin/Rscript'
launcher.write_text(
    '#!/usr/bin/env bash\nset -euo pipefail\n'
    f'export R_HOME={shlex.quote(str(tree / "usr/lib/R"))}\n'
    'export RHOME="$R_HOME"\n'
    f'export R_LIBS_SITE={shlex.quote(str(tree / "usr/lib/R/site-library"))}\n'
    f'export LD_LIBRARY_PATH={shlex.quote(":".join(map(str, r_library_dirs)))}:"${{LD_LIBRARY_PATH:-}}"\n'
    f'exec {shlex.quote(str(tree / "usr/bin/Rscript"))} "$@"\n'
)
launcher.chmod(0o755)
for executable in ('blastn', 'makeblastdb'):
    launcher = prefix / 'bin' / executable
    launcher.write_text(
        '#!/usr/bin/env bash\nset -euo pipefail\n'
        f'export LD_LIBRARY_PATH={shlex.quote(str(tree / "usr/lib/ncbi-blast+") + ":" + str(tree / "usr/lib" / multiarch))}:"${{LD_LIBRARY_PATH:-}}"\n'
        f'export BLASTMAT={shlex.quote(str(tree / "usr/share/ncbi/data"))}\n'
        f'exec {shlex.quote(str(tree / "usr/bin" / executable))} "$@"\n'
    )
    launcher.chmod(0o755)
PY
"$runtime_prefix/bin/Rscript" --vanilla -e \
  'stopifnot(getRversion() >= "4.3.0"); library(Biostrings); cat("R", as.character(getRversion()), "Biostrings", as.character(packageVersion("Biostrings")), "\n")'
"$runtime_prefix/bin/blastn" -version
"$runtime_prefix/bin/makeblastdb" -version
for archive in "$apt_tree"/cache/archives/*.deb; do
  dpkg-deb --show "$archive"
done > "$runtime_prefix/native-packages.tsv"
printf 'R runtime ready: %s\n' "$runtime_prefix/bin/Rscript"
