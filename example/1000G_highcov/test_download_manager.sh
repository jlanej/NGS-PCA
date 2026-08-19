#!/usr/bin/env bash
# =============================================================================
# test_download_manager.sh — hermetic end-to-end test of the download manager
# =============================================================================
#
# Exercises 01_download_and_mosdepth.sh and 01b_mosdepth_sample.sh with every
# external dependency shimmed: sbatch runs 01b inline, squeue reports nothing,
# apptainer fakes mosdepth, wget serves files from a local store, and curl and
# aria2c fail so nothing touches the network. Five samples, one of them
# corrupt, cover the full cycle: download, verify, dispatch, mosdepth, clean
# up, fail the corrupt one after a re-download, then heal it on a resweep.
#
# Usage: test_download_manager.sh <scratch_dir>
# =============================================================================

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
T="${1:?usage: test_download_manager.sh <scratch_dir>}"
rm -rf "${T}"
mkdir -p "${T}/bin" "${T}/remote" "${T}/work/reference"

# ── Shims ────────────────────────────────────────────────────────────────────

cat > "${T}/bin/md5sum" <<'EOF'
#!/usr/bin/env bash
python3 -c 'import hashlib,sys; print(hashlib.md5(open(sys.argv[1],"rb").read()).hexdigest()+"  "+sys.argv[1])' "$1"
EOF

cat > "${T}/bin/wget" <<EOF
#!/usr/bin/env bash
# serves \$FAKE_REMOTE/<basename of url>; understands -O <dest>
dest=""; url=""
while [[ \$# -gt 0 ]]; do
  case "\$1" in
    -O) dest="\$2"; shift 2 ;;
    -*) shift ;;
    *) url="\$1"; shift ;;
  esac
done
src="${T}/remote/\$(basename "\${url}")"
[[ -f "\${src}" ]] || exit 1
cp "\${src}" "\${dest}"
EOF

cat > "${T}/bin/curl" <<'EOF'
#!/usr/bin/env bash
exit 1
EOF
cat > "${T}/bin/aria2c" <<'EOF'
#!/usr/bin/env bash
exit 1
EOF

cat > "${T}/bin/squeue" <<'EOF'
#!/usr/bin/env bash
exit 0
EOF

cat > "${T}/bin/sbatch" <<'EOF'
#!/usr/bin/env bash
# runs the submitted script inline and prints a fake job id, like --parsable
args=()
for a in "$@"; do
  [[ "${a}" == -* ]] && continue
  args+=("${a}")
done
bash "${args[@]}" > /dev/null 2>&1 || true
echo $(( RANDOM + 1000 ))
EOF

cat > "${T}/bin/apptainer" <<'EOF'
#!/usr/bin/env bash
# fakes "apptainer exec ... mosdepth ...": --version prints one, a real run
# writes <prefix>.regions.bed.gz under the directory bound to /mosdepth
if [[ "$*" == *"--version"* ]]; then
  echo "mosdepth 9.9.9-test"
  exit 0
fi
mosdepth_host=""
prefix=""
prev=""
for a in "$@"; do
  if [[ "${prev}" == "--bind" && "${a}" == *:/mosdepth ]]; then
    mosdepth_host="${a%:/mosdepth}"
  fi
  [[ "${a}" == /mosdepth/* ]] && prefix="${a#/mosdepth/}"
  prev="${a}"
done
[[ -n "${mosdepth_host}" && -n "${prefix}" ]] || exit 1
echo "fake-coverage" | gzip > "${mosdepth_host}/${prefix}.regions.bed.gz"
EOF

chmod +x "${T}/bin/"*
export PATH="${T}/bin:${PATH}"

# ── Fixture: five samples, S5 corrupt (manifest md5 of different content) ────

python3 - "${T}" <<'EOF'
import hashlib, os, sys
T = sys.argv[1]
rows = []
for i in range(1, 6):
    sample = f"S{i}"
    content = (f"cram-bytes-{i}-" * 200).encode()
    with open(f"{T}/remote/{sample}.final.cram", "wb") as out:
        out.write(content)
    with open(f"{T}/remote/{sample}.final.cram.crai", "wb") as out:
        out.write(b"crai" + str(i).encode())
    md5 = hashlib.md5(content).hexdigest()
    if sample == "S5":
        md5 = hashlib.md5(b"what the manifest thinks S5 is").hexdigest()
    base = f"ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR{i}/{sample}.final.cram"
    rows.append(f"{sample}\t{base}\t{base}.crai\t{md5}\textra")
with open(f"{T}/work/manifest.tsv", "w") as out:
    out.write("SAMPLE\tCRAM\tCRAI\tMD5\tNOTES\n")
    out.write("\n".join(rows) + "\n")
print("fixture written")
EOF

touch "${T}/work/ngs-pca.sif" "${T}/work/reference/ref.fa"

run_manager() {
  WORK_DIR="${T}/work" \
  SIF_IMAGE="${T}/work/ngs-pca.sif" \
  REF_DIR="${T}/work/reference" \
  REF_FASTA="${T}/work/reference/ref.fa" \
  MANIFEST="${T}/work/manifest.tsv" \
  MIN_MANIFEST_SAMPLES=1 \
  EXPECTED_MANIFEST_SAMPLES=5 \
  USE_ASPERA=0 \
  COMPARE_FAST_MODE=1 \
  DOWNLOADER_LOCAL=1 \
  DOWNLOAD_SLOTS=3 \
  MAX_LOCAL_CRAMS="${MAX_LOCAL_CRAMS:-10}" \
  FAKE_REMOTE="${T}/remote" \
  bash "${HERE}/01_download_and_mosdepth.sh"
}

run_generator() {
  WORK_DIR="${T}/work" \
  MANIFEST="${T}/work/manifest.tsv" \
  COMPARE_FAST_MODE=1 \
  bash "${HERE}/globus_batch_from_manifest.sh"
}

check() { # label condition-description: pass a test command as the rest
  local label="$1"; shift
  if "$@"; then echo "PASS: ${label}"; else echo "FAIL: ${label}"; exit 1; fi
}

echo "=== Sweep 1: four good samples, one corrupt ==="
rc=0
run_manager > "${T}/sweep1.log" 2>&1 || rc=$?
check "manager exits non-zero when a sample fails" test "${rc}" -ne 0
grep -q "4 dispatched to mosdepth, 1 failed" "${T}/sweep1.log" \
  || { echo "FAIL: summary line"; cat "${T}/sweep1.log"; exit 1; }
echo "PASS: summary counts 4 ok / 1 failed"

for i in 1 2 3 4; do
  check "S${i} normal tree output" test -s "${T}/work/mosdepth_output/S${i}.by1000.regions.bed.gz"
  check "S${i} fast tree output" test -s "${T}/work/mosdepth_output_fast/S${i}.by1000.regions.bed.gz"
  check "S${i} timing rows" test -s "${T}/work/qc_output/mosdepth_timing/S${i}.normal.tsv" \
    -a -s "${T}/work/qc_output/mosdepth_timing/S${i}.fast.tsv"
done
check "no CRAMs left on disk" test -z "$(find "${T}/work/crams" -name '*.cram' 2>/dev/null)"
check "S5 failed" grep -qx "fail" "${T}/work/download_state/results/S5"
check "S5 was re-downloaded before failing" \
  grep -q "re-downloading over HTTPS" "${T}/work/logs/download_S5.log"
check "distrust flag set by the mismatch" test -f "${T}/work/download_state/aspera_disabled"

echo "=== Sweep 2: S5 pre-staged (the Globus workflow), remote gone ==="
batch="$(run_generator 2> "${T}/generator.log")"
echo "${batch}" | grep -q "^\"/vol1/run/ERR5/S5.final.cram\" \"${T}/work/crams/S5.cram\"$" \
  || { echo "FAIL: generator batch line"; echo "${batch}"; exit 1; }
check "generator emits only the needed sample" test "$(echo "${batch}" | wc -l | tr -d '[:space:]')" = "2"
grep -q "Emitted 1 samples" "${T}/generator.log" || { echo "FAIL: generator summary"; exit 1; }
echo "PASS: generator emits S5's two files, quoted and renamed, skips the done four"

# stage what Globus would deliver, with the correct content, and remove the
# remote copies entirely: a download attempt for S5 would now fail loudly
python3 - "${T}" <<'EOF'
import os, sys
T = sys.argv[1]
with open(f"{T}/work/crams/S5.cram", "wb") as out:
    out.write(b"what the manifest thinks S5 is")
with open(f"{T}/work/crams/S5.cram.crai", "wb") as out:
    out.write(b"crai5")
os.remove(f"{T}/remote/S5.final.cram")
os.remove(f"{T}/remote/S5.final.cram.crai")
EOF
# MAX_LOCAL_CRAMS=1 with one staged CRAM on disk: the old backpressure would
# deadlock here; gating only real downloads must let the staged sample through
MAX_LOCAL_CRAMS=1 run_manager > "${T}/sweep2.log" 2>&1 \
  || { echo "FAIL: sweep 2 should succeed"; cat "${T}/sweep2.log" "${T}/work/logs/download_S5.log"; exit 1; }
grep -q "Already done:        4" "${T}/sweep2.log" || { echo "FAIL: sweep 2 skip count"; exit 1; }
echo "PASS: sweep 2 skips the four finished samples"
grep -q "CRAM already present" "${T}/work/logs/download_S5.log" \
  || { echo "FAIL: staged CRAM should skip the download"; exit 1; }
grep -q "CRAM: .* download complete" "${T}/work/logs/download_S5.log" \
  && { echo "FAIL: staged CRAM must not be re-downloaded"; exit 1; }
grep -q "MD5 verified" "${T}/work/logs/download_S5.log" \
  || { echo "FAIL: staged CRAM must still be verified"; exit 1; }
echo "PASS: staged S5 verified and dispatched without downloading, despite MAX_LOCAL_CRAMS=1"
check "S5 healed via staging" test -s "${T}/work/mosdepth_output/S5.by1000.regions.bed.gz"

echo "=== Sweep 3: nothing left ==="
run_manager > "${T}/sweep3.log" 2>&1 || { echo "FAIL: sweep 3 should succeed"; exit 1; }
grep -q "Nothing to download." "${T}/sweep3.log" || { echo "FAIL: sweep 3 should be a no-op"; exit 1; }
echo "PASS: sweep 3 is a no-op"

echo "=== Generator: a manifest row with whitespace in a field ==="
mkdir -p "${T}/work2/crams"
printf 'SAMPLE\tCRAM\tCRAI\tMD5\n' > "${T}/work2/manifest.tsv"
printf 'G1\tftp://ftp.sra.ebi.ac.uk/vol1/run/E/G1.final.cram\tftp://ftp.sra.ebi.ac.uk/vol1/run/E/G1.final.cram.crai\tabc\n' >> "${T}/work2/manifest.tsv"
printf 'B AD\tftp://ftp.sra.ebi.ac.uk/vol1/run/E/BAD.final.cram\tftp://ftp.sra.ebi.ac.uk/vol1/run/E/BAD.final.cram.crai\tdef\n' >> "${T}/work2/manifest.tsv"
batch2="$(WORK_DIR="${T}/work2" MANIFEST="${T}/work2/manifest.tsv" COMPARE_FAST_MODE=0 \
  bash "${HERE}/globus_batch_from_manifest.sh" 2> "${T}/generator2.log")"
check "good row emitted" test "$(echo "${batch2}" | grep -c "G1.cram")" = "2"
echo "${batch2}" | grep -q "BAD" && { echo "FAIL: whitespace row must not be emitted"; exit 1; }
grep -q "skipping manifest line 3: whitespace inside a field (sample 'B AD')" "${T}/generator2.log" \
  || { echo "FAIL: whitespace warning should name the manifest line"; cat "${T}/generator2.log"; exit 1; }
echo "PASS: whitespace row skipped with a warning naming manifest line 3"

echo ""
echo "=== All download-manager tests passed ==="
