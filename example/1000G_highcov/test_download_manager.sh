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
# fails unless FAKE_S3 is set and holds the requested basename; the URL is the
# last argument and --dir/--out name the destination, like the real thing
[[ -n "${FAKE_S3:-}" ]] || exit 1
dir=""; out=""; url=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dir) dir="$2"; shift 2 ;;
    --out) out="$2"; shift 2 ;;
    -*) shift ;;
    *) url="$1"; shift ;;
  esac
done
src="${FAKE_S3}/$(basename "${url}")"
[[ -f "${src}" ]] || exit 1
cp "${src}" "${dir}/${out}"
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

run_generator() { # <output_prefix>
  WORK_DIR="${T}/work" \
  MANIFEST="${T}/work/manifest.tsv" \
  COMPARE_FAST_MODE=1 \
  bash "${HERE}/globus_batch_from_manifest.sh" "$1"
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
run_generator "${T}/gbatch" 2> "${T}/generator.log"
grep -q "^\"/vol1/run/ERR5/S5.final.cram\" \"${T}/work/crams/S5.cram\"$" "${T}/gbatch_cram.txt" \
  || { echo "FAIL: CRAM batch line"; cat "${T}/gbatch_cram.txt"; exit 1; }
grep -q "^\"/vol1/run/ERR5/S5.final.cram.crai\" \"${T}/work/crams/S5.cram.crai\"$" "${T}/gbatch_crai.txt" \
  || { echo "FAIL: CRAI batch line"; cat "${T}/gbatch_crai.txt"; exit 1; }
check "CRAI batch holds only the needed sample" test "$(wc -l < "${T}/gbatch_crai.txt" | tr -d '[:space:]')" = "1"
check "CRAM batch holds only the needed sample" test "$(wc -l < "${T}/gbatch_cram.txt" | tr -d '[:space:]')" = "1"
grep -q "Emitted 1 samples" "${T}/generator.log" || { echo "FAIL: generator summary"; exit 1; }
echo "PASS: generator splits S5 into CRAI and CRAM batches, quoted and renamed, skips the done four"

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
WORK_DIR="${T}/work2" MANIFEST="${T}/work2/manifest.tsv" COMPARE_FAST_MODE=0 \
  bash "${HERE}/globus_batch_from_manifest.sh" "${T}/gbatch2" 2> "${T}/generator2.log"
check "good row in the CRAI batch" test "$(grep -c "G1.cram.crai" "${T}/gbatch2_crai.txt")" = "1"
check "good row in the CRAM batch" test "$(grep -c "G1.cram\"" "${T}/gbatch2_cram.txt")" = "1"
grep -q "BAD" "${T}/gbatch2_crai.txt" "${T}/gbatch2_cram.txt" \
  && { echo "FAIL: whitespace row must not be emitted"; exit 1; }
grep -q "skipping manifest line 3: whitespace inside a field (sample 'B AD')" "${T}/generator2.log" \
  || { echo "FAIL: whitespace warning should name the manifest line"; cat "${T}/generator2.log"; exit 1; }
echo "PASS: whitespace row skipped with a warning naming manifest line 3"

echo "=== S3 mirror: URL mapping and source preference ==="
sed -n '/^s3_url_for() {/,/^}/p' "${HERE}/01_download_and_mosdepth.sh" > "${T}/s3fn.sh"
(
  source "${T}/s3fn.sh"
  S3_HTTPS_BASE="https://1000genomes.s3.amazonaws.com"
  u="ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3239480/NA12718.final.cram"
  [[ "$(s3_url_for "${u}" 2504)" == "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/ERR3239480/NA12718.final.cram" ]] \
    || { echo "FAIL: 2504 mapping"; exit 1; }
  [[ "$(s3_url_for "${u}.crai" 698)" == "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3239480/NA12718.final.cram.crai" ]] \
    || { echo "FAIL: 698 mapping"; exit 1; }
  [[ -z "$(s3_url_for "${u}" "")" ]] || { echo "FAIL: unknown batch must not map"; exit 1; }
  [[ -z "$(s3_url_for "ftp://elsewhere.org/x.cram" 2504)" ]] || { echo "FAIL: foreign URL must not map"; exit 1; }
  S3_HTTPS_BASE=""
  [[ -z "$(s3_url_for "${u}" 2504)" ]] || { echo "FAIL: empty base must disable"; exit 1; }
) || exit 1
echo "PASS: s3_url_for maps both batches, refuses the unmappable, and can be disabled"

mkdir -p "${T}/work5/crams" "${T}/s3store"
python3 - "${T}" <<'EOF'
import hashlib, sys
T = sys.argv[1]
content = b"x1-bytes" * 120
open(f"{T}/s3store/X1.final.cram", "wb").write(content)
open(f"{T}/s3store/X1.final.cram.crai", "wb").write(b"cx1")
base = "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERRX/ERRX123/X1.final.cram"
with open(f"{T}/work5/manifest.tsv", "w") as out:
    out.write("SAMPLE\tCRAM\tCRAI\tMD5\tRELEASE_BATCH\n")
    out.write(f"X1\t{base}\t{base}.crai\t{hashlib.md5(content).hexdigest()}\t2504\n")
EOF
WORK_DIR="${T}/work5" \
SIF_IMAGE="${T}/work/ngs-pca.sif" \
REF_DIR="${T}/work/reference" \
REF_FASTA="${T}/work/reference/ref.fa" \
MANIFEST="${T}/work5/manifest.tsv" \
MIN_MANIFEST_SAMPLES=1 EXPECTED_MANIFEST_SAMPLES=1 \
USE_ASPERA=0 COMPARE_FAST_MODE=1 DOWNLOADER_LOCAL=1 DOWNLOAD_SLOTS=1 \
FAKE_S3="${T}/s3store" FAKE_REMOTE="${T}/remote" \
bash "${HERE}/01_download_and_mosdepth.sh" > "${T}/s3run.log" 2>&1 \
  || { echo "FAIL: S3-served sweep should succeed"; cat "${T}/s3run.log" "${T}/work5/logs/download_X1.log" 2>/dev/null; exit 1; }
grep -q "aria2c download complete (S3," "${T}/work5/logs/download_X1.log" \
  || { echo "FAIL: S3 should have served the download"; cat "${T}/work5/logs/download_X1.log"; exit 1; }
check "S3-served sample reached both trees" test -s "${T}/work5/mosdepth_output/X1.by1000.regions.bed.gz" \
  -a -s "${T}/work5/mosdepth_output_fast/X1.by1000.regions.bed.gz"
echo "PASS: S3 mirror preferred and served end to end"

echo "=== Stage watcher: dispatch on arrival, re-dispatch after premature MD5, park the hopeless ==="
mkdir -p "${T}/work3/crams"
python3 - "${T}" <<'EOF'
import hashlib, sys
T = sys.argv[1]
rows = []
for sample, content in (("W1", b"w1-bytes" * 100), ("W2", b"w2-bytes" * 100), ("W3", b"w3-bytes" * 100)):
    base = f"ftp://ftp.sra.ebi.ac.uk/vol1/run/E/{sample}.final.cram"
    md5 = hashlib.md5(content).hexdigest()
    if sample == "W3":
        md5 = hashlib.md5(b"never what lands on disk").hexdigest()
    rows.append(f"{sample}\t{base}\t{base}.crai\t{md5}\tx")
with open(f"{T}/work3/manifest.tsv", "w") as out:
    out.write("SAMPLE\tCRAM\tCRAI\tMD5\tNOTES\n")
    out.write("\n".join(rows) + "\n")
# W1 fully staged before the watcher starts; W3 staged but forever wrong
open(f"{T}/work3/crams/W1.cram", "wb").write(b"w1-bytes" * 100)
open(f"{T}/work3/crams/W1.cram.crai", "wb").write(b"c1")
open(f"{T}/work3/crams/W3.cram", "wb").write(b"corrupt forever")
open(f"{T}/work3/crams/W3.cram.crai", "wb").write(b"c3")
EOF

# W2 arrives mid-run: partial first, completed a few seconds later
(
  sleep 2
  python3 -c "open('${T}/work3/crams/W2.cram','wb').write(b'w2-byt')"
  python3 -c "open('${T}/work3/crams/W2.cram.crai','wb').write(b'c2')"
  sleep 3
  python3 -c "open('${T}/work3/crams/W2.cram','wb').write(b'w2-bytes' * 100)"
) &
STAGER_PID=$!

watch_rc=0
WORK_DIR="${T}/work3" \
SIF_IMAGE="${T}/work/ngs-pca.sif" \
REF_DIR="${T}/work/reference" \
REF_FASTA="${T}/work/reference/ref.fa" \
MANIFEST="${T}/work3/manifest.tsv" \
COMPARE_FAST_MODE=1 \
WATCHER_LOCAL=1 \
WATCH_INTERVAL=1 \
WATCH_IDLE_EXIT=0 \
WATCH_DISPATCH_ATTEMPTS=2 \
bash "${HERE}/01c_dispatch_staged.sh" > "${T}/watch.log" 2>&1 || watch_rc=$?
wait "${STAGER_PID}" 2>/dev/null || true

check "watcher exits non-zero when a sample is parked" test "${watch_rc}" -ne 0
grep -q "Fast-mode comparison: ON" "${T}/watch.log" \
  || { echo "FAIL: watcher must announce the comparison mode"; exit 1; }
echo "PASS: watcher announces the comparison is ON"
check "W1 dispatched and completed" test -s "${T}/work3/mosdepth_output/W1.by1000.regions.bed.gz"
check "W2 completed after arriving mid-run" test -s "${T}/work3/mosdepth_output/W2.by1000.regions.bed.gz"
check "W2 fast tree too" test -s "${T}/work3/mosdepth_output_fast/W2.by1000.regions.bed.gz"
check "W3 parked" test -f "${T}/work3/download_state/watch/W3.parked"
check "W3 produced no output" test ! -e "${T}/work3/mosdepth_output/W3.by1000.regions.bed.gz"
check "W3's staged file was not deleted" test -s "${T}/work3/crams/W3.cram"
check "completed samples cleaned their CRAMs" test ! -e "${T}/work3/crams/W1.cram"
check "W2's CRAM and CRAI cleaned after its verified run" \
  test ! -e "${T}/work3/crams/W2.cram" -a ! -e "${T}/work3/crams/W2.cram.crai"
grep -q "parking W3" "${T}/watch.log" || { echo "FAIL: parking warning missing"; exit 1; }
echo "PASS: watcher log explains the parking"

echo "=== Stage watcher: a growing unpaired file holds the idle clock open ==="
mkdir -p "${T}/work4/crams"
printf 'SAMPLE\tCRAM\tCRAI\tMD5\n' > "${T}/work4/manifest.tsv"
printf 'V1\tftp://ftp.sra.ebi.ac.uk/vol1/run/E/V1.final.cram\tftp://ftp.sra.ebi.ac.uk/vol1/run/E/V1.final.cram.crai\tabc\n' >> "${T}/work4/manifest.tsv"

# a CRAM that grows for ~6 s and never gets its CRAI - the CRAM-first phase
(
  for i in 1 2 3 4 5 6; do
    python3 -c "open('${T}/work4/crams/V1.cram','ab').write(b'chunk')"
    sleep 1
  done
) &
GROWER_PID=$!

T0=$(date +%s)
idle_rc=0
WORK_DIR="${T}/work4" \
MANIFEST="${T}/work4/manifest.tsv" \
COMPARE_FAST_MODE=1 \
WATCHER_LOCAL=1 \
WATCH_INTERVAL=1 \
WATCH_IDLE_EXIT=2 \
bash "${HERE}/01c_dispatch_staged.sh" > "${T}/idle.log" 2>&1 || idle_rc=$?
ELAPSED=$(( $(date +%s) - T0 ))
wait "${GROWER_PID}" 2>/dev/null || true

check "idle exit is clean when nothing is parked" test "${idle_rc}" -eq 0
check "growth held the 2 s idle clock open past the 6 s of writes" test "${ELAPSED}" -ge 7
grep -q "still missing - the transfer is finished or stalled" "${T}/idle.log" \
  || { echo "FAIL: idle-exit message missing"; cat "${T}/idle.log"; exit 1; }
echo "PASS: idle exit fired only after the directory went quiet"

echo ""
echo "=== All download-manager tests passed ==="
