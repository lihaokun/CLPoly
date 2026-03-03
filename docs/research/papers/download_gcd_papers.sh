#!/bin/bash
# Download GCD algorithm reference papers
# Run from any directory: bash docs/research/papers/download_gcd_papers.sh

DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

echo "=== Downloading GCD algorithm papers to $DIR ==="
echo ""

# 1. Zippel 1979 - Probabilistic algorithms for sparse polynomials (EUROSAM)
#    Freely available from Monagan's teaching page at SFU
echo "[1/5] Zippel 1979 - Probabilistic algorithms for sparse polynomials..."
wget -q --show-progress -O "Zippel-1979.pdf" \
  "http://www.cecm.sfu.ca/~monaganm/teaching/TopicsinCA15/zippel79.pdf" 2>&1
if file "Zippel-1979.pdf" | grep -q PDF; then
  echo "  -> OK"
else
  echo "  -> FAILED (not a valid PDF)"
fi
echo ""

# 2. Brown 1971 - On Euclid's algorithm and the computation of polynomial GCDs (JACM)
#    Freely available from Monagan's teaching page at SFU
echo "[2/5] Brown 1971 - On Euclid's algorithm and polynomial GCDs..."
wget -q --show-progress -O "Brown-Euclid-Algorithm-1971.pdf" \
  "http://www.cecm.sfu.ca/~monaganm/teaching/TopicsinCA15/Brown1971.pdf" 2>&1
if file "Brown-Euclid-Algorithm-1971.pdf" | grep -q PDF; then
  echo "  -> OK"
else
  echo "  -> FAILED (not a valid PDF)"
fi
echo ""

# 3. Thull & Yap 1990 - A Unified Approach to HGCD Algorithms
#    NYU CS Technical Report TR1990-505
echo "[3/5] Thull & Yap 1990 - A Unified Approach to HGCD..."
wget -q --show-progress -O "Thull-Yap-HGCD-1990.pdf" \
  "https://cs.nyu.edu/media/publications/TR1990-505.pdf" 2>&1
if file "Thull-Yap-HGCD-1990.pdf" | grep -q PDF; then
  echo "  -> OK"
else
  echo "  -> FAILED - trying Springer..."
  wget -q --show-progress -O "Thull-Yap-HGCD-1990.pdf" \
    "https://link.springer.com/content/pdf/10.1007/3-540-54458-5_76.pdf" 2>&1
  if file "Thull-Yap-HGCD-1990.pdf" | grep -q PDF; then
    echo "  -> OK (from Springer)"
  else
    echo "  -> FAILED (may need institutional access)"
  fi
fi
echo ""

# 4. Hu & Monagan 2021 - Speeding up polynomial GCD (Maple Transactions)
#    Open access journal
echo "[4/5] Hu & Monagan 2021 - Speeding up polynomial GCD..."
wget -q --show-progress -O "Hu-Monagan-Speeding-Up-Polynomial-GCD-2021.pdf" \
  "https://mapletransactions.org/index.php/maple/article/download/14452/11895/36425" 2>&1
if file "Hu-Monagan-Speeding-Up-Polynomial-GCD-2021.pdf" | grep -q PDF; then
  echo "  -> OK"
else
  echo "  -> FAILED (not a valid PDF)"
fi
echo ""

# 5. GCDHEU - Char, Geddes, Gonnet 1989 (JSC)
#    ScienceDirect is usually paywalled; try multiple sources
echo "[5/5] GCDHEU - Char, Geddes, Gonnet 1989..."
# Try ScienceDirect first (may need institutional access)
wget -q --show-progress -O "GCDHEU-Char-Geddes-Gonnet-1989.pdf" \
  "https://www.sciencedirect.com/science/article/pii/S0747717189800045/pdf" \
  --header="User-Agent: Mozilla/5.0" 2>&1
if file "GCDHEU-Char-Geddes-Gonnet-1989.pdf" | grep -q PDF; then
  echo "  -> OK (from ScienceDirect)"
else
  echo "  -> ScienceDirect failed, trying Springer (EUROSAM proceedings version)..."
  wget -q --show-progress -O "GCDHEU-Char-Geddes-Gonnet-1989.pdf" \
    "https://link.springer.com/content/pdf/10.1007/BFb0032851.pdf" \
    --header="User-Agent: Mozilla/5.0" 2>&1
  if file "GCDHEU-Char-Geddes-Gonnet-1989.pdf" | grep -q PDF; then
    echo "  -> OK (from Springer, EUROSAM proceedings)"
  else
    echo "  -> FAILED (paywalled - see manual download instructions below)"
    rm -f "GCDHEU-Char-Geddes-Gonnet-1989.pdf"
    echo ""
    echo "  Manual download options for GCDHEU:"
    echo "  a) ScienceDirect (with institutional access):"
    echo "     https://www.sciencedirect.com/science/article/pii/S0747717189800045"
    echo "  b) Springer (EUROSAM 1984 proceedings version):"
    echo "     https://link.springer.com/chapter/10.1007/BFb0032851"
    echo "  c) Academia.edu (may require account):"
    echo "     https://www.academia.edu/73430211/"
    echo "  d) ETH Zurich (Gonnet's page, algorithm description):"
    echo "     https://people.inf.ethz.ch/gonnet/CAII/HeuristicAlgorithms/node1.html"
  fi
fi
echo ""

echo "=== Download summary ==="
echo "Files in $DIR:"
ls -lh *.pdf 2>/dev/null | while read line; do echo "  $line"; done
echo ""
echo "Done."
