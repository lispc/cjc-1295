#!/bin/bash
# Quick catalytic geometry analysis for WT and D-Ala2 MD
set -e
source /home/scroll/miniforge3/bin/activate gmx
cd /home/scroll/personal/cjc-1295/workspace/step3

SUMMARY="catalytic_summary.txt"
echo "Catalytic Geometry Analysis — $(date)" > $SUMMARY
echo "============================================================" | tee -a $SUMMARY

analyze_xvg() {
    local f=$1 label=$2 unit=$3
    awk -v label="$label" -v unit="$unit" '
    !/^[@#]/ && NF>1 {
        val=$2; sum+=val; sumsq+=val*val; n++
    } END {
        if(n>0) {
            mean=sum/n; std=sqrt(sumsq/n - mean*mean)
            printf "  %-35s %6.2f ± %5.2f %s (n=%d)\n", label, mean, std, unit, n
        }
    }' $f | tee -a $SUMMARY
}

echo ""
echo "=== WT MD: Ser630 OG -> Ala2 C ==="
gmx distance -s md.tpr -f md.trr md.part0002.trr \
    -o dist_wt_d1.xvg -len 1000 \
    -select 'atom OG and resnr 630 plus atom C and resnr 2' 2>&1 | tail -1

echo "=== WT MD: Tyr1 N -> Glu205 OE1 ==="
gmx distance -s md.tpr -f md.trr md.part0002.trr \
    -o dist_wt_d2.xvg -len 1000 \
    -select 'atom N and resnr 1 plus atom OE1 and resnr 205' 2>&1 | tail -1

echo "=== WT MD: Tyr1 N -> Glu206 OE1 ==="
gmx distance -s md.tpr -f md.trr md.part0002.trr \
    -o dist_wt_d3.xvg -len 1000 \
    -select 'atom N and resnr 1 plus atom OE1 and resnr 206' 2>&1 | tail -1

echo "=== WT MD: Ala2 O -> Tyr547 OH ==="
gmx distance -s md.tpr -f md.trr md.part0002.trr \
    -o dist_wt_d4.xvg -len 1000 \
    -select 'atom O and resnr 2 plus atom OH and resnr 547' 2>&1 | tail -1

echo "=== WT MD: Attack angle OG-C-N ==="
echo 'atom OG and resnr 630 atom C and resnr 2 atom N and resnr 2' | \
    gmx angle -type angle -s md.tpr -f md.trr md.part0002.trr -ov ang_wt.xvg 2>&1 | tail -1

echo ""
echo "--- WT MD Results (53.5 ns) ---"
analyze_xvg dist_wt_d1.xvg "Ser630 OG -> Ala2 C" "Å"
analyze_xvg dist_wt_d2.xvg "Tyr1 N -> Glu205 OE1" "Å"
analyze_xvg dist_wt_d3.xvg "Tyr1 N -> Glu206 OE1" "Å"
analyze_xvg dist_wt_d4.xvg "Ala2 O -> Tyr547 OH" "Å"
awk '!/^[@#]/ && NF>1 {v=$2; s+=v; s2+=v*v; n++} END {m=s/n; sd=sqrt(s2/n-m*m); printf "  %-35s %6.1f ± %5.1f %s (n=%d)\n", "Attack angle OG-C-N", m, sd, "deg", n}' ang_wt.xvg | tee -a $SUMMARY

echo ""
echo "=== D-Ala2 MD: Ser630 OG -> Ala2 C ==="
gmx distance -s DAla2_production.tpr -f DAla2_md.xtc \
    -o dist_dala2_d1.xvg -len 1000 \
    -select 'atom OG and resnr 630 plus atom C and resnr 2' 2>&1 | tail -1

echo "=== D-Ala2 MD: Tyr1 N -> Glu205 OE1 ==="
gmx distance -s DAla2_production.tpr -f DAla2_md.xtc \
    -o dist_dala2_d2.xvg -len 1000 \
    -select 'atom N and resnr 1 plus atom OE1 and resnr 205' 2>&1 | tail -1

echo "=== D-Ala2 MD: Tyr1 N -> Glu206 OE1 ==="
gmx distance -s DAla2_production.tpr -f DAla2_md.xtc \
    -o dist_dala2_d3.xvg -len 1000 \
    -select 'atom N and resnr 1 plus atom OE1 and resnr 206' 2>&1 | tail -1

echo "=== D-Ala2 MD: Ala2 O -> Tyr547 OH ==="
gmx distance -s DAla2_production.tpr -f DAla2_md.xtc \
    -o dist_dala2_d4.xvg -len 1000 \
    -select 'atom O and resnr 2 plus atom OH and resnr 547' 2>&1 | tail -1

echo "=== D-Ala2 MD: Attack angle OG-C-N ==="
echo 'atom OG and resnr 630 atom C and resnr 2 atom N and resnr 2' | \
    gmx angle -type angle -s DAla2_production.tpr -f DAla2_md.xtc -ov ang_dala2.xvg 2>&1 | tail -1

echo ""
echo "--- D-Ala2 MD Results (17.9 ns) ---"
analyze_xvg dist_dala2_d1.xvg "Ser630 OG -> Ala2 C" "Å"
analyze_xvg dist_dala2_d2.xvg "Tyr1 N -> Glu205 OE1" "Å"
analyze_xvg dist_dala2_d3.xvg "Tyr1 N -> Glu206 OE1" "Å"
analyze_xvg dist_dala2_d4.xvg "Ala2 O -> Tyr547 OH" "Å"
awk '!/^[@#]/ && NF>1 {v=$2; s+=v; s2+=v*v; n++} END {m=s/n; sd=sqrt(s2/n-m*m); printf "  %-35s %6.1f ± %5.1f %s (n=%d)\n", "Attack angle OG-C-N", m, sd, "deg", n}' ang_dala2.xvg | tee -a $SUMMARY

echo ""
echo "============================================================"
echo "  COMPARISON: Start vs Current"
echo "============================================================"
echo ""
printf "%-30s %12s %12s %12s %12s\n" "Metric" "Start" "WT(53.5ns)" "D-Ala2(17.9ns)" "Ideal"
printf "%-30s %12s %12s %12s %12s\n" "------------------------------" "--------" "--------" "--------" "--------"

# Get means
WT1=$(grep "Ser630 OG" $SUMMARY | head -1 | awk '{print $(NF-3)}')
WT2=$(grep "Tyr1 N -> Glu205" $SUMMARY | head -1 | awk '{print $(NF-3)}')
WT3=$(grep "Tyr1 N -> Glu206" $SUMMARY | head -1 | awk '{print $(NF-3)}')
WT4=$(grep "Ala2 O -> Tyr547" $SUMMARY | head -1 | awk '{print $(NF-3)}')
WTA=$(grep "Attack angle" $SUMMARY | head -1 | awk '{print $(NF-3)}')

D1=$(grep "Ser630 OG" $SUMMARY | tail -1 | awk '{print $(NF-3)}')
D2=$(grep "Tyr1 N -> Glu205" $SUMMARY | tail -1 | awk '{print $(NF-3)}')
D3=$(grep "Tyr1 N -> Glu206" $SUMMARY | tail -1 | awk '{print $(NF-3)}')
D4=$(grep "Ala2 O -> Tyr547" $SUMMARY | tail -1 | awk '{print $(NF-3)}')
DA=$(grep "Attack angle" $SUMMARY | tail -1 | awk '{print $(NF-3)}')

printf "%-30s %12s %12s %12s %12s\n" "Ser630 OG -> Ala2 C" "2.82 Å" "${WT1} Å" "${D1} Å" "2.5-4.0 Å"
printf "%-30s %12s %12s %12s %12s\n" "Tyr1 N -> Glu205 OE1" "2.02 Å" "${WT2} Å" "${D2} Å" "< 4.0 Å"
printf "%-30s %12s %12s %12s %12s\n" "Tyr1 N -> Glu206 OE1" "~4.5 Å" "${WT3} Å" "${D3} Å" "< 4.5 Å"
printf "%-30s %12s %12s %12s %12s\n" "Ala2 O -> Tyr547 OH" "3.10 Å" "${WT4} Å" "${D4} Å" "< 4.0 Å"
printf "%-30s %12s %12s %12s %12s\n" "Attack angle OG-C-N" "113.6°" "${WTA}°" "${DA}°" "80-120°"

echo ""
echo "Done."
