"""Reproduce every computed table, statistic and data figure of the manuscript.

    python code/run_all.py            # about 2-3 minutes on a laptop
    python code/run_all.py --timing   # also the run-time experiment of Section 5.3
"""
import sys
import time
import build_dataset, analysis_alkanes, analysis_dataset, boiling_point, sensitivity, networks, verify_theory, make_figures, timing

STEPS = [("Data set and descriptors (data/)", build_dataset.main),
         ("Table 3 and Section 7.6", analysis_dataset.main),
         ("Tables 5-8", analysis_alkanes.main),
         ("Table 9", boiling_point.main),
         ("Table 10", sensitivity.main),
         ("Table 11", networks.main),
         ("Section 6 checks", verify_theory.main),
         ("Figures 10-16", make_figures.main)]

if __name__ == "__main__":
    if "--timing" in sys.argv:
        STEPS.append(("Section 5.3 timing", timing.main))
    t = time.time()
    for name, fn in STEPS:
        print(f"\n=== {name} ===")
        fn()
    print(f"\nDone in {time.time() - t:.0f} s. Results are in results/ and figures/.")
