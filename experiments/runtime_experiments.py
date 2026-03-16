# experiments/runtime_experiment.py
# Empirically validates O(mn) time complexity for both algorithms.
# Generates runtime data and saves a plot for the project report.
#
# Run from project root:
#   python experiments/runtime_experiment.py

import sys, os, time, random, csv
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.needleman_wunsch import needleman_wunsch
from src.smith_waterman   import smith_waterman

# ── Configuration ─────────────────────────────────────────────
ALPHABET = "ACGT"
SIZES    = [50, 100, 200, 300, 400, 500, 750, 1000, 1250, 1500]
REPEATS  = 3  # average over this many random pairs per size

def random_seq(n):
    return ''.join(random.choices(ALPHABET, k=n))

# ── Run experiments ────────────────────────────────────────────
print("Running runtime experiments...")
print(f"{'Size':>6}  {'NW (s)':>10}  {'SW (s)':>10}  {'n²':>10}")
print("-" * 42)

results = []

for n in SIZES:
    nw_total = 0.0
    sw_total = 0.0

    for _ in range(REPEATS):
        s1 = random_seq(n)
        s2 = random_seq(n)

        t0 = time.perf_counter()
        needleman_wunsch(s1, s2)
        nw_total += time.perf_counter() - t0

        t0 = time.perf_counter()
        smith_waterman(s1, s2)
        sw_total += time.perf_counter() - t0

    nw_avg = nw_total / REPEATS
    sw_avg = sw_total / REPEATS
    n_sq   = n * n
    results.append((n, nw_avg, sw_avg, n_sq))
    print(f"{n:>6}  {nw_avg:>10.4f}  {sw_avg:>10.4f}  {n_sq:>10}")

# ── Save CSV ───────────────────────────────────────────────────
os.makedirs("experiments/output", exist_ok=True)
csv_path = "experiments/output/runtime_results.csv"
with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["sequence_length", "nw_time_sec", "sw_time_sec", "n_squared"])
    writer.writerows(results)
print(f"\nCSV saved → {csv_path}")

# ── Plot ───────────────────────────────────────────────────────
try:
    import matplotlib.pyplot as plt

    ns       = [r[0] for r in results]
    nw_times = [r[1] for r in results]
    sw_times = [r[2] for r in results]
    n_sq     = [r[3] for r in results]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Empirical Runtime Validation: O(mn) Behavior",
                 fontsize=14, fontweight='bold')

    # ── Left: time vs n ───────────────────────────────────────
    axes[0].plot(ns, nw_times, 'o-', color='red',
                 label='Needleman-Wunsch', linewidth=2)
    axes[0].plot(ns, sw_times, 's-', color='green',
                 label='Smith-Waterman', linewidth=2)
    axes[0].set_xlabel("Sequence Length n")
    axes[0].set_ylabel("Average Runtime (seconds)")
    axes[0].set_title("Runtime vs Sequence Length")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # ── Right: time vs n² (linear = confirms O(n²)) ───────────
    axes[1].plot(n_sq, nw_times, 'o-', color='red',
                 label='Needleman-Wunsch', linewidth=2)
    axes[1].plot(n_sq, sw_times, 's-', color='green',
                 label='Smith-Waterman', linewidth=2)
    axes[1].set_xlabel("n² (Sequence Length Squared)")
    axes[1].set_ylabel("Average Runtime (seconds)")
    axes[1].set_title("Runtime vs n²  (Linear = confirms O(n²))")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plot_path = "experiments/output/runtime_plot.png"
    plt.savefig(plot_path, dpi=150)
    print(f"Plot saved  → {plot_path}")
    plt.show()

except ImportError:
    print("matplotlib not installed — run: pip install matplotlib")