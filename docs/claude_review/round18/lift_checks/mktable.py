import json
rows = ["| `N` | `w` | contacts with no divisor in window | contacts used | splits tested | max skew / `w` | mean skew / `w` | mover `a`: `h/w`, `log lambda_1/w`, index`/w`, skew`/w` |",
        "|---|---|---|---|---|---|---|---|"]
for r in json.load(open('intlift.json')):
    if r['N'] < 3: continue
    rows.append("| %d | %.1f | %d | %d | %d | %.3f | %.3f | %.2f, %.2f, %.2f, %.2f |" % (
        r['N'], r['w'], r['nonsplittable_sides'], r['n_sides_used'], r['n_tested'], r['max_skew_over_w'],
        r['mean_skew_over_w'], r['mover_a']['h_a_over_w'], r['mover_a']['lambda1_over_w'], r['mover_a']['index_over_w'],
        r['mover_a_skew_over_w']))
print("\n".join(rows))
