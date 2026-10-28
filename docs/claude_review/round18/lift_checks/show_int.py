import json
for r in json.load(open('intlift.json')):
    print(r['N'], 'w=%.1f'%r['w'], 'unsplit', r['nonsplittable_sides'], 'used', r.get('n_sides_used'), 'tested', r.get('n_tested'), 'min ell/w %.2f'%r.get('min_ell_over_w',0), 'mean %.2f'%r.get('mean_ell_over_w',0), 'idx/w %.2f'%r.get('mean_index_over_w',0), 'mover a: h/w %.2f lam1/w %.2f idx/w %.2f'%(r['mover_a']['h_a_over_w'], r['mover_a']['lambda1_over_w'], r['mover_a']['index_over_w']), 'ft %.1f'%r['factor_time'])
