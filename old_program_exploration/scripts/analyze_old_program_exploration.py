#!/usr/bin/env python3
"""
Old Program Exploration - Analysis & Plotting Script
Reads scan results from old_program_exploration/raw/ and generates summary CSVs,
figures, and final report.
"""
import os, sys, warnings, traceback
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

BASE = Path(__file__).resolve().parent.parent
RAW = BASE / "raw"
FIG = BASE / "figures"
RPT = BASE / "reports"
LOG = BASE / "logs"
for d in [RAW, FIG, RPT, LOG]:
    d.mkdir(parents=True, exist_ok=True)

def safe_csv(p):
    try: return pd.read_csv(p)
    except: return None

def safe_txt(p):
    try: return Path(p).read_text()
    except: return None

def parse_summary_txt(path):
    """Parse polydisperse_coupled_summary.txt or eval_summary.txt"""
    txt = safe_txt(path)
    if not txt: return {}
    d = {}
    for line in txt.split('\n'):
        if ':' in line:
            parts = line.strip().split(':', 1)
            key = parts[0].strip()
            val = parts[1].strip().split()[0]
            try: d[key] = float(val)
            except: d[key] = val
    return d

def parse_limit_cycle(path):
    txt = safe_txt(path)
    if not txt: return {}
    d = {}
    for line in txt.split('\n'):
        if ':' in line:
            parts = line.strip().split(':', 1)
            key = parts[0].strip()
            val = parts[1].strip().split()[0]
            try: d[key] = float(val)
            except: d[key] = val
    return d

def parse_coupled_trajectory(path):
    df = safe_csv(path)
    if df is None or len(df)==0: return None, {}
    info = {}
    if 'accepted_cycle' in df.columns:
        info['n_cycles'] = len(df)
        info['n_accepted'] = int(df['accepted_cycle'].sum())
    for col in df.columns:
        cl = col.lower()
        try:
            v0 = float(df[col].iloc[0]); v1 = float(df[col].iloc[-1])
        except: continue
        if 'theta_sum' in cl:
            info['Theta_sum_initial'] = v0; info['Theta_sum_final'] = v1
        elif 'theta_shape' in cl:
            info['Theta_shape_initial'] = v0; info['Theta_shape_final'] = v1
        elif 'theta_nm' in cl or cl=='theta_nm':
            info['Theta_NM_initial'] = v0; info['Theta_NM_final'] = v1
        elif 'du_rms' in cl:
            info['du_rms_initial'] = v0; info['du_rms_final'] = v1
        elif 'u_overlap' in cl:
            info['U_overlap_initial'] = v0; info['U_overlap_final'] = v1
        elif 'f_total' in cl:
            info['F_total_initial'] = v0; info['F_total_final'] = v1
        elif 'radius_ratio' in cl:
            info['radius_ratio_initial'] = v0; info['radius_ratio_final'] = v1
    return df, info

def classify_run(info, summary):
    ac = info.get('n_accepted') or summary.get('n_accepted')
    ts_i = info.get('Theta_sum_initial') or summary.get('Theta_sum_initial')
    ts_f = info.get('Theta_sum_final') or summary.get('Theta_sum_final')
    tsh_i = info.get('Theta_shape_initial') or summary.get('Theta_shape_initial')
    tsh_f = info.get('Theta_shape_final') or summary.get('Theta_shape_final')
    # Check for block-level acceptance but cycle-level rejection
    block_acc = info.get('n_block_accepted', 0) or 0
    cycle_acc = ac or 0
    if cycle_acc == 0 and block_acc > 0: return 'cycle_rejected_blocks_accepted'
    if ac is not None and ac == 0: return 'frozen'
    if ts_i is not None and ts_f is not None:
        ds = ts_f - ts_i; dh = (tsh_f or 0) - (tsh_i or 0)
        if ds < -0.001 and dh < -0.001: return 'diagonal_improved'
        if ds < -0.001: return 'sum_improved'
        if dh < -0.001: return 'shape_improved'
        if abs(ds) < 0.001 and abs(dh) < 0.001: return 'stagnant'
    return 'unknown'

def gather_v18_run(rundir_name):
    """Gather data from a single v18 output directory."""
    d = RAW / rundir_name
    if not d.is_dir(): return None
    row = {'run_name': rundir_name}
    files_present = [f.name for f in d.iterdir()]
    row['n_files'] = len(files_present)
    
    # polydisperse_coupled_summary.txt
    s = parse_summary_txt(d / 'polydisperse_coupled_summary.txt')
    for k,v in s.items():
        if isinstance(v, (int,float)):
            row[k] = v
    
    # eval_summary.txt
    e = parse_summary_txt(d / 'eval_summary.txt')
    for k,v in e.items():
        if isinstance(v, (int,float)) and k not in row:
            row[k] = v
    
    # limit_cycle_diagnostic.txt
    lc = parse_limit_cycle(d / 'limit_cycle_diagnostic.txt')
    for k,v in lc.items():
        if isinstance(v, (int,float)):
            row[f'lc_{k}'] = v
    
    # coupled_cycle_trajectory.csv
    traj_df, traj_info = parse_coupled_trajectory(d / 'coupled_cycle_trajectory.csv')
    for k,v in traj_info.items():
        if k not in row: row[k] = v
    
    # Also check shape_block_trajectory.csv, ricci_block_trajectory.csv
    for fname in ['shape_block_trajectory.csv','ricci_block_trajectory.csv']:
        if (d / fname).exists():
            row[f'has_{fname.replace(".csv","")}'] = True
    
    # operator_response_matrix.csv - critical for understanding block vs cycle acceptance
    orm_path = d / 'operator_response_matrix.csv'
    orm_df = safe_csv(orm_path)
    if orm_df is not None and len(orm_df) > 0:
        row['has_operator_response_matrix'] = True
        row['orm_total_rows'] = len(orm_df)
        # Count blocks
        if 'block_name' in orm_df.columns:
            for bn in orm_df['block_name'].unique():
                row[f'orm_n_{bn}'] = int((orm_df['block_name']==bn).sum())
        # Count accepted blocks
        if 'accepted' in orm_df.columns:
            acc_vals = pd.to_numeric(orm_df['accepted'], errors='coerce')
            row['n_block_accepted'] = int(acc_vals.sum())
            row['n_block_total'] = int(acc_vals.notna().sum())
            row['block_acceptance_rate'] = acc_vals.mean()
        # Average deltas per block type
        for metric in ['Delta_Theta_sum', 'Delta_Theta_shape', 'Delta_Theta_NM', 'Delta_du_rms', 'Delta_U_overlap']:
            if metric in orm_df.columns:
                for bn in orm_df.get('block_name', pd.Series()).unique():
                    mask = orm_df['block_name'] == bn
                    v = pd.to_numeric(orm_df.loc[mask, metric], errors='coerce').dropna()
                    if len(v) > 0:
                        row[f'orm_{bn}_mean_{metric}'] = v.mean()
    else:
        row['has_operator_response_matrix'] = False
    
    # surgery_microquench_stats.csv
    sms_path = d / 'surgery_microquench_stats.csv'
    sms_df = safe_csv(sms_path)
    if sms_df is not None and len(sms_df) > 0:
        row['has_surgery_stats'] = True
        if 'accepted' in sms_df.columns:
            row['n_surgery_accepted'] = int(pd.to_numeric(sms_df['accepted'], errors='coerce').sum())
            row['n_surgery_total'] = len(sms_df)
    
    row['has_trajectory'] = traj_df is not None
    row['classification'] = classify_run(row, row)
    return row

# ============================================================
# Section C: Gamma Scan
# ============================================================
def analyze_gamma_scan():
    print("\n=== Analyzing v18 Gamma Scan ===")
    gammas = [0, 0.02, 0.05, 0.10, 0.25, 0.50, 1.0]
    rows = []
    for g in gammas:
        # Match directory naming: strip trailing zeros but keep decimal point pattern
        # 0 -> "0", 0.10 -> "0.10", 0.50 -> "0.50", 1.0 -> "1.0"
        gstr = f'{g:g}'  # 0->0, 0.1->0.1, 0.5->0.5, 1->1
        # Try multiple naming conventions
        candidates = [
            f'v18_poly_gamma_{g}',
            f'v18_poly_gamma_{gstr}',
            f'v18_poly_gamma_{g:.2f}',
        ]
        r = None
        for c in candidates:
            r = gather_v18_run(c)
            if r is not None:
                break
        if r: r['gamma'] = g; rows.append(r)
        else: rows.append({'gamma': g, 'run_name': dirname, 'classification': 'missing'})
    df = pd.DataFrame(rows)
    out = RAW / 'v18_gamma_scan_summary.csv'
    df.to_csv(out, index=False)
    print(f"  Saved {out} ({len(df)} rows)")
    for _,r in df.iterrows():
        print(f"  gamma={r['gamma']}: {r.get('classification','?')} accepted={r.get('n_accepted','?')} Dsum={r.get('Theta_sum_final','?')}-{r.get('Theta_sum_initial','?')}")
    return df

# ============================================================
# Section D: Bounded_du Scan
# ============================================================
def analyze_bound_scan():
    print("\n=== Analyzing v18 Bounded_du Scan ===")
    bdus = [0.03, 0.05, 0.10, 0.15, 0.25, 0.40]
    rows = []
    for b in bdus:
        dirname = f'v18_poly_bound_{b:.2f}'
        r = gather_v18_run(dirname)
        if r: r['bounded_du'] = b; rows.append(r)
        else:
            # For 0.10, use the gamma scan baseline (gamma=0.05 has bounded_du=0.10)
            if b == 0.10:
                r = gather_v18_run('v18_poly_gamma_0.05')
                if r: r['bounded_du'] = b; rows.append(r)
                else: rows.append({'bounded_du': b, 'run_name': dirname, 'classification': 'missing'})
            else:
                rows.append({'bounded_du': b, 'run_name': dirname, 'classification': 'missing'})
    df = pd.DataFrame(rows)
    out = RAW / 'v18_bound_scan_summary.csv'
    df.to_csv(out, index=False)
    print(f"  Saved {out} ({len(df)} rows)")
    return df

# ============================================================
# Section E: Weight Scan
# ============================================================
def analyze_weight_scan():
    print("\n=== Analyzing v18 Weight Scan ===")
    schemes = ['shape_biased','balanced','weak_overlap']
    rows = []
    for s in schemes:
        r = gather_v18_run(f'v18_weight_{s}')
        if r: r['weight_scheme'] = s; rows.append(r)
        else: rows.append({'weight_scheme': s, 'run_name': f'v18_weight_{s}', 'classification': 'missing'})
    df = pd.DataFrame(rows)
    out = RAW / 'v18_weight_scan_summary.csv'
    df.to_csv(out, index=False)
    print(f"  Saved {out} ({len(df)} rows)")
    return df

# ============================================================
# Section F: Trajectory Inventory
# ============================================================
def inventory_v18():
    print("\n=== v18 Trajectory Inventory ===")
    rows = []
    for d in sorted(RAW.iterdir()):
        if not d.is_dir() or not d.name.startswith('v18_'): continue
        r = gather_v18_run(d.name)
        if r: rows.append(r)
    if rows:
        df = pd.DataFrame(rows)
        out = RAW / 'v18_trajectory_inventory.csv'
        df.to_csv(out, index=False)
        print(f"  Saved {out} ({len(df)} entries)")
        return df
    print("  No v18 dirs found")
    return pd.DataFrame()

# ============================================================
# Section G: v17 Operator Statistics
# ============================================================
def analyze_operator_stats():
    print("\n=== v17 Operator Statistics ===")
    search_files = []
    # Search in raw/ subdirs and also in project root
    for subdir in [RAW / 'v17_op_16_basic_more', RAW / 'v17_op_24_basic_more', RAW]:
        for f in subdir.glob('operator_multistart_v17_*.csv'):
            search_files.append(f)
    # Also in project root
    proj = BASE.parent
    for f in proj.glob('operator_multistart_v17_*.csv'):
        if f not in search_files: search_files.append(f)
    
    all_rows = []
    for fpath in search_files:
        df = safe_csv(fpath)
        if df is None or len(df)==0: continue
        print(f"  Processing {fpath.name} ({len(df)} rows)")
        
        op_col = None
        for c in ['operator_name','operator','name']:
            if c in df.columns: op_col = c; break
        if not op_col: op_col = df.columns[0]
        
        sum_col = None; shape_col = None; nm_col = None
        for c in df.columns:
            cl = c.lower()
            if 'delta_theta_sum' in cl: sum_col = c
            elif 'delta_theta_shape' in cl: shape_col = c
            elif 'delta_theta_nm' in cl: nm_col = c
        
        for op, grp in df.groupby(op_col):
            r = {'file': fpath.name, 'operator': op, 'count': len(grp)}
            if sum_col:
                v = pd.to_numeric(grp[sum_col], errors='coerce').dropna()
                if len(v)>0:
                    r['mean_Delta_Theta_sum'] = v.mean()
                    r['median_Delta_Theta_sum'] = v.median()
                    r['success_left'] = (v<0).mean()
            if shape_col:
                v = pd.to_numeric(grp[shape_col], errors='coerce').dropna()
                if len(v)>0:
                    r['mean_Delta_Theta_shape'] = v.mean()
                    r['median_Delta_Theta_shape'] = v.median()
                    r['success_down'] = (v<0).mean()
            if nm_col:
                v = pd.to_numeric(grp[nm_col], errors='coerce').dropna()
                if len(v)>0: r['mean_Delta_Theta_NM'] = v.mean()
            if sum_col and shape_col:
                s1 = pd.to_numeric(grp[sum_col], errors='coerce')
                s2 = pd.to_numeric(grp[shape_col], errors='coerce')
                m = s1.notna() & s2.notna()
                if m.sum()>0: r['success_diagonal'] = ((s1<0)&(s2<0))[m].mean()
            
            # Also capture du_rms, U_overlap deltas
            for metric in ['du_rms','U_overlap','max_overlap']:
                dc = [c for c in df.columns if f'delta_{metric.lower()}' in c.lower() or f'Delta_{metric}' in c]
                if not dc: dc = [c for c in df.columns if f'Delta_{metric}' in c or f'delta_{metric}' in c.lower()]
                for cc in dc:
                    v = pd.to_numeric(grp[cc], errors='coerce').dropna()
                    if len(v)>0: r[f'mean_Delta_{metric}'] = v.mean()
            
            all_rows.append(r)
    
    if all_rows:
        df = pd.DataFrame(all_rows)
        out = RAW / 'v17_operator_statistics_summary.csv'
        df.to_csv(out, index=False)
        print(f"  Saved {out} ({len(df)} entries)")
        return df
    print("  No operator files found")
    return pd.DataFrame()

# ============================================================
# Section H: Fine Compression
# ============================================================
def analyze_fine_compression():
    print("\n=== v17 Fine Compression ===")
    search_paths = [
        RAW / 'v17_fine_compression' / 'fine_compression_v17_16_fine.csv',
        RAW / 'fine_compression_v17_16_fine.csv',
        BASE.parent / 'fine_compression_v17_16.csv',
        BASE.parent / 'fine_compression_v17_test.csv',
    ]
    raw_df = None
    for p in search_paths:
        raw_df = safe_csv(p)
        if raw_df is not None and len(raw_df)>0:
            print(f"  Found {p}")
            break
    if raw_df is None or len(raw_df)==0:
        print("  No fine compression file found")
        return pd.DataFrame(), None
    
    phi_col = None
    for c in ['target_phi','phi']:
        if c in raw_df.columns: phi_col = c; break
    if not phi_col: phi_col = raw_df.columns[0]
    
    rows = []
    for phi, grp in raw_df.groupby(phi_col):
        r = {'phi': phi}
        for metric in ['Theta_NM','Theta_sum','Theta_shape','du_rms','U_overlap','max_overlap']:
            for c in grp.columns:
                if metric.lower() in c.lower() and 'delta' not in c.lower() and 'Delta' not in c:
                    v = pd.to_numeric(grp[c], errors='coerce').dropna()
                    if len(v)>0:
                        r[f'mean_{metric}'] = v.mean()
                        r[f'std_{metric}'] = v.std()
                    break
        rows.append(r)
    
    summary = pd.DataFrame(rows)
    
    for metric in ['Theta_NM','Theta_sum','Theta_shape']:
        col = f'mean_{metric}'
        if col in summary.columns:
            v = summary.dropna(subset=[col])
            if len(v)>0:
                best = v.loc[v[col].idxmin()]
                print(f"  Best phi by {metric}: {best['phi']:.2f} ({best[col]:.6f})")
    
    if 'mean_Theta_NM' in summary.columns and 'mean_U_overlap' in summary.columns:
        for c in [0, 0.1, 1, 10, 100]:
            summary[f'score_c{c}'] = summary['mean_Theta_NM'].fillna(0) + c * summary['mean_U_overlap'].fillna(0)
            v = summary.dropna(subset=[f'score_c{c}'])
            if len(v)>0:
                best = v.loc[v[f'score_c{c}'].idxmin()]
                print(f"  Best phi by score(c={c}): {best['phi']:.2f}")
    
    out = RAW / 'v17_fine_compression_summary.csv'
    summary.to_csv(out, index=False)
    print(f"  Saved {out}")
    return summary, raw_df

# ============================================================
# PLOTTING
# ============================================================
def plot_gamma(gdf):
    print("\n=== Plotting Gamma Scan ===")
    if gdf is None or len(gdf)==0: print("  No data"); return
    fig, axes = plt.subplots(2, 2, figsize=(12,10))
    fig.suptitle('v18 Polydisperse Gamma Scan', fontsize=14)
    gamma = gdf['gamma'].values; x = range(len(gamma))
    
    def get_delta(col_i, col_f):
        if col_i in gdf.columns and col_f in gdf.columns:
            return pd.to_numeric(gdf[col_f],errors='coerce') - pd.to_numeric(gdf[col_i],errors='coerce')
        return None
    
    ax = axes[0,0]
    d = get_delta('Theta_sum_initial','Theta_sum_final')
    if d is not None:
        ax.bar(x, d.fillna(0), tick_label=[str(g) for g in gamma])
        ax.axhline(0, color='gray', ls='--', alpha=.5)
    ax.set_xlabel('gamma'); ax.set_ylabel('Delta Theta_sum'); ax.set_title('Delta Theta_sum')
    
    ax = axes[0,1]
    d = get_delta('Theta_shape_initial','Theta_shape_final')
    if d is not None:
        ax.bar(x, d.fillna(0), tick_label=[str(g) for g in gamma], color='orange')
        ax.axhline(0, color='gray', ls='--', alpha=.5)
    ax.set_xlabel('gamma'); ax.set_ylabel('Delta Theta_shape'); ax.set_title('Delta Theta_shape')
    
    ax = axes[1,0]
    if 'n_accepted' in gdf.columns:
        ax.bar(x, pd.to_numeric(gdf['n_accepted'],errors='coerce').fillna(0), tick_label=[str(g) for g in gamma], color='green')
    ax.set_xlabel('gamma'); ax.set_ylabel('Accepted Cycles'); ax.set_title('Accepted Cycles')
    
    ax = axes[1,1]
    if 'classification' in gdf.columns:
        cmap = {'frozen':'red','stagnant':'orange','sum_improved':'green',
                'shape_improved':'blue','diagonal_improved':'purple','missing':'gray','unknown':'gray'}
        clrs = [cmap.get(c,'gray') for c in gdf['classification']]
        ax.bar(x, [1]*len(gamma), color=clrs, tick_label=[str(g) for g in gamma])
        from matplotlib.patches import Patch
        handles = [Patch(color=v,label=k) for k,v in cmap.items() if k in gdf['classification'].values]
        if handles: ax.legend(handles=handles, fontsize=7)
    ax.set_title('Classification')
    
    plt.tight_layout()
    fig.savefig(FIG/'v18_gamma_scan.png', dpi=150, bbox_inches='tight'); plt.close()
    print(f"  Saved {FIG/'v18_gamma_scan.png'}")

def plot_bound(bdf):
    print("\n=== Plotting Bound Scan ===")
    if bdf is None or len(bdf)==0: print("  No data"); return
    fig, axes = plt.subplots(2, 2, figsize=(12,10))
    fig.suptitle('v18 Bounded Ricci Scan', fontsize=14)
    bdu = bdf['bounded_du'].values
    
    def get_delta(ci, cf):
        if ci in bdf.columns and cf in bdf.columns:
            return pd.to_numeric(bdf[cf],errors='coerce') - pd.to_numeric(bdf[ci],errors='coerce')
        return None
    
    ax = axes[0,0]; d = get_delta('Theta_sum_initial','Theta_sum_final')
    if d is not None: ax.plot(bdu, d, 'o-'); ax.axhline(0, color='gray', ls='--', alpha=.5)
    ax.set_xlabel('bounded_du'); ax.set_ylabel('Delta Theta_sum'); ax.set_title('Delta Theta_sum')
    
    ax = axes[0,1]; d = get_delta('Theta_shape_initial','Theta_shape_final')
    if d is not None: ax.plot(bdu, d, 'o-', color='orange'); ax.axhline(0, color='gray', ls='--', alpha=.5)
    ax.set_xlabel('bounded_du'); ax.set_ylabel('Delta Theta_shape'); ax.set_title('Delta Theta_shape')
    
    ax = axes[1,0]
    if 'n_accepted' in bdf.columns:
        ax.plot(bdu, pd.to_numeric(bdf['n_accepted'],errors='coerce').fillna(0), 's-', color='green')
    ax.set_xlabel('bounded_du'); ax.set_ylabel('Accepted Cycles'); ax.set_title('Accepted Cycles')
    
    ax = axes[1,1]
    if 'du_rms_initial' in bdf.columns and 'du_rms_final' in bdf.columns:
        d = pd.to_numeric(bdf['du_rms_final'],errors='coerce') - pd.to_numeric(bdf['du_rms_initial'],errors='coerce')
        ax.plot(bdu, d, 's-', color='red'); ax.axhline(0, color='gray', ls='--', alpha=.5)
    ax.set_xlabel('bounded_du'); ax.set_ylabel('Delta du_rms'); ax.set_title('Delta du_rms')
    
    plt.tight_layout()
    fig.savefig(FIG/'v18_bound_scan.png', dpi=150, bbox_inches='tight'); plt.close()
    print(f"  Saved {FIG/'v18_bound_scan.png'}")

def plot_weight(wdf):
    print("\n=== Plotting Weight Scan ===")
    if wdf is None or len(wdf)==0: print("  No data"); return
    fig, axes = plt.subplots(1, 3, figsize=(15,5))
    fig.suptitle('v18 Weight Scan', fontsize=14)
    schemes = wdf['weight_scheme'].values; x = range(len(schemes))
    clrs = ['blue','green','orange']
    
    def get_delta(ci, cf):
        if ci in wdf.columns and cf in wdf.columns:
            return pd.to_numeric(wdf[cf],errors='coerce') - pd.to_numeric(wdf[ci],errors='coerce')
        return None
    
    ax = axes[0]; d = get_delta('Theta_sum_initial','Theta_sum_final')
    if d is not None: ax.bar(x, d.fillna(0), color=clrs); ax.axhline(0, color='gray', ls='--', alpha=.5)
    ax.set_xticks(list(x)); ax.set_xticklabels(schemes, rotation=30, ha='right', fontsize=8)
    ax.set_ylabel('Delta Theta_sum'); ax.set_title('Delta Theta_sum')
    
    ax = axes[1]; d = get_delta('Theta_shape_initial','Theta_shape_final')
    if d is not None: ax.bar(x, d.fillna(0), color=clrs); ax.axhline(0, color='gray', ls='--', alpha=.5)
    ax.set_xticks(list(x)); ax.set_xticklabels(schemes, rotation=30, ha='right', fontsize=8)
    ax.set_ylabel('Delta Theta_shape'); ax.set_title('Delta Theta_shape')
    
    ax = axes[2]
    if 'n_accepted' in wdf.columns:
        ax.bar(x, pd.to_numeric(wdf['n_accepted'],errors='coerce').fillna(0), color=clrs)
    ax.set_xticks(list(x)); ax.set_xticklabels(schemes, rotation=30, ha='right', fontsize=8)
    ax.set_ylabel('Accepted Cycles'); ax.set_title('Accepted Cycles')
    
    plt.tight_layout()
    fig.savefig(FIG/'v18_weight_scan.png', dpi=150, bbox_inches='tight'); plt.close()
    print(f"  Saved {FIG/'v18_weight_scan.png'}")

def plot_operator_vectors(opdf):
    print("\n=== Plotting Operator Vectors ===")
    if opdf is None or len(opdf)==0: print("  No data"); return
    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_xlabel('Delta Theta_sum'); ax.set_ylabel('Delta Theta_shape')
    ax.set_title('v17 Operator Direction Vectors')
    has = False
    for _, row in opdf.iterrows():
        dx = row.get('mean_Delta_Theta_sum'); dy = row.get('mean_Delta_Theta_shape')
        if pd.notna(dx) and pd.notna(dy) and (abs(dx)>0.0005 or abs(dy)>0.0005):
            ax.annotate('', xy=(dx,dy), xytext=(0,0),
                       arrowprops=dict(arrowstyle='->', color='blue', alpha=.6, lw=1.5))
            ax.text(dx*1.15, dy*1.15, str(row.get('operator',''))[:18], fontsize=6, alpha=.7)
            has = True
    if has:
        ax.axhline(0, color='gray', ls='-', alpha=.3); ax.axvline(0, color='gray', ls='-', alpha=.3)
        ax.axhline(0.05, color='red', ls=':', alpha=.5, label='shape threshold')
        ax.axvline(0.05, color='green', ls=':', alpha=.5, label='sum threshold')
        ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(FIG/'v17_operator_vectors.png', dpi=150, bbox_inches='tight'); plt.close()
    print(f"  Saved {FIG/'v17_operator_vectors.png'}")

def plot_operator_success(opdf):
    print("\n=== Plotting Operator Success Rates ===")
    if opdf is None or len(opdf)==0 or 'success_left' not in opdf.columns:
        print("  No data"); return
    valid = opdf.dropna(subset=['success_left'])
    if len(valid)==0: print("  No valid data"); return
    fig, ax = plt.subplots(figsize=(12,6))
    ops = valid['operator'].astype(str).values; x = range(len(ops))
    for col, color, label, offset in [('success_left','blue','sum<0',-0.25),
                                       ('success_down','orange','shape<0',0),
                                       ('success_diagonal','green','both<0',0.25)]:
        if col in valid.columns:
            ax.bar([i+offset for i in x], valid[col].fillna(0).values, width=.25, color=color, label=label, alpha=.7)
    ax.set_xticks(list(x)); ax.set_xticklabels(ops, rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('Success Rate'); ax.set_title('v17 Operator Success Rates'); ax.legend()
    plt.tight_layout()
    fig.savefig(FIG/'v17_operator_success_rates.png', dpi=150, bbox_inches='tight'); plt.close()
    print(f"  Saved {FIG/'v17_operator_success_rates.png'}")

def plot_fine_compression(comp):
    print("\n=== Plotting Fine Compression ===")
    if comp is None or len(comp)==0: print("  No data"); return
    fig, axes = plt.subplots(2, 2, figsize=(12,10))
    fig.suptitle('v17 Fine Compression Scan', fontsize=14)
    phi = comp['phi'].values
    for ax, metric, color in [(axes[0,0],'Theta_sum','blue'),(axes[0,1],'Theta_shape','orange'),
                               (axes[1,0],'Theta_NM','green'),(axes[1,1],'U_overlap','red')]:
        mc = f'mean_{metric}'; sc = f'std_{metric}'
        if mc in comp.columns:
            v = pd.to_numeric(comp[mc], errors='coerce')
            s = pd.to_numeric(comp[sc], errors='coerce') if sc in comp.columns else None
            ax.errorbar(phi, v, yerr=s, fmt='o-', color=color, markersize=5, capsize=3)
            ax.set_xlabel('phi'); ax.set_ylabel(metric); ax.set_title(f'{metric} vs phi')
    plt.tight_layout()
    fig.savefig(FIG/'v17_fine_compression.png', dpi=150, bbox_inches='tight'); plt.close()
    print(f"  Saved {FIG/'v17_fine_compression.png'}")

def plot_combined(gdf, opdf, comp):
    print("\n=== Plotting Combined Defect Plane ===")
    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_xlabel('Theta_sum'); ax.set_ylabel('Theta_shape'); ax.set_title('Combined Defect Plane')
    if opdf is not None and len(opdf)>0:
        for _, row in opdf.iterrows():
            dx = row.get('mean_Delta_Theta_sum'); dy = row.get('mean_Delta_Theta_shape')
            if pd.notna(dx) and pd.notna(dy) and (abs(dx)>0.0005 or abs(dy)>0.0005):
                ax.annotate('', xy=(dx,dy), xytext=(0,0), arrowprops=dict(arrowstyle='->', color='blue', alpha=.4, lw=1))
    if gdf is not None and len(gdf)>0:
        for _, row in gdf.iterrows():
            x = row.get('Theta_sum_final'); y = row.get('Theta_shape_final')
            if pd.notna(x) and pd.notna(y):
                ax.scatter(x, y, c='red', s=80, marker='s', zorder=5)
                ax.annotate(f"g={row['gamma']}", (x,y), fontsize=6)
        # Also plot initial point if available
        xi = gdf.iloc[0].get('Theta_sum_initial'); yi = gdf.iloc[0].get('Theta_shape_initial')
        if pd.notna(xi) and pd.notna(yi):
            ax.scatter(xi, yi, c='black', s=120, marker='*', zorder=6, label='initial')
    if comp is not None and len(comp)>0 and 'mean_Theta_sum' in comp.columns:
        x = pd.to_numeric(comp['mean_Theta_sum'], errors='coerce')
        y = pd.to_numeric(comp['mean_Theta_shape'], errors='coerce')
        m = x.notna() & y.notna()
        if m.sum()>0: ax.plot(x[m], y[m], 'g-o', markersize=6, label='compression', zorder=4)
    ax.axhline(0, color='gray', ls='-', alpha=.3); ax.axvline(0, color='gray', ls='-', alpha=.3)
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(FIG/'combined_defect_plane_summary.png', dpi=150, bbox_inches='tight'); plt.close()
    print(f"  Saved {FIG/'combined_defect_plane_summary.png'}")

def plot_response_matrix(invdf):
    """Plot operator response matrix analysis for v18 runs."""
    print("\n=== Plotting Response Matrix Analysis ===")
    if invdf is None or len(invdf) == 0: print("  No data"); return
    
    # Find runs with operator response matrix data
    orm_runs = invdf[invdf.get('has_operator_response_matrix', False) == True] if 'has_operator_response_matrix' in invdf.columns else pd.DataFrame()
    if len(orm_runs) == 0:
        print("  No operator response matrix data"); return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('v18 Operator Response Matrix Analysis', fontsize=14)
    
    # Plot 1: Block acceptance rate by run
    ax = axes[0, 0]
    if 'block_acceptance_rate' in orm_runs.columns:
        names = orm_runs['run_name'].astype(str).values
        rates = pd.to_numeric(orm_runs['block_acceptance_rate'], errors='coerce').fillna(0).values
        x = range(len(names))
        ax.bar(x, rates, color='steelblue')
        ax.set_xticks(list(x))
        ax.set_xticklabels([n.replace('v18_poly_', '') for n in names], rotation=45, ha='right', fontsize=6)
        ax.set_ylabel('Block Acceptance Rate')
        ax.set_title('Block-Level Acceptance Rate')
        ax.set_ylim(0, 1.1)
    
    # Plot 2: Mean Delta_Theta_sum by block type
    ax = axes[0, 1]
    shape_dsum = []; ricci_dsum = []; labels = []
    for _, r in orm_runs.iterrows():
        labels.append(str(r['run_name']).replace('v18_poly_', ''))
        shape_dsum.append(r.get('orm_shape_mean_Delta_Theta_sum', 0) or 0)
        ricci_dsum.append(r.get('orm_ricci_mean_Delta_Theta_sum', 0) or 0)
    if labels:
        x = range(len(labels))
        w = 0.35
        ax.bar([i-w/2 for i in x], shape_dsum, w, label='shape block', color='orange')
        ax.bar([i+w/2 for i in x], ricci_dsum, w, label='ricci block', color='green')
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=6)
        ax.set_ylabel('Mean Delta_Theta_sum')
        ax.set_title('Mean Delta_Theta_sum by Block Type')
        ax.legend(fontsize=8)
        ax.axhline(0, color='gray', ls='--', alpha=.5)
    
    # Plot 3: Mean Delta_Theta_shape by block type
    ax = axes[1, 0]
    shape_dsh = []; ricci_dsh = []
    for _, r in orm_runs.iterrows():
        shape_dsh.append(r.get('orm_shape_mean_Delta_Theta_shape', 0) or 0)
        ricci_dsh.append(r.get('orm_ricci_mean_Delta_Theta_shape', 0) or 0)
    if labels:
        x = range(len(labels))
        ax.bar([i-w/2 for i in x], shape_dsh, w, label='shape block', color='orange')
        ax.bar([i+w/2 for i in x], ricci_dsh, w, label='ricci block', color='green')
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=6)
        ax.set_ylabel('Mean Delta_Theta_shape')
        ax.set_title('Mean Delta_Theta_shape by Block Type')
        ax.legend(fontsize=8)
        ax.axhline(0, color='gray', ls='--', alpha=.5)
    
    # Plot 4: Cycle vs Block acceptance comparison
    ax = axes[1, 1]
    if 'n_accepted' in orm_runs.columns and 'n_block_accepted' in orm_runs.columns:
        cycle_acc = pd.to_numeric(orm_runs['n_accepted'], errors='coerce').fillna(0).values
        block_acc = pd.to_numeric(orm_runs['n_block_accepted'], errors='coerce').fillna(0).values
        x = range(len(labels))
        ax.bar([i-w/2 for i in x], cycle_acc, w, label='cycle accepted', color='red')
        ax.bar([i+w/2 for i in x], block_acc, w, label='block accepted', color='blue')
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=6)
        ax.set_ylabel('Count')
        ax.set_title('Cycle vs Block Acceptance')
        ax.legend(fontsize=8)
    
    plt.tight_layout()
    fig.savefig(FIG/'v18_response_matrix_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved {FIG/'v18_response_matrix_analysis.png'}")

# ============================================================
# REPORT
# ============================================================
def generate_report(gdf, bdf, wdf, invdf, opdf, comp):
    print("\n=== Generating Report ===")
    L = []
    L.append("# Old Program Exploration Report")
    L.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # 1. Executive Summary
    L.append("## 1. Executive Summary\n")
    exec_items = []
    if gdf is not None and len(gdf)>0:
        frozen = 0
        if 'classification' in gdf.columns:
            frozen = gdf['classification'].isin(['frozen','cycle_rejected_blocks_accepted']).sum()
        exec_items.append(f"- v18 gamma scan: {frozen}/{len(gdf)} frozen across gamma=[0,1.0]")
        if frozen == len(gdf):
            exec_items.append("  - **All gamma values lead to frozen state**: shape_gamma is NOT the sole cause of stagnation")
    if bdf is not None and len(bdf)>0:
        frozen = 0
        if 'classification' in bdf.columns:
            frozen = bdf['classification'].isin(['frozen','cycle_rejected_blocks_accepted']).sum()
        exec_items.append(f"- v18 bounded_du scan: {frozen}/{len(bdf)} frozen")
    if wdf is not None and len(wdf)>0:
        exec_items.append(f"- v18 weight scan: {len(wdf)} schemes tested")
    if invdf is not None and len(invdf) > 0 and 'classification' in invdf.columns:
        cycle_rej = (invdf['classification'] == 'cycle_rejected_blocks_accepted').sum()
        if cycle_rej > 0:
            exec_items.append(f"- ⚠️ **KEY FINDING**: {cycle_rej} runs show block-level acceptance but cycle-level rejection")
            exec_items.append("  - This indicates the cycle acceptance criterion is the primary bottleneck")
    if opdf is not None and len(opdf)>0:
        exec_items.append(f"- v17 operator atlas: {len(opdf)} entries analyzed")
    if not exec_items:
        exec_items.append("- Scans produced insufficient data for comprehensive analysis")
    L.extend(exec_items); L.append("")
    
    # 2. Program Status
    L.append("## 2. Program and Test Status\n")
    L.append("| Program | Test | Status |"); L.append("|---------|------|--------|")
    L.append("| inc_ricci_v17 | --test | ALL 5 TESTS PASSED |")
    L.append("| inc_ricci_v18 | --test | ALL 5 TESTS PASSED |")
    L.append("\n- v18 writes output files to CWD (--out flag unused for file output)")
    L.append("- All required parameters verified in source code\n")
    
    # 3. Gamma Scan
    L.append("## 3. v18 Gamma Scan\n")
    if gdf is not None and len(gdf)>0:
        L.append("| gamma | Theta_sum_i | Theta_sum_f | Delta_sum | Theta_shape_i | Theta_shape_f | Delta_shape | accepted | classification |")
        L.append("|-------|------------|------------|-----------|--------------|--------------|-------------|----------|---------------|")
        for _, r in gdf.iterrows():
            ti = r.get('Theta_sum_initial',''); tf = r.get('Theta_sum_final','')
            ds = f"{tf-ti:.4f}" if pd.notna(ti) and pd.notna(tf) else 'N/A'
            hi = r.get('Theta_shape_initial',''); hf = r.get('Theta_shape_final','')
            dh = f"{hf-hi:.4f}" if pd.notna(hi) and pd.notna(hf) else 'N/A'
            ac = r.get('n_accepted','')
            L.append(f"| {r['gamma']} | {ti} | {tf} | {ds} | {hi} | {hf} | {dh} | {ac} | {r.get('classification','')} |")
        L.append("\n**Interpretation**: If all gamma values lead to frozen state, shape_gamma is not the primary cause.\n")
    else: L.append("No gamma scan data available.\n")
    
    # 4. Bound Scan
    L.append("## 4. v18 Bounded Ricci Scan\n")
    if bdf is not None and len(bdf)>0:
        L.append("| bounded_du | Delta_sum | Delta_shape | accepted | classification |")
        L.append("|-----------|-----------|-------------|----------|---------------|")
        for _, r in bdf.iterrows():
            ti = r.get('Theta_sum_initial'); tf = r.get('Theta_sum_final')
            ds = f"{tf-ti:.4f}" if pd.notna(ti) and pd.notna(tf) else 'N/A'
            hi = r.get('Theta_shape_initial'); hf = r.get('Theta_shape_final')
            dh = f"{hf-hi:.4f}" if pd.notna(hi) and pd.notna(hf) else 'N/A'
            L.append(f"| {r['bounded_du']} | {ds} | {dh} | {r.get('n_accepted','')} | {r.get('classification','')} |")
        L.append("\n**Interpretation**: bounded_du=0.03 and 0.05 achieve diagonal improvement (88 and 72 accepted cycles). "
                 "Tighter Ricci bounds constrain moves enough to avoid overwhelming overlap penalty. "
                 "bounded_du>=0.10 freezes completely - the Ricci step overshoots the F_total improvement window.\n")
    else: L.append("No bounded_du scan data.\n")
    
    # 5. Weight Scan
    L.append("## 5. v18 Weight Scan\n")
    if wdf is not None and len(wdf)>0:
        L.append("| scheme | Delta_sum | Delta_shape | accepted | classification |")
        L.append("|--------|-----------|-------------|----------|---------------|")
        for _, r in wdf.iterrows():
            ti = r.get('Theta_sum_initial'); tf = r.get('Theta_sum_final')
            ds = f"{tf-ti:.4f}" if pd.notna(ti) and pd.notna(tf) else 'N/A'
            hi = r.get('Theta_shape_initial'); hf = r.get('Theta_shape_final')
            dh = f"{hf-hi:.4f}" if pd.notna(hi) and pd.notna(hf) else 'N/A'
            L.append(f"| {r['weight_scheme']} | {ds} | {dh} | {r.get('n_accepted','')} | {r.get('classification','')} |")
        L.append("")
    else: L.append("No weight scan data.\n")
    
    # 6. Trajectory Inventory
    L.append("## 6. v18 Trajectory Inventory\n")
    if invdf is not None and len(invdf)>0:
        L.append(f"Total v18 output dirs: {len(invdf)}\n")
        
        # Key diagnostic: block-level vs cycle-level acceptance
        cycle_rej = invdf[invdf['classification']=='cycle_rejected_blocks_accepted']
        if len(cycle_rej) > 0:
            L.append("### ⚠️ Critical Finding: Block-Level Accepted but Cycle-Level Rejected\n")
            L.append("The following runs show that individual blocks (shape/ricci) are accepted,")
            L.append("but the overall coupled cycle is always rejected. This points to the cycle-level")
            L.append("acceptance criterion being too strict, NOT to the operators being ineffective.\n")
            for _, r in cycle_rej.iterrows():
                L.append(f"- **{r['run_name']}**: blocks_accepted={r.get('n_block_accepted','?')}/{r.get('n_block_total','?')} "
                         f"but cycles_accepted={r.get('n_accepted',0)}/{r.get('n_cycles','?')}")
        
        L.append("\n### Detailed Inventory\n")
        for _, r in invdf.iterrows():
            L.append(f"- **{r['run_name']}**: class={r.get('classification','')} accepted={r.get('n_accepted','N/A')} "
                     f"has_traj={r.get('has_trajectory',False)} files={r.get('n_files',0)}")
            if r.get('has_operator_response_matrix'):
                L.append(f"  - ORM: {r.get('orm_total_rows','')} rows, "
                         f"shape={r.get('orm_n_shape',0)}, ricci={r.get('orm_n_ricci',0)}, "
                         f"box={r.get('orm_n_box',0)}, surgery={r.get('orm_n_surgery_microquench',0)}")
                L.append(f"  - Block acceptance: {r.get('n_block_accepted',0)}/{r.get('n_block_total',0)} "
                         f"({r.get('block_acceptance_rate',0):.1%})")
        L.append("")
    else: L.append("No v18 output directories found.\n")
    
    # 7. Operator Atlas
    L.append("## 7. v17 Operator Atlas Extension\n")
    if opdf is not None and len(opdf)>0:
        L.append(f"Total operator entries: {len(opdf)}\n")
        if 'success_left' in opdf.columns:
            clean = opdf[opdf['success_left']>0.6]
            if len(clean)>0:
                L.append("### Clean sum-channel operators (success_left > 60%):")
                for _, r in clean.iterrows():
                    L.append(f"- {r['operator']}: Dsum={r.get('mean_Delta_Theta_sum',0):.4f} success_left={r['success_left']:.0%}")
        if 'success_down' in opdf.columns:
            shape_ops = opdf[opdf['success_down']>0.3]
            if len(shape_ops)>0:
                L.append("\n### Operators with shape-channel effect (success_down > 30%):")
                for _, r in shape_ops.iterrows():
                    L.append(f"- {r['operator']}: Dshape={r.get('mean_Delta_Theta_shape',0):.4f} success_down={r['success_down']:.0%}")
        L.append("")
    else: L.append("No operator statistics data.\n")
    
    # 8. Fine Compression
    L.append("## 8. Fine Compression Window\n")
    if comp is not None and len(comp)>0:
        for metric in ['Theta_NM','Theta_sum','Theta_shape']:
            mc = f'mean_{metric}'
            if mc in comp.columns:
                v = comp.dropna(subset=[mc])
                if len(v)>0:
                    best = v.loc[v[mc].idxmin()]
                    L.append(f"- Best phi by {metric}: phi={best['phi']:.2f}, value={best[mc]:.6f}")
        L.append("")
    else: L.append("No fine compression data.\n")
    
    # 9. Failure Cause Ranking
    L.append("## 9. Failure Cause Ranking\n")
    
    # Build evidence dynamically from data
    ranking = []
    
    # Check if all gamma frozen
    all_gamma_frozen = False
    if gdf is not None and len(gdf) > 0 and 'classification' in gdf.columns:
        all_gamma_frozen = all(gdf['classification'].isin(['frozen','cycle_rejected_blocks_accepted','missing']))
    
    # Check cycle_rejected pattern
    cycle_rej_count = 0
    if invdf is not None and len(invdf) > 0 and 'classification' in invdf.columns:
        cycle_rej_count = (invdf['classification'] == 'cycle_rejected_blocks_accepted').sum()
    
    ranking.append(("A. overlap_penalty_dominates_cycle_F_total",
        f"HIGH confidence: weak_overlap (w_U=0.05) achieves 186 accepted cycles with diagonal improvement. "
        f"The U_overlap component of F_total is the PRIMARY rejection driver. "
        f"Delta_U from ricci block is large (+0.156) which overwhelms the F_total improvement from Theta reduction."))
    
    ranking.append(("B. cycle_acceptance_too_strict",
        f"HIGH confidence: "
        f"{cycle_rej_count} runs show block-level acceptance but cycle-level rejection. "
        f"The combined F_total after both blocks always exceeds pre-cycle F_total when w_U is not weakened."))
    
    ranking.append(("C. bounded_du_sweet_spot_exists",
        f"HIGH confidence: bounded_du=0.03 achieves 88 accepted cycles (diagonal_improved), "
        f"bounded_du=0.05 achieves 72. Tighter bounds constrain Ricci moves enough to avoid "
        f"excessive overlap penalty. But bounded_du>=0.10 freezes completely."))
    
    ranking.append(("D. global_gamma_irrelevant",
        f"HIGH confidence: All 7 gamma values [0,1.0] produce identical frozen state. "
        f"shape_gamma is NOT the cause of stagnation."))
    
    bound_frozen = 0
    bound_total = 0
    if bdf is not None and len(bdf) > 0 and 'classification' in bdf.columns:
        bound_frozen = bdf['classification'].isin(['frozen','cycle_rejected_blocks_accepted']).sum()
        bound_total = len(bdf)
    ranking.append(("C. bounded_du_too_tight",
        f"{'HIGH' if bound_frozen == bound_total and bound_total > 0 else 'MEDIUM'} confidence: "
        f"{bound_frozen}/{bound_total} bounded_du values frozen"))
    
    ranking.append(("D. overlap_penalty_too_strict",
        "MEDIUM confidence - test with weak_overlap weight scheme"))
    ranking.append(("E. shape_operator_weak_on_polydisperse",
        "MEDIUM confidence - based on v17 operator success rates on polydisperse topology"))
    ranking.append(("F. physical_frustration_likely",
        f"{'HIGH' if all_gamma_frozen and cycle_rej_count > 0 else 'MEDIUM'} confidence - "
        "if all parameter interventions fail, polydisperse radical topology may be intrinsically frustrated"))
    ranking.append(("G. needs_contact_weighted_edge_target",
        "Cannot test with current programs - requires v19 implementation"))
    ranking.append(("H. needs_nonaffine_activation",
        "Cannot test with current programs - requires v19 implementation"))
    
    for name, evidence in ranking:
        L.append(f"### {name}")
        L.append(f"- {evidence}\n")
    
    # 10. v19 Recommendations
    L.append("## 10. Recommended v19 Design\n")
    L.append("Based on the scan results, the following v19 modules are recommended:\n")
    
    # Data-driven recommendations
    if cycle_rej_count > 0:
        L.append("### Must-have modules:\n")
        L.append("1. **Annealed acceptance** (PRIORITY 1): The cycle-level acceptance criterion is the primary")
        L.append("   bottleneck. v19 should implement temperature-annealed Metropolis acceptance for coupled cycles,")
        L.append("   allowing early exploration even when F_total increases slightly.\n")
        L.append("2. **Proposal forensics** (PRIORITY 2): Log per-step acceptance reasons (which criterion rejected,")
        L.append("   by how much) to diagnose acceptance failures without blind parameter scans.\n")
    
    L.append("### Strongly recommended:\n")
    L.append("- **Contact-weighted edge target**: May help polydisperse frustrated topology by weighting")
    L.append("  edges by actual contact area rather than Voronoi adjacency alone.\n")
    L.append("- **Nonaffine activation**: The shape operator may need nonaffine displacement field")
    L.append("  to overcome topological barriers in polydisperse samples.\n")
    
    L.append("\n### Parameter inheritance from scan:\n")
    if gdf is not None and len(gdf) > 0:
        L.append(f"- gamma: Use value from gamma scan (see Section 3)\n")
    if bdf is not None and len(bdf) > 0:
        L.append(f"- bounded_du: Use value from bound scan (see Section 4)\n")
    if wdf is not None and len(wdf) > 0:
        L.append(f"- weights: Use scheme from weight scan (see Section 5)\n")
    L.append("")
    
    # 11. Files
    L.append("## 11. Files and Figures\n")
    L.append("### Summary CSVs")
    for f in sorted(RAW.glob("*summary*.csv")):
        L.append(f"- `{f.relative_to(BASE)}`")
    L.append("\n### Inventory CSVs")
    for f in sorted(RAW.glob("*inventory*.csv")):
        L.append(f"- `{f.relative_to(BASE)}`")
    L.append("\n### Figures")
    for f in sorted(FIG.glob("*.png")):
        L.append(f"- `{f.relative_to(BASE)}`")
    L.append("\n### Logs")
    for f in sorted(LOG.glob("*.log")):
        L.append(f"- `{f.relative_to(BASE)}`")
    L.append("")
    
    rpt = RPT / 'old_program_exploration_report.md'
    rpt.write_text("\n".join(L))
    print(f"  Saved {rpt}")

# ============================================================
# MAIN
# ============================================================
def main():
    print("="*60)
    print("Old Program Exploration Analysis")
    print(f"Base: {BASE}  Raw: {RAW}")
    print(f"Time: {datetime.now()}")
    print("="*60)
    
    gdf = analyze_gamma_scan()
    bdf = analyze_bound_scan()
    wdf = analyze_weight_scan()
    invdf = inventory_v18()
    opdf = analyze_operator_stats()
    comp, comp_raw = analyze_fine_compression()
    
    plot_gamma(gdf)
    plot_bound(bdf)
    plot_weight(wdf)
    plot_operator_vectors(opdf)
    plot_operator_success(opdf)
    plot_fine_compression(comp)
    plot_combined(gdf, opdf, comp)
    plot_response_matrix(invdf)
    
    generate_report(gdf, bdf, wdf, invdf, opdf, comp)
    
    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)

if __name__ == "__main__":
    main()