#!/usr/bin/env python3
"""
CuBind.py — Place Cu(I)/Cu(II) at His-N sites in any PDB.
Usage:  python3 CuBind.py input.pdb -1|-2 [options]
        -1 Cu(I) 2×His  |  -2 Cu(II) 3×His
Options: --cu-n-ideal Å[2.00] --cu-n-max Å[2.60] --don-cutoff Å[7.00]
         --cu-sep Å[5.00] --clash Å[1.80]
Output:  <stem>_1Cu.pdb  or  <stem>_2Cu.pdb
"""
import sys, math, argparse, textwrap
from itertools import combinations
import numpy as np
from scipy.optimize import minimize

def parse_pdb(path):
    atoms, atom_lines, header = [], [], []
    try:    lines = open(path).readlines()
    except OSError as e: sys.exit(f"[ERROR] {e}")
    for line in lines:
        rec = line[:6].strip()
        if rec in ('ATOM','HETATM'):
            try:
                atoms.append({'record':rec, 'atname':line[12:16],
                    'resname':line[17:20].strip(),
                    'chain':line[21] if len(line)>21 else ' ',
                    'resnum':int(line[22:26]),
                    'x':float(line[30:38]),'y':float(line[38:46]),'z':float(line[46:54])})
            except (ValueError,IndexError): pass
            atom_lines.append(line)
        elif rec not in ('TER','END','ENDMDL'):
            header.append(line)
    if not atoms: sys.exit(f"[ERROR] No ATOM/HETATM records in {path}")
    return atoms, atom_lines, header

def cu_line(serial, resnum, xyz):
    x,y,z = xyz
    return (f"HETATM{serial:5d}  CU  CU  {resnum:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          Cu\n")

def renumber(lines):
    n, out = 1, []
    for l in lines:
        if l[:6].strip() in ('ATOM','HETATM'):
            l = l[:6]+f"{n:5d}"+l[11:]; n+=1
        out.append(l)
    return out

def write_pdb(path, header, atom_lines, cu_records, sites, ox):
    ox_str = 'I' if ox==1 else 'II'
    rem = ([f"REMARK CuBind_2  Cu({ox_str})  sites={len(sites)}\n",
             "REMARK Donors: HIE=ND1, HID=NE2\n"] +
           [f"REMARK Cu{i:02d}: "
            +"|".join(f"{d['resname']}{d['resnum']}({d['chain']})-{d['atname'].strip()}"
                      for d in s['donors'])
            +"  Cu-N="+"/".join(f"{v:.3f}" for v in s['dists'])
            +f"A  xyz=({s['xyz'][0]:.2f},{s['xyz'][1]:.2f},{s['xyz'][2]:.2f})\n"
            for i,s in enumerate(sites,1)])
    open(path,'w').writelines(
        renumber(header+rem+atom_lines+cu_records+['TER\n','END\n']))

def d3(a,b): return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

def best_cu(pts, ideal=2.00):
    pts = [np.asarray(p) for p in pts]
    r = minimize(lambda p: sum((np.linalg.norm(p-c)-ideal)**2 for c in pts),
                 np.mean(pts,axis=0), method='Nelder-Mead',
                 options={'xatol':1e-6,'fatol':1e-6,'maxiter':100_000})
    return tuple(float(v) for v in r.x), float(r.fun)

def get_donors(atoms):
    return [a for a in atoms
            if (a['resname']=='HIE' and a['atname'].strip()=='ND1') or
               (a['resname']=='HID' and a['atname'].strip()=='NE2')]

def find_sites(donors, atoms, cfg, cn):
    coords  = [np.array((d['x'],d['y'],d['z'])) for d in donors]
    all_xyz = [(a['x'],a['y'],a['z']) for a in atoms]
    dkey    = lambda d: (d['chain'],d['resnum'])
    sites   = []
    for combo in combinations(range(len(donors)),cn):
        pts = [coords[i] for i in combo]
        if any(np.linalg.norm(pts[a]-pts[b])>cfg.don_cutoff
               for a in range(cn) for b in range(a+1,cn)): continue
        cu,score = best_cu(pts,cfg.cu_n_ideal)
        dists = [float(np.linalg.norm(np.asarray(cu)-p)) for p in pts]
        if max(dists)>cfg.cu_n_max: continue
        exempt = {dkey(donors[i]) for i in combo}
        if any(d3(cu,xyz)<cfg.clash for a,xyz in zip(atoms,all_xyz)
               if dkey(a) not in exempt): continue
        sites.append({'donors':[donors[i] for i in combo],
                      'xyz':cu,'dists':dists,'score':score,'cn':cn})
    sites.sort(key=lambda s:s['score'])
    used_d,used_cu,final = set(),[],[]
    for s in sites:
        keys = tuple(dkey(d) for d in s['donors'])
        if (any(d3(s['xyz'],uc)<cfg.cu_sep for uc in used_cu) or
                any(k in used_d for k in keys)): continue
        final.append(s); used_cu.append(s['xyz']); used_d.update(keys)
    return final

def print_summary(inp, out, atoms, donors, hip, sites, ox):
    ox_str = 'I' if ox==1 else 'II'
    hie_n  = sum(1 for d in donors if d['resname']=='HIE')
    hid_n  = sum(1 for d in donors if d['resname']=='HID')

    # fixed column widths  Rank CuID Donors  Dist  Score
    C1,C2,C3,C4,C5 = 6, 8, 40, 22, 10

    def fmt_donors(s):
        # fixed per-donor width so columns stay rigid regardless of resnum digits
        parts = [f"{d['chain']}:{d['resname']}{d['resnum']}" for d in s['donors']]
        return "  ".join(f"{p:<10}" for p in parts)

    def fmt_dists(s):
        return "  ".join(f"{v:.2f}" for v in s['dists'])

    def fmt_score(v):
        if v >= 0.001: return f"{v:.3f}"
        m,e = f"{v:.3e}".split('e')
        return f"{float(m):.3f}e{int(e):+d}"

    SEP = '-' * (C1+C2+C3+C4+C5+8)

    print(f"\n  {'CuBind_2':<20}  Cu({ox_str}) Binding Site Identification")
    print(f"  {'Input':<20}  {inp}")
    print(f"  {'Output':<20}  {out}\n")
    print(f"  {'Total atoms':<20}  {len(atoms)}")
    print(f"  {'Histidine donors':<20}  {len(donors)}")
    print(f"  {'  HIE':<20}  {hie_n}")
    print(f"  {'  HID':<20}  {hid_n}")
    print(f"  {'  HIP (skipped)':<20}  {hip}")
    print(f"  {'  Cu sites found':<20}  {len(sites)}")
    print(f"\n  Top 5 Ranked Binding Sites")
    print(f"  {SEP}")
    print(f"  {'Rank':<{C1}} {'Cu ID':<{C2}} {'Donor Residues':<{C3}} {'Cu-N Distances (A)':<{C4}} {'Score':<{C5}}")
    print(f"  {SEP}")
    for i,s in enumerate(sites[:5],1):
        print(f"  {i:<{C1}} {s['resnum']:<{C2}} {fmt_donors(s):<{C3}} {fmt_dists(s):<{C4}} {fmt_score(s['score']):<{C5}}")
    print(f"  {SEP}\n")

def main():
    p = argparse.ArgumentParser(prog='CuBind_2.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            Place Cu(I)/Cu(II) at His-N sites in any PDB.
              python3 CuBind_2.py protein.pdb -1   # Cu(I)  2×His
              python3 CuBind_2.py protein.pdb -2   # Cu(II) 3×His"""))
    p.add_argument('input')
    p.add_argument('-1',dest='ox1',action='store_true')
    p.add_argument('-2',dest='ox2',action='store_true')
    p.add_argument('--cu-n-ideal',type=float,default=2.00,metavar='Å')
    p.add_argument('--cu-n-max',  type=float,default=2.60,metavar='Å')
    p.add_argument('--don-cutoff',type=float,default=7.00,metavar='Å')
    p.add_argument('--cu-sep',    type=float,default=5.00,metavar='Å')
    p.add_argument('--clash',     type=float,default=1.80,metavar='Å')
    args = p.parse_args()

    if not (args.ox1^args.ox2): p.error("Specify exactly one of -1 (CuI) or -2 (CuII)")
    ox,cn = (1,2) if args.ox1 else (2,3)
    cfg   = argparse.Namespace(cu_n_ideal=args.cu_n_ideal, cu_n_max=args.cu_n_max,
                               don_cutoff=args.don_cutoff, cu_sep=args.cu_sep, clash=args.clash)
    stem  = args.input.removesuffix('.pdb') if args.input.lower().endswith('.pdb') else args.input

    atoms,atom_lines,header = parse_pdb(args.input)
    donors = get_donors(atoms)
    hip    = len({(a['chain'],a['resnum']) for a in atoms if a['resname']=='HIP'})

    if len(donors)<cn: sys.exit(f"[ERROR] Need ≥{cn} His-N donors, found {len(donors)}")
    sites = find_sites(donors,atoms,cfg,cn)
    if not sites: sys.exit("  [!] No valid sites found. Try --cu-n-max 2.8 or --don-cutoff 8.0")

    rn_max = max((a['resnum'] for a in atoms if a['record']=='ATOM'),default=0)
    for i,s in enumerate(sites,1): s['resnum'] = rn_max+i

    out = f"{stem}_{ox}Cu.pdb"
    write_pdb(out,header,atom_lines,
              [cu_line(0,s['resnum'],s['xyz']) for s in sites],sites,ox)
    print_summary(args.input,out,atoms,donors,hip,sites,ox)

if __name__=='__main__': main()