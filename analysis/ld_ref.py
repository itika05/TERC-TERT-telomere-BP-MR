"""1000 Genomes phase 3 EUR (503 individuals, 1006 phased haplotypes) reference panel utilities.
Genotypes were read directly from the phase 3 release VCFs (ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/, GRCh37) with tabix
range requests; each record is stored as [chr, pos, id, ref, alt, alt_count, n_missing, base64 haplotype bits, multiallelic flag, VT].
LD r is computed from phased haplotypes, so it is signed; signs are aligned to the KP alternate (effect) allele."""
import base64, json, numpy as np, pandas as pd

NHAP = 1006


def decode(b64):
    b = np.frombuffer(base64.b64decode(b64), dtype=np.uint8)
    return np.unpackbits(b, bitorder='little')[:NHAP].astype(np.float32)


def load(path):
    recs = json.load(open(path))
    df = pd.DataFrame(recs, columns=['chr', 'pos', 'id', 'ref', 'alt', 'ac', 'miss', 'bits', 'multi', 'vt'])
    df['chr'] = df.chr.astype(str); df['pos'] = df.pos.astype(int)
    return df


def match(varids, ref):
    """Match KP varIds (chr:pos:ref:alt, GRCh37) to panel records. Returns DataFrame with status, sign and haplotype vector index."""
    key = {}
    for i, r in enumerate(ref.itertuples()):
        key.setdefault((r.chr, r.pos), []).append(i)
    out = []
    for v in varids:
        c, p, a0, a1 = v.split(':'); p = int(p)
        cand = key.get((c, p), [])
        st, sign, idx = 'absent from panel', 0, -1
        for i in cand:
            r = ref.iloc[i]
            if r.ref == a0 and r.alt == a1: st, sign, idx = 'matched', 1, i; break
            if r.ref == a1 and r.alt == a0: st, sign, idx = 'matched (alleles swapped)', -1, i; break
        else:
            if cand: st = 'allele mismatch'
        if idx >= 0:
            ac = ref.iloc[idx].ac
            if ac == 0 or ac == NHAP: st = 'monomorphic in EUR'
            elif ref.iloc[idx].miss > 0: st = st + ' (missing calls)'
        pal = {a0, a1} in ({'A', 'T'}, {'C', 'G'})
        out.append(dict(varId=v, panel_status=st, sign=sign, idx=idx, palindromic=pal,
                        eur_af_alt=(ref.iloc[idx].ac / NHAP if sign == 1 else 1 - ref.iloc[idx].ac / NHAP) if idx >= 0 else np.nan))
    return pd.DataFrame(out)


def hap_matrix(m, ref):
    """rows = variants in m (usable ones), haplotypes coded as effect (KP alt) allele count."""
    H = np.vstack([decode(ref.iloc[i].bits) for i in m.idx])
    H[m.sign.values == -1] = 1 - H[m.sign.values == -1]
    return H


def corr(H):
    Hc = H - H.mean(1, keepdims=True)
    sd = np.sqrt((Hc ** 2).sum(1)); sd[sd == 0] = np.nan
    return (Hc @ Hc.T) / np.outer(sd, sd)
