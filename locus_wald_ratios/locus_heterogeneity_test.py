"""Region-excluded IVW comparators, overdispersion factors and heterogeneity-adjusted locus tests (v24).
Run in mr_v4_inputs/ of the analysis package (snp_level_*.csv). Pure Python."""
import csv, math
def load(f): return list(csv.DictReader(open(f)))
def ivw(rows):
    bx=[float(r['bx']) for r in rows]; by=[float(r['by']) for r in rows]; sy=[float(r['sy']) for r in rows]
    W=sum(x*x/s/s for x,s in zip(bx,sy)); b=sum(x*y/s/s for x,y,s in zip(bx,by,sy))/W
    k=len(rows); Q=sum((y-b*x)**2/s/s for x,y,s in zip(bx,by,sy)); phi=max(1,math.sqrt(Q/(k-1)))
    return b,phi/math.sqrt(W),Q,k,phi
def p2(z): return math.erfc(abs(z)/math.sqrt(2))
W=[('TERC rs2293607 SBP','S',0.351354,0.032029),('TERC rs2293607 DBP','D',0.097830,0.030190),('TERC rs2293607 HTN MVP','H',0.312862,0.059946),
   ('TERC rs2293607 HTN FinnGen','F',0.223280,0.060579),('TERC rs12638862 SBP','S',0.383990,0.034617),('TERT rs7705526 SBP','S',0.240715,0.036123),
   ('TERT rs7705526 DBP','D',0.166619,0.034912),('TERT rs7705526 HTN MVP','H',0.186911,0.062306),('TERT rs7705526 HTN FinnGen','F',0.378866,0.070425)]
files={'S':'snp_level_S_mvpEU.csv','D':'snp_level_D_mvpEU.csv','H':'snp_level_H_mvpEU.csv','F':'snp_level_H_fg12.csv'}
res={}
for c,f in files.items():
    rows=load(f); res[c]=ivw([r for r in rows if r['terc_tert']=='False'])
    print(c,'excl b=%.5f se=%.5f Q=%.1f k=%d phi=%.3f'%res[c])
for lab,c,r,se in W:
    b,seb,Q,k,phi=res[c]
    z1=(r-b)/math.sqrt(se*se+seb*seb); z2=(r-b)/math.sqrt((phi*se)**2+seb*seb)
    print('%s comparator %.4f z=%.2f P=%.2g z_het=%.2f P_het=%.2g'%(lab,b,z1,p2(z1),z2,p2(z2)))
