import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.edgecolor':'#c3c2b7','xtick.color':'#52514e','ytick.color':'#0b0b0b','axes.labelcolor':'#52514e'})
INK='#0b0b0b'; MUTED='#898781'; GRID='#e1e0d9'; BLUE='#2a78d6'; ORANGE='#eb6834'; AQUA='#1baf7a'; YEL='#eda100'
F=pd.read_csv('step3_mr_forward_results.csv'); Rv=pd.read_csv('step3_mr_reverse_results.csv'); S=pd.read_csv('step3_mr_forward_snp_level.csv')
fig=plt.figure(figsize=(11,8.2)); gs=fig.add_gridspec(2,3,height_ratios=[1.25,1],hspace=0.42,wspace=0.55)
gs2=fig.add_gridspec(2,5,height_ratios=[1.25,1],hspace=0.42,wspace=0.3)
# (a) forward forest: methods x outcomes (MVP EUR primary)
meths=['IVW (MRE)','MR-Egger','Weighted median','Weighted mode','MR-PRESSO (outlier-corrected)','IVW, Steiger-filtered','IVW, excluding TERC & TERT loci']
for j,(oc,xl) in enumerate([('SBP','SBP, SD per SD longer LTL'),('DBP','DBP, SD per SD longer LTL'),('Hypertension','Hypertension, log OR per SD longer LTL')]):
    ax=fig.add_subplot(gs[0,j]); rows=[]
    srcs=[s for s in ['MVP European (Verma 2024)','FinnGen','Genes & Health (Huang 2021)','Biobank Japan','MVP African American','MVP Hispanic'] if ((F.outcome==oc)&(F.source==s)).any()]
    y=0; yt=[];yl=[]
    for s in srcs:
        ms=meths if s=='MVP European (Verma 2024)' else ['IVW (MRE)']
        for mth in ms:
            r=F[(F.outcome==oc)&(F.source==s)&(F.method==mth)].iloc[0]
            c=INK if mth=='IVW (MRE)' else '#52514e'
            col={'MVP European (Verma 2024)':BLUE,'FinnGen':BLUE,'Genes & Health (Huang 2021)':ORANGE,'Biobank Japan':AQUA,'MVP African American':YEL,'MVP Hispanic':'#e87ba4'}[s]
            ax.plot([r.b-1.96*r.se,r.b+1.96*r.se],[-y,-y],color=col,lw=1.3); ax.plot(r.b,-y,'s' if mth=='IVW (MRE)' else 'o',ms=4.5 if mth=='IVW (MRE)' else 3.5,color=col,mec='white',mew=0.5)
            lab=(s.replace(' (Verma 2024)','').replace(' (Huang 2021)','')+': ' if mth=='IVW (MRE)' else '   ')+mth.replace(' (outlier-corrected)','').replace(' (MRE)','')
            yt.append(-y); yl.append(lab); y+=1
        y+=0.5
    ax.axvline(0,color='#c3c2b7',lw=0.8); ax.set_yticks(yt); ax.set_yticklabels(yl if j==0 else yl,fontsize=7); ax.tick_params(axis='y',length=0)
    for sp in ['top','right','left']: ax.spines[sp].set_visible(False)
    ax.grid(axis='x',color=GRID,lw=0.6); ax.set_axisbelow(True); ax.set_xlabel(xl,fontsize=7.5); ax.set_title(oc,loc='left',fontweight='bold',fontsize=9.5,color=INK)
fig.text(0.01,0.975,'a  Forward MR: genetically longer LTL (Codd 2021, 126 instruments) → blood pressure',fontweight='bold',fontsize=9.5)
# (b) scatter SBP MVP
ax=fig.add_subplot(gs2[1,0:2]); s=S[S.code=='S_mvpEU'].copy(); sg=np.sign(s.bx); s['bx']*=sg; s['by']*=sg
ax.errorbar(s.bx,s.by,xerr=1.96*s.sx,yerr=1.96*s.sy,fmt='none',ecolor='#e1e0d9',lw=0.6,zorder=1)
ax.scatter(s.bx,s.by,s=10,color=np.where(s.presso_outlier,MUTED,INK),zorder=2,lw=0)
hi=s[s.gene.isin(['TERT','MYNN','TERC'])]
for r in hi.itertuples(): ax.annotate(r.gene+' ('+r.rsid+')',(r.bx,r.by),fontsize=6.5,xytext=(3,3),textcoords='offset points',color=INK)
xx=np.linspace(0,s.bx.max()*1.05,10)
for mth,col,ls in [('IVW (MRE)',BLUE,'-'),('Weighted median',AQUA,'--')]:
    r=F[(F.code=='S_mvpEU')&(F.method==mth)].iloc[0]; ax.plot(xx,r.b*xx,color=col,ls=ls,lw=1.6,label=mth.replace(' (MRE)',''))
e=F[(F.code=='S_mvpEU')&(F.method=='MR-Egger')].iloc[0]; ax.plot(xx,e.intercept+e.b*xx,color=ORANGE,ls=':',lw=1.6,label='MR-Egger')
ax.axhline(0,color='#c3c2b7',lw=0.8); ax.legend(frameon=False,fontsize=7,loc='upper left')
ax.set_xlabel('SNP effect on LTL (Codd 2021)',fontsize=7.5); ax.set_ylabel('SNP effect on SBP (MVP European)',fontsize=7.5)
for sp in ['top','right']: ax.spines[sp].set_visible(False)
ax.set_title('b  SNP-level, SBP (MVP European)',loc='left',fontweight='bold',fontsize=9.5)
# (c) reverse MR
ax=fig.add_subplot(gs2[1,3:]); y=0; yt=[];yl=[]
for ex in ['SBP','DBP','Hypertension']:
    for oc,col in [('LTL (Codd 2021, UKB)',BLUE),('LTL South Asian (Nakao 2026)',ORANGE)]:
        for mth in ['IVW (MRE)','Weighted median','MR-PRESSO (outlier-corrected)']:
            r=Rv[(Rv.exposure==ex)&(Rv.outcome==oc)&(Rv.method==mth)].iloc[0]
            ax.plot([r.b-1.96*r.se,r.b+1.96*r.se],[-y,-y],color=col,lw=1.3); ax.plot(r.b,-y,'s' if mth.startswith('IVW') else 'o',ms=4,color=col,mec='white',mew=0.5)
            yt.append(-y); yl.append((f"{ex} → {'LTL (Codd, EUR)' if 'Codd' in oc else 'LTL (Nakao, SAS)'}: " if mth.startswith('IVW') else '')+mth.replace(' (outlier-corrected)','').replace(' (MRE)','')); y+=1
        y+=0.4
ax.axvline(0,color='#c3c2b7',lw=0.8); ax.set_yticks(yt); ax.set_yticklabels(yl,fontsize=6.8); ax.tick_params(axis='y',length=0)
for sp in ['top','right','left']: ax.spines[sp].set_visible(False)
ax.grid(axis='x',color=GRID,lw=0.6); ax.set_axisbelow(True); ax.set_xlabel('SD change in LTL per SD BP (or per log-odds hypertension)',fontsize=7.5)
ax.set_title('c  Reverse MR: BP (MVP European instruments) → LTL',loc='left',fontweight='bold',fontsize=9.5)
fig.text(0.01,0.005,'Points = estimate, bars = 95% CI. Exposure and outcome samples do not overlap (forward: UK Biobank → MVP/FinnGen/G&H/BBJ; reverse: MVP → UK Biobank). Grey points in (b) = MR-PRESSO outliers.',fontsize=6.8,color=MUTED)
for ext in ['png','pdf']: plt.savefig(f'figures/Fig4_mendelian_randomization.{ext}',dpi=300,bbox_inches='tight')
