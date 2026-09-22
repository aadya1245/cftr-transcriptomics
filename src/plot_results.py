"""Render three scientific figures from saved results, without rerunning DE."""
from pathlib import Path
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
ROOT=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
COLORS={'Nasal':'#007F73','Bronchial':'#C76539'}
def save(fig,name):
    fig.savefig(ROOT/'figures'/f'{name}.png',bbox_inches='tight')
    fig.savefig(ROOT/'figures'/f'{name}.svg',bbox_inches='tight')
    plt.close(fig)
def main():
    meta=pd.read_csv(ROOT/'data/processed/samples.tsv',sep='\t',index_col=0)
    norm=pd.read_csv(ROOT/'data/processed/normalized_counts.tsv.gz',sep='\t',index_col=0).loc[meta.index]
    de=pd.read_csv(ROOT/'results/differential_expression.tsv',sep='\t',index_col=0)
    log=np.log2(norm+1)
    top=log.var().nlargest(500).index
    model=PCA(n_components=2,svd_solver='full')
    xy=pd.DataFrame(model.fit_transform(log[top]),index=meta.index,columns=['PC1','PC2'])
    xy.to_csv(ROOT/'results/pca_scores.tsv',sep='\t')
    fig,ax=plt.subplots(figsize=(8,6))
    for donor,group in meta.groupby('donor'):
        points=xy.loc[group.index]
        ax.plot(points.PC1,points.PC2,color='#ccd3d6',lw=1,zorder=0)
    for tissue,color in COLORS.items():
        points=xy.loc[meta.index[meta.tissue==tissue]]
        ax.scatter(points.PC1,points.PC2,c=color,label=tissue,s=60,edgecolors='white')
    ax.set(xlabel=f'PC1 ({model.explained_variance_ratio_[0]:.1%} variance)',ylabel=f'PC2 ({model.explained_variance_ratio_[1]:.1%} variance)',
           title='Cultured airway samples: expression structure')
    ax.legend(frameon=False)
    fig.text(.12,.01,'Top 500 variable genes; log2(normalized counts + 1). Lines connect donors.\nExploratory PCA; no donor or batch residualization.',fontsize=9,color='#555555')
    fig.subplots_adjust(bottom=.18)
    save(fig,'01_pca')
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12,5.8),gridspec_kw={'width_ratios':[1.6,1]})
    valid=de.dropna(subset=['padj','log2FoldChange'])
    selected=(valid.padj<.05)&(valid.log2FoldChange.abs()>=1)
    colors=np.where(selected,np.where(valid.log2FoldChange>0,COLORS['Nasal'],COLORS['Bronchial']),'#c2c9cf')
    ax.scatter(valid.log2FoldChange,-np.log10(valid.padj.clip(lower=1e-300)),s=9,c=colors,alpha=.65,rasterized=True)
    ax.axhline(-np.log10(.05),color='#777777',ls='--',lw=.8)
    for v in [-1,1]: ax.axvline(v,color='#777777',ls=':',lw=.8)
    row=de.loc['CFTR']
    ax.scatter([row.log2FoldChange],[-np.log10(row.padj)],s=65,c='black',marker='*',zorder=3)
    ax.annotate('CFTR',(row.log2FoldChange,-np.log10(row.padj)),xytext=(8,12),textcoords='offset points')
    ax.set(xlabel='log2 fold change (nasal / bronchial)',ylabel='−log10(BH adjusted p)',title='A  Transcriptome-wide tissue contrast')
    for donor,group in meta.groupby('donor'):
        ordered=group.set_index('tissue').loc[['Bronchial','Nasal']]
        vals=norm.loc[ordered['gsm'].map(dict(zip(meta.gsm,meta.index))), 'CFTR'].to_numpy()
        bx.plot([0,1],np.log2(vals+1),color='#adb7be',lw=1)
        bx.scatter([0,1],np.log2(vals+1),c=[COLORS['Bronchial'],COLORS['Nasal']],s=35,zorder=3)
    bx.set(xticks=[0,1],xticklabels=['Bronchial','Nasal'],ylabel='log2(normalized CFTR count + 1)',title='B  CFTR within each donor')
    fig.text(.08,.015,f'Paired DESeq2 model: donor + tissue. Unshrunk effect estimates. CFTR log2FC = {row.log2FoldChange:.2f}; BH q = {row.padj:.3g}.\nColored genes: q < 0.05 and |log2FC| ≥ 1. The fold-change cutoff is descriptive, not a threshold-based hypothesis test.',fontsize=9)
    fig.tight_layout(rect=[0,.1,1,1])
    save(fig,'02_differential_expression')
    go=pd.read_csv(ROOT/'results/go_enrichment.tsv',sep='\t')
    fig,axes=plt.subplots(1,2,figsize=(14,7))
    for ax,direction,tissue in zip(axes,['Higher bronchial','Higher nasal'],['Bronchial','Nasal']):
        rows=go[(go.direction==direction)&(go.padj<.05)].head(5).iloc[::-1]
        if rows.empty:
            ax.text(.5,.5,'No terms pass q < 0.05',ha='center',transform=ax.transAxes)
        else:
            labels=['\n'.join(textwrap.wrap(t,38))+f'  [{k}/{K}]' for t,k,K in zip(rows.term,rows.k,rows.K)]
            ax.barh(labels,-np.log10(rows.padj.clip(lower=1e-300)),color=COLORS[tissue])
        ax.set(title=direction,xlabel='−log10(BH adjusted p)')
        ax.tick_params(axis='y',labelsize=9)
    fig.suptitle('GO biological process enrichment',fontsize=16)
    fig.text(.02,.02,'Top five significant terms per direction; brackets show overlap / term genes in background.\nBackground: genes with finite DE adjusted p. GO BP 2023; BH across all eligible terms and both directions. Terms may overlap.',fontsize=9)
    fig.tight_layout(rect=[0,.09,1,.95],w_pad=3)
    save(fig,'03_go_enrichment')
if __name__=='__main__': main()
