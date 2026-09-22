"""Reproducible paired airway RNA-seq reanalysis. Run from repository root."""
from pathlib import Path
import gzip, json, re, hashlib, platform
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, false_discovery_control
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

ROOT = Path(__file__).resolve().parents[1]
PANEL = ['CFTR','SLC26A9','SLC6A14','SCNN1A','SCNN1B','SCNN1G','SLC9A3','MUC5AC','MUC5B','FOXJ1','KRT5']

def validate_counts(counts):
    x = counts.to_numpy()
    if not counts.index.is_unique or not counts.columns.is_unique:
        raise ValueError('Duplicate sample or gene identifiers')
    if not np.isfinite(x).all() or (x < 0).any() or (x != np.floor(x)).any():
        raise ValueError('Counts must be finite nonnegative integers')

def load_data():
    raw = pd.read_csv(ROOT/'data/raw/GSE172232_Raw_Gene_Count_Matrix.txt.gz', sep='\t')
    if not raw.Gene.equals(raw.Feature):
        raise ValueError('Gene/Feature annotations differ')
    counts = raw.drop(columns='Feature').set_index('Gene').T
    validate_counts(counts)
    blocks = gzip.open(ROOT/'data/raw/GSE172232_family.soft.gz','rt').read().split('^SAMPLE = ')[1:]
    records=[]
    for block in blocks:
        gsm=block.splitlines()[0]
        title=re.search(r'!Sample_title = (.*)',block)[1]
        m=re.fullmatch(r'Culture (Bronchial|Nasal) AWCTX_(\d+)',title)
        if not m: continue
        tissue, donor=m.groups()
        base=f'AWCTX_{donor}{tissue[0]}'
        matching=[c for c in counts.index if c in [base,base+'_V2']]
        if len(matching)!=1:
            raise ValueError(f'Ambiguous or missing match for {gsm}: {matching}')
        records.append(dict(sample=matching[0],gsm=gsm,donor='D'+donor,tissue=tissue,
                            mapping_basis='donor/tissue label convention; uppercase cultured; suffix V2 retained'))
    meta=pd.DataFrame(records).set_index('sample').sort_values(['donor','tissue'])
    pairs=meta.groupby('donor').tissue.agg(set)
    if not all(v=={'Nasal','Bronchial'} for v in pairs) or len(meta)!=2*len(pairs):
        raise ValueError('Incomplete or duplicated pairs')
    audit=pd.DataFrame({'sample':counts.index})
    audit['included']=audit['sample'].isin(meta.index)
    audit['reason']=audit['sample'].map(lambda c:'paired culture, donor-coded' if c in meta.index else
        'naive tissue outside contrast' if re.fullmatch(r'AWCTX_\d+[bn]',c) else
        'numeric label: mapping unresolved, excluded')
    audit.to_csv(ROOT/'data/processed/sample_audit.tsv',sep='\t',index=False)
    meta.to_csv(ROOT/'data/processed/samples.tsv',sep='\t')
    return counts.loc[meta.index].astype('int64'),meta

def enrich(result, libraries):
    universe=set(result.index[result.padj.notna()])
    rows=[]
    for direction, sign in [('Higher nasal',1),('Higher bronchial',-1)]:
        hits=set(result.index[(result.padj<.05)&(result.log2FoldChange*sign>=1)])
        for term,genes in libraries.items():
            members=genes&universe
            if not 15<=len(members)<=500: continue
            overlap=members&hits
            rows.append(dict(direction=direction,term=term,M=len(universe),N=len(hits),
                K=len(members),k=len(overlap),pvalue=float(hypergeom.sf(len(overlap)-1,len(universe),len(members),len(hits)))) if hits else
                dict(direction=direction,term=term,M=len(universe),N=0,K=len(members),k=0,pvalue=1.))
            rows[-1]['overlap_genes']=';'.join(sorted(overlap))
    out=pd.DataFrame(rows)
    out['padj']=false_discovery_control(out.pvalue.to_numpy(),method='bh')
    return out.sort_values('padj')

def main():
    counts,meta=load_data()
    keep=(counts>=10).sum(axis=0)>=len(meta)//2
    filtered=counts.loc[:,keep]
    inference=DefaultInference(n_cpus=2)
    dds=DeseqDataSet(counts=filtered,metadata=meta[['donor','tissue']],design='~ donor + tissue',
                    refit_cooks=True,inference=inference)
    design=np.asarray(dds.obsm['design_matrix'])
    if np.linalg.matrix_rank(design)!=design.shape[1]: raise ValueError('Rank-deficient design')
    dds.deseq2()
    stats=DeseqStats(dds,contrast=['tissue','Nasal','Bronchial'],alpha=.05,
                    cooks_filter=True,independent_filter=False,inference=inference)
    stats.summary()
    result=stats.results_df
    result.index.name='gene'
    result['ci95_low']=result.log2FoldChange-1.96*result.lfcSE
    result['ci95_high']=result.log2FoldChange+1.96*result.lfcSE
    result.to_csv(ROOT/'results/differential_expression.tsv',sep='\t')
    result.reindex(PANEL).to_csv(ROOT/'results/cftr_context_panel.tsv',sep='\t')
    pd.DataFrame(dds.layers['normed_counts'],index=meta.index,columns=filtered.columns).to_csv(
        ROOT/'data/processed/normalized_counts.tsv.gz',sep='\t',compression='gzip')
    qc=pd.DataFrame({'library_size':counts.sum(axis=1),'detected_genes':(counts>0).sum(axis=1),
                     'size_factor':dds.obs['size_factors']})
    qc.to_csv(ROOT/'results/sample_qc.tsv',sep='\t')
    libraries={}
    for line in (ROOT/'data/raw/GO_Biological_Process_2023.gmt').read_text().splitlines():
        term,_,*genes=line.split('\t')
        libraries[term]=set(genes)
    go=enrich(result,libraries)
    go.to_csv(ROOT/'results/go_enrichment.tsv',sep='\t',index=False)
    summary=dict(samples=len(meta),donors=meta.donor.nunique(),input_genes=counts.shape[1],
                 filtered_genes=filtered.shape[1],tested_genes=int(result.padj.notna().sum()),
                 significant_genes=int((result.padj<.05).sum()),
                 higher_nasal=int(((result.padj<.05)&(result.log2FoldChange>=1)).sum()),
                 higher_bronchial=int(((result.padj<.05)&(result.log2FoldChange<=-1)).sum()),
                 enrichment_tests=len(go),significant_go_tests=int((go.padj<.05).sum()),
                 design_rank=int(np.linalg.matrix_rank(design)),python=platform.python_version())
    (ROOT/'results/summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
