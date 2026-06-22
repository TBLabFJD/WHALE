process MERGE_DMR_FILES {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/biocontainers/pandas:1.5.2'

    input:
    tuple val(meta), path(tsv), path(dmr_bed), path(promoters), path(enhancers)

    output:
    tuple val(meta), path("${prefix}_merged.tsv"), emit: merged_tsv
    path "versions.yml"                          , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    cat << 'EOF' > merge_script.py
    import sys
    import pandas as pd
    import os

    tsv_file   = sys.argv[1]
    dmr_file   = sys.argv[2]
    prom_file  = sys.argv[3]
    enh_file   = sys.argv[4]
    out_prefix = sys.argv[5]

    ## load annotSV file
    df_tsv = pd.read_csv(tsv_file, sep='\t')
    
    ## keep only full annotation mode to avoid duplicates
    if 'Annotation_mode' in df_tsv.columns:
        df_tsv = df_tsv[df_tsv['Annotation_mode'] == 'full']

    ## change format coordinates for joining (annotsv coordinates are 0-based)
    df_tsv['join_chr'] = df_tsv['SV_chrom'].astype(str).apply(lambda x: x if x.startswith('chr') else 'chr' + x)
    df_tsv['join_start'] = df_tsv['SV_start'].astype(int) - 1 
    df_tsv['join_end'] = df_tsv['SV_end'].astype(int)

    def read_bed(filepath, suffix, keep_last_n=None):
        if os.path.getsize(filepath) > 0:
            df = pd.read_csv(filepath, sep='\t', header=None)
            
            df.rename(columns={0: 'join_chr', 1: 'join_start', 2: 'join_end'}, inplace=True)
            df['join_chr'] = df['join_chr'].astype(str).apply(lambda x: x if x.startswith('chr') else 'chr' + x)
            
            cols_to_keep = ['join_chr', 'join_start', 'join_end']
            
            if keep_last_n is not None:
                last_n_cols = list(df.columns)[-keep_last_n:]
                for c in last_n_cols:
                    if c not in cols_to_keep:
                        cols_to_keep.append(c)

            df = df[cols_to_keep]
            ## Prefix feature columns
            df.columns = [c if str(c).startswith('join_') else f"{suffix}_col{c}" for c in df.columns]
            return df
        else:
            return pd.DataFrame(columns=['join_chr', 'join_start', 'join_end'])

    ## reduce multiple overlaps into comma-separated strings
    def collapse_duplicates(df):
        if df.empty:
            return df
        
        group_cols = ['join_chr', 'join_start', 'join_end']
        val_cols = [c for c in df.columns if c not in group_cols]
        
        ## Group by region and join unique values with commas
        df_collapsed = df.groupby(group_cols)[val_cols].agg(
            lambda x: ','.join(sorted(list(set(x.dropna().astype(str)))))
        ).reset_index()
        
        return df_collapsed

    ## load and remove duplicates from the BED files
    df_dmr  = read_bed(dmr_file, 'dmr', keep_last_n=3)
    df_dmr  = collapse_duplicates(df_dmr)
    
    df_prom = read_bed(prom_file, 'promoter', keep_last_n=3)
    df_prom = collapse_duplicates(df_prom)
    
    df_enh  = read_bed(enh_file, 'enhancer', keep_last_n=3)
    df_enh  = collapse_duplicates(df_enh)

    ## left joins
    df_merged = pd.merge(df_tsv, df_dmr, on=['join_chr', 'join_start', 'join_end'], how='left')
    df_merged = pd.merge(df_merged, df_prom, on=['join_chr', 'join_start', 'join_end'], how='left')
    df_merged = pd.merge(df_merged, df_enh, on=['join_chr', 'join_start', 'join_end'], how='left')

    df_merged.drop(columns=['join_chr', 'join_start', 'join_end'], inplace=True)
    df_merged.to_csv(f"{out_prefix}_merged.tsv", sep='\t', index=False, na_rep='')
    EOF

    python3 merge_script.py ${tsv} ${dmr_bed} ${promoters} ${enhancers} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}