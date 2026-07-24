#!/usr/bin/env python3

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy import stats
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False


def log(x):
    print(x, flush=True)


def clean_string(x):
    out = x.astype('string').str.strip()
    return out.replace({'': pd.NA, 'NA': pd.NA, 'NaN': pd.NA, 'nan': pd.NA,
                        'None': pd.NA, 'null': pd.NA, '<NA>': pd.NA})


def normalize_study(x):
    raw = clean_string(x)
    upper = raw.str.upper()
    out = raw.copy()
    out.loc[upper.str.startswith('AIBL', na=False)] = 'AIBL'
    out.loc[upper.str.startswith('BLSA', na=False)] = 'BLSA'
    out.loc[upper.str.startswith('OASIS', na=False)] = 'OASIS'
    return out


def normalize_apoe(x):
    raw = clean_string(x).str.upper()
    raw = raw.str.replace('APOE', '', regex=False)
    raw = raw.str.replace(' ', '', regex=False)
    raw = raw.str.replace('-', '/', regex=False)
    raw = raw.str.replace('_', '/', regex=False)
    replacements = {
        '2/2': 'E2/E2', '2/3': 'E2/E3', '3/2': 'E2/E3',
        '2/4': 'E2/E4', '4/2': 'E2/E4', '3/3': 'E3/E3',
        '3/4': 'E3/E4', '4/3': 'E3/E4', '4/4': 'E4/E4',
        'E2E2': 'E2/E2', 'E2E3': 'E2/E3', 'E3E2': 'E2/E3',
        'E2E4': 'E2/E4', 'E4E2': 'E2/E4', 'E3E3': 'E3/E3',
        'E3E4': 'E3/E4', 'E4E3': 'E3/E4', 'E4E4': 'E4/E4',
        '22': 'E2/E2', '23': 'E2/E3', '32': 'E2/E3',
        '24': 'E2/E4', '42': 'E2/E4', '33': 'E3/E3',
        '34': 'E3/E4', '43': 'E3/E4', '44': 'E4/E4',
    }
    raw = raw.replace(replacements)
    keep = ['E2/E2', 'E2/E3', 'E2/E4', 'E3/E3', 'E3/E4', 'E4/E4']
    out = pd.Series(pd.NA, index=raw.index, dtype='string')
    out.loc[raw.isin(keep)] = raw
    return out


def normalize_dx(x):
    raw = clean_string(x).str.upper()
    out = pd.Series(pd.NA, index=raw.index, dtype='string')
    out.loc[raw.isin(['CN', 'NC', 'NORMAL', 'COGNITIVELY NORMAL',
                     'COGNITIVE NORMAL', 'CONTROL', 'HEALTHY CONTROL',
                     'HC', '0'])] = 'CN'
    out.loc[raw.isin(['MCI', 'LMCI', 'EMCI', 'EARLY MCI',
                     'MILD COGNITIVE IMPAIRMENT', '1'])] = 'MCI'
    out.loc[raw.isin(['AD', 'DEMENTIA', 'ALZHEIMER', "ALZHEIMER'S DISEASE",
                     'ALZHEIMERS DISEASE', 'ALZHEIMER DISEASE', '2'])] = 'AD'
    unresolved = out.isna() & raw.notna()
    out.loc[unresolved & raw.str.contains(r'(^|[^A-Z])MCI([^A-Z]|$)', regex=True, na=False)] = 'MCI'
    unresolved = out.isna() & raw.notna()
    out.loc[unresolved & raw.str.contains('ALZHEIMER|DEMENTIA', regex=True, na=False)] = 'AD'
    unresolved = out.isna() & raw.notna()
    out.loc[unresolved & raw.str.contains('COGNITIVELY NORMAL|COGNITIVE NORMAL', regex=True, na=False)] = 'CN'
    return out


def find_col(columns, preferred, regex, label, required=True):
    for col in preferred:
        if col in columns:
            return col
    matches = [c for c in columns if regex and re.search(regex, c, re.I)]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(f'Multiple candidate columns for {label}: {matches}')
    if required:
        raise ValueError(f'Could not find {label}')
    return None


def p_from_t(t, df):
    if not np.isfinite(t) or not np.isfinite(df) or df <= 0:
        return np.nan
    if HAVE_SCIPY:
        return float(2 * stats.t.sf(abs(t), df))
    return float(math.erfc(abs(t) / math.sqrt(2)))


def bh_fdr(series):
    p = pd.to_numeric(series, errors='coerce')
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna() & np.isfinite(p)
    if not valid.any():
        return out
    vals = p.loc[valid].to_numpy(float)
    order = np.argsort(vals)
    ranked = vals[order]
    m = len(ranked)
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    restored = np.empty(m)
    restored[order] = adj
    out.loc[valid] = restored
    return out


def select_baseline(sample):
    d = sample.copy()
    d['has_dx'] = d['baseline_diagnosis'].notna()
    d['date_missing'] = d['baseline_date'].isna()
    d['age_missing'] = d['baseline_age'].isna()
    d = d.sort_values(
        ['study', 'participant_id', 'has_dx', 'date_missing', 'baseline_date',
         'age_missing', 'baseline_age', 'sample_source_row'],
        ascending=[True, True, False, True, True, True, True, True],
        kind='mergesort'
    )
    d = d.groupby(['study', 'participant_id'], as_index=False, sort=False).head(1).copy()
    d['baseline_selection_source'] = np.select(
        [d['has_dx'] & d['baseline_date'].notna(),
         d['has_dx'] & d['baseline_date'].isna() & d['baseline_age'].notna(),
         ~d['has_dx'] & d['baseline_date'].notna(),
         ~d['has_dx'] & d['baseline_date'].isna() & d['baseline_age'].notna()],
        ['earliest mapped-diagnosis Date', 'youngest mapped-diagnosis Age',
         'fallback earliest APOE-record Date', 'fallback youngest APOE-record Age'],
        default='fallback first APOE record'
    )
    return d.drop(columns=['has_dx', 'date_missing', 'age_missing'])


def match_prediction(baseline, pred):
    d = pred.merge(baseline, on=['study', 'participant_id'], how='inner', validate='many_to_one')
    d['date_difference_days'] = (d['prediction_date'] - d['baseline_date']).abs().dt.total_seconds() / 86400
    d['age_difference_years'] = (d['prediction_age'] - d['baseline_age']).abs()
    d['priority'] = np.select(
        [d['date_difference_days'].notna(), d['age_difference_years'].notna(),
         d['prediction_date'].notna(), d['prediction_age'].notna()],
        [1, 2, 3, 4], default=5
    )
    d['distance'] = np.select(
        [d['priority'].eq(1), d['priority'].eq(2), d['priority'].eq(3), d['priority'].eq(4)],
        [d['date_difference_days'], d['age_difference_years'],
         d['prediction_date'].map(lambda x: x.toordinal() if pd.notna(x) else np.nan),
         d['prediction_age']],
        default=d['prediction_source_row'].astype(float)
    )
    d['prediction_match_source'] = np.select(
        [d['priority'].eq(1), d['priority'].eq(2), d['priority'].eq(3), d['priority'].eq(4)],
        ['closest scored scan by Date', 'closest scored scan by Age',
         'earliest scored scan by Date', 'youngest scored scan by Age'],
        default='first available scored scan'
    )
    d = d.sort_values(['study', 'participant_id', 'priority', 'distance',
                       'prediction_date', 'prediction_age', 'prediction_source_row'],
                      kind='mergesort')
    return d.groupby(['study', 'participant_id'], as_index=False, sort=False).head(1).copy()


def make_design(data, exposure, numeric_covs, categorical_covs):
    parts = [pd.DataFrame({exposure: pd.to_numeric(data[exposure], errors='coerce')}, index=data.index)]
    terms = [exposure]
    for cov in numeric_covs:
        if cov in data:
            x = pd.to_numeric(data[cov], errors='coerce')
            if x.dropna().nunique() >= 2:
                parts.append(pd.DataFrame({cov: x}, index=data.index))
                terms.append(cov)
    for cov in categorical_covs:
        if cov in data:
            x = clean_string(data[cov])
            if x.dropna().nunique() >= 2:
                dummy = pd.get_dummies(x, prefix=cov, drop_first=True, dtype=float)
                if dummy.shape[1] > 0:
                    parts.append(dummy)
                    terms.extend(dummy.columns.tolist())
    return pd.concat(parts, axis=1), terms


def fit_ols(data, outcome, exposure, numeric_covs, categorical_covs, standardized, model_name):
    y0 = pd.to_numeric(data[outcome], errors='coerce')
    X0, terms = make_design(data, exposure, numeric_covs, categorical_covs)
    d = pd.concat([y0.rename(outcome), X0], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    n = len(d)
    n2 = int((d[exposure] == 0).sum()) if exposure in d else 0
    n4 = int((d[exposure] == 1).sum()) if exposure in d else 0
    res = dict(model=model_name, outcome=outcome, standardized_outcome=standardized,
               n=n, n_e2e2=n2, n_e4e4=n4, model_terms=','.join(terms),
               mean_e2e2=np.nan, mean_e4e4=np.nan,
               mean_difference_e4e4_minus_e2e2=np.nan,
               beta_e4e4_vs_e2e2=np.nan, se=np.nan, t=np.nan, p=np.nan,
               ci_low=np.nan, ci_high=np.nan, r2=np.nan, df_resid=np.nan,
               status='not_fit')
    if n == 0 or n2 == 0 or n4 == 0:
        res['status'] = 'missing one APOE group'
        return res
    y_raw = d[outcome].astype(float)
    expo = d[exposure].astype(float)
    res['mean_e2e2'] = float(y_raw[expo == 0].mean())
    res['mean_e4e4'] = float(y_raw[expo == 1].mean())
    res['mean_difference_e4e4_minus_e2e2'] = res['mean_e4e4'] - res['mean_e2e2']
    y = y_raw.copy()
    if standardized:
        sd = y.std(ddof=1)
        if not np.isfinite(sd) or sd <= 0:
            res['status'] = 'invalid outcome SD'
            return res
        y = (y - y.mean()) / sd
    X = d[X0.columns].astype(float).to_numpy()
    X = np.column_stack([np.ones(n), X])
    df_resid = n - X.shape[1]
    res['df_resid'] = df_resid
    if df_resid <= 0:
        res['status'] = 'nonpositive residual df'
        return res
    try:
        bhat = np.linalg.lstsq(X, y.to_numpy(float), rcond=None)[0]
        resid = y.to_numpy(float) - X @ bhat
        rss = float(resid.T @ resid)
        tss = float(((y.to_numpy(float) - y.mean()) ** 2).sum())
        cov = (rss / df_resid) * np.linalg.pinv(X.T @ X)
        se = float(np.sqrt(np.diag(cov))[1])
        beta = float(bhat[1])
        tval = beta / se if se > 0 else np.nan
        pval = p_from_t(tval, df_resid)
        crit = float(stats.t.ppf(0.975, df_resid)) if HAVE_SCIPY else 1.96
        res.update(beta_e4e4_vs_e2e2=beta, se=se, t=tval, p=pval,
                   ci_low=beta - crit * se, ci_high=beta + crit * se,
                   r2=1 - rss / tss if tss > 0 else np.nan, status='ok')
    except Exception as exc:
        res['status'] = f'OLS failed: {exc}'
    return res


def read_sample(path, studies):
    raw = pd.read_csv(path, sep='\t', low_memory=False)
    cols = raw.columns.tolist()
    idc = find_col(cols, ['PTID', 'participant_id', 'IID'], r'(^PTID$|participant.*id|^IID$)', 'sample ID')
    studyc = find_col(cols, ['Study', 'STUDY'], r'(^|_)study$', 'Study')
    apoec = find_col(cols, ['APOE_Genotype', 'APOE_genotype'], r'APOE.*genotype', 'APOE genotype')
    agec = find_col(cols, ['Age', 'AGE'], r'^age$', 'Age')
    sexc = find_col(cols, ['Sex', 'SEX', 'sex'], r'^sex$', 'Sex', False)
    datec = find_col(cols, ['Date', 'scan_date', 'MRI_Date'], r'(^|_)date$', 'Date', False)
    sitec = find_col(cols, ['SITE', 'Site', 'external_SITE'], r'(^|_)site$', 'SITE', False)
    dxc = find_col(cols, ['DX_Binary', 'Dx_binary', 'dx_binary', 'Diagnosis'], r'(^|_)dx[_\.]*binary$', 'diagnosis', False)
    d = pd.DataFrame({
        'sample_source_row': np.arange(len(raw)),
        'participant_id': clean_string(raw[idc]),
        'study': normalize_study(raw[studyc]),
        'APOE_genotype': normalize_apoe(raw[apoec]),
        'baseline_age': pd.to_numeric(raw[agec], errors='coerce'),
        'baseline_sex': clean_string(raw[sexc]) if sexc else pd.Series(pd.NA, index=raw.index, dtype='string'),
        'baseline_date': pd.to_datetime(raw[datec], errors='coerce') if datec else pd.Series(pd.NaT, index=raw.index),
        'baseline_site': clean_string(raw[sitec]) if sitec else pd.Series('Unknown', index=raw.index, dtype='string'),
        'baseline_diagnosis': normalize_dx(raw[dxc]) if dxc else pd.Series(pd.NA, index=raw.index, dtype='string'),
    })
    d = d[d['study'].isin(studies) & d['participant_id'].notna() & d['APOE_genotype'].isin(['E2/E2', 'E4/E4'])].copy()
    d['APOE_e4e4_vs_e2e2'] = np.where(d['APOE_genotype'].eq('E4/E4'), 1.0, 0.0)
    return d


def read_predictions(path, studies):
    raw = pd.read_csv(path, sep='\t', low_memory=False)
    cols = raw.columns.tolist()
    idc = find_col(cols, ['PTID', 'participant_id', 'IID'], r'(^PTID$|participant.*id|^IID$)', 'prediction ID')
    studyc = find_col(cols, ['external_Study', 'Study', 'STUDY'], r'(^|_)study$', 'prediction Study')
    agec = find_col(cols, ['Age', 'age_at_scan_used_for_model', 'AGE'], r'(^|_)age($|_at_scan)', 'prediction Age')
    datec = find_col(cols, ['Date', 'scan_date', 'MRI_Date'], r'(^|_)date$', 'prediction Date', False)
    sitec = find_col(cols, ['external_SITE', 'SITE', 'Site'], r'(^|_)site$', 'prediction SITE', False)
    yrc = find_col(cols, ['adni_brain_mri_ad_epoch_acceleration_years', 'adni_brain_mri_ad_lepoch_acceleration_years'], r'acceleration[_\.]*years$', 'acceleration years')
    zc = find_col(cols, ['adni_brain_mri_ad_epoch_acceleration_z', 'adni_brain_mri_ad_lepoch_acceleration_z'], r'acceleration[_\.]*z$', 'acceleration z', False)
    rc = find_col(cols, ['adni_brain_mri_ad_epoch_risk_score', 'adni_brain_mri_ad_lepoch_risk_score'], r'risk[_\.]*score$', 'risk score', False)
    d = pd.DataFrame({
        'prediction_source_row': np.arange(len(raw)),
        'participant_id': clean_string(raw[idc]),
        'study': normalize_study(raw[studyc]),
        'prediction_age': pd.to_numeric(raw[agec], errors='coerce'),
        'prediction_date': pd.to_datetime(raw[datec], errors='coerce') if datec else pd.Series(pd.NaT, index=raw.index),
        'prediction_site': clean_string(raw[sitec]) if sitec else pd.Series('Unknown', index=raw.index, dtype='string'),
        'ad_epoch_acceleration_years': pd.to_numeric(raw[yrc], errors='coerce'),
    })
    outcomes = ['ad_epoch_acceleration_years']
    if zc:
        d['ad_epoch_acceleration_z'] = pd.to_numeric(raw[zc], errors='coerce')
        outcomes.append('ad_epoch_acceleration_z')
    if rc:
        d['ad_epoch_risk_score'] = pd.to_numeric(raw[rc], errors='coerce')
        outcomes.append('ad_epoch_risk_score')
    d = d[d['study'].isin(studies) & d['participant_id'].notna() & d['ad_epoch_acceleration_years'].notna()].copy()
    return d, outcomes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample-file', required=True)
    ap.add_argument('--prediction-file', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--analysis-group', choices=['AIBL', 'BLSA', 'OASIS', 'POOLED'], required=True)
    ap.add_argument('--studies', nargs='+', default=['AIBL', 'BLSA', 'OASIS'])
    ap.add_argument('--primary-outcome', default='ad_epoch_acceleration_years')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    studies = list(args.studies)

    sample = read_sample(Path(args.sample_file), studies)
    pred, outcomes = read_predictions(Path(args.prediction_file), studies)
    baseline = select_baseline(sample)
    matched = match_prediction(baseline, pred)
    matched['analysis_site'] = matched['prediction_site'].fillna(matched['baseline_site'])
    matched['analysis_age'] = matched['prediction_age'].fillna(matched['baseline_age'])
    matched['analysis_sex'] = matched['baseline_sex']

    if args.analysis_group == 'POOLED':
        data = matched.copy()
    else:
        data = matched[matched['study'].eq(args.analysis_group)].copy()
    if data.empty:
        raise ValueError(f'No matched participants for {args.analysis_group}')

    prefix = f'{args.analysis_group}_APOE_E4E4_vs_E2E2_baseline_AD_EPOCH'
    exposure = 'APOE_e4e4_vs_e2e2'
    specs = [
        ('unadjusted', [], []),
        ('age_sex_adjusted', ['analysis_age'], ['analysis_sex']),
        ('extended_adjusted', ['analysis_age'], ['analysis_sex', 'baseline_diagnosis', 'analysis_site'] + (['study'] if args.analysis_group == 'POOLED' else [])),
    ]
    rows = []
    for outcome in outcomes:
        for name, numeric, categorical in specs:
            for standardized in [False, True]:
                r = fit_ols(data, outcome, exposure, numeric, categorical, standardized, name)
                r.update(analysis_group=args.analysis_group,
                         included_studies=','.join(sorted(data['study'].dropna().unique())),
                         reference_group='APOE E2/E2', comparison_group='APOE E4/E4',
                         effect_direction='positive beta = higher AD EPOCH in E4/E4')
                rows.append(r)
    results = pd.DataFrame(rows)
    if not results.empty:
        results['p_fdr_bh_within_job'] = results.groupby('standardized_outcome')['p'].transform(bh_fdr)

    data.to_csv(outdir / f'{prefix}_matched_baseline_data.tsv', sep='\t', index=False)
    results.to_csv(outdir / f'{prefix}_association_results.tsv', sep='\t', index=False)
    data.groupby(['study', 'APOE_genotype'], dropna=False).size().reset_index(name='N').to_csv(
        outdir / f'{prefix}_genotype_counts.tsv', sep='\t', index=False)
    data.groupby(['study', 'prediction_match_source'], dropna=False).agg(
        n_participants=('participant_id', 'nunique'),
        median_date_difference_days=('date_difference_days', 'median'),
        median_age_difference_years=('age_difference_years', 'median')
    ).reset_index().to_csv(outdir / f'{prefix}_matching_QC.tsv', sep='\t', index=False)
    data.groupby(['study', 'baseline_selection_source'], dropna=False).size().reset_index(name='N').to_csv(
        outdir / f'{prefix}_baseline_selection_QC.tsv', sep='\t', index=False)

    summary = {
        'analysis_group': args.analysis_group,
        'n_matched': int(len(data)),
        'n_e2e2': int(data['APOE_genotype'].eq('E2/E2').sum()),
        'n_e4e4': int(data['APOE_genotype'].eq('E4/E4').sum()),
        'outcomes': outcomes,
    }
    (outdir / f'{prefix}_run_summary.json').write_text(json.dumps(summary, indent=2))

    log('=' * 72)
    log(json.dumps(summary, indent=2))
    primary = results[(results['outcome'] == args.primary_outcome) &
                      (results['model'] == 'extended_adjusted') &
                      (results['standardized_outcome'] == True)]
    log(primary.to_string(index=False) if not primary.empty else 'Primary result unavailable')
    log('=' * 72)


if __name__ == '__main__':
    main()
