"""Render explicit current-run batch-correction evidence for the cohort report; no results-tree search.

One bundle per corrected rank (written by Modules/batch_correction.nf): the quality, settings, methods,
sample-status, PERMANOVA and status files plus a manifest of their SHA-256. The panel shows what the run
did (corrected, or not and why), batch and covariate R2 before and after on the same samples, distance
and evaluation mode, and the covariates that need a look. It makes no claim that biological signal was
preserved: there is no validated criterion for that yet.

"Review needed" on a covariate (a provisional rule): its share of the non-batch variation, R2 / (1 - batch
R2), fell by at least 20% while p < 0.05 before correction (the existing pipeline warning level, not a
validated threshold), or rose by at least 20% and 0.01 with p < 0.05 after, or its p-value moved from
< 0.05 before correction to >= 0.05 after. Using the share
avoids reading the rise that removing batch variation causes in every other R2 as a change. Every
covariate's change is shown either way.
"""
import csv
import hashlib
import html
import math
from pathlib import Path

FILES = ['batch_correction_quality.tsv', 'batch_correction_run_settings.tsv', 'batch_correction_methods.tsv',
         'batch_correction_sample_status.tsv', 'permanova_summary.tsv', 'batch_correction_effect_sizes.tsv',
         'batch_correction_dispersion.tsv', 'BATCH_CORRECTION_STATUS.tsv']
AFTER_SOURCES = ['ConQuR_tuned', 'UNCORRECTED_by_decision', 'UNCORRECTED_passthrough']
REVIEW_DROP = 0.20      # provisional; the pipeline's existing covariate-warning level
REVIEW_ALPHA = 0.05
REVIEW_RISE_ABS = 0.01   # provisional: a rise in share must also be at least this large
STRONG_LOG2 = 1.0       # provisional: taxa with at least a two-fold difference before correction
MIXED_WARNING = ('In a synthetic check (2026-11-17), leaving sample type unprotected while it overlapped with batch '
                 'shrank a known tissue difference, and the covariate R2 did not show it.')
DEPTH_NOTE = ('If batches differ in depth, correction can remove part of a true biological difference. Adding depth as a '
              'covariate (--batch_corr_depth_covariate true) keeps more of it but leaves more batch effect when depth '
              'tracks batch (synthetic check, 2026-11-19: 9–11% of the batch R2 left at a 1.77x depth gap, 36–45% at 2.65x).')
NAMED_COMPARISON ='Direction and size of a named clinical comparison (none is specified for this run).'
NOT_ASSESSED = [
    'Whether biological signal was preserved (no validated criterion yet; the flags above only point to covariates worth a look).',
    'Other correction methods or settings on this data.',
    'Compositional metrics (for example Aitchison distance) and stratification by library size.',
    'Patient-restricted permutations.',
    NAMED_COMPARISON,
]


def read(path):
    with Path(path).open() as handle:
        return [{k: ('NA' if v in ('', None) else v) for k, v in row.items()} for row in csv.DictReader(handle, delimiter='\t')]


def value(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def esc(x):
    return html.escape(str(x))


def fmt(x, digits=4):
    return 'NA' if x is None else f'{x:.{digits}f}'


def load(path):
    path = Path(path)
    pairs = read(path / 'manifest.tsv')
    m = {r['key']: r['value'] for r in pairs}
    if len(m) != len(pairs) or m.get('report_schema') != '1' or m.get('rank') in (None, 'NA'):
        raise ValueError('Invalid batch-correction evidence manifest')
    data = {}
    for name in FILES:
        f = path / name
        if hashlib.sha256(f.read_bytes()).hexdigest() != m.get('report_sha256.' + name):
            raise ValueError(f'Batch-correction evidence hash mismatch: {name}')
        data[name] = read(f)
    quality = data['batch_correction_quality.tsv']
    status = data['BATCH_CORRECTION_STATUS.tsv']
    if len(quality) != 1 or len(status) != 1:
        raise ValueError('Expected one quality row and one status row')
    settings = {r['key']: r['value'] for r in data['batch_correction_run_settings.tsv']}
    included = [r for r in data['batch_correction_sample_status.tsv'] if r['status'] == 'included']
    n = value(quality[0]['n_samples'])
    if n is None or int(n) != len(included):
        raise ValueError('Quality sample count disagrees with the sample ledger')
    counts = [r for r in data['permanova_summary.tsv'] if r['Transform'] == 'counts']
    if any(int(value(r['N']) or -1) != len(included) for r in counts):
        raise ValueError('PERMANOVA sample count disagrees with the sample ledger')
    return m, data, settings


def term_rows(perm, batch_column, covariates, transform='counts'):
    sources = {r['Source'] for r in perm if r['Transform'] == transform}
    after_source = next((s for s in AFTER_SOURCES if s in sources), None)
    pick = lambda src, term: next((r for r in perm if r['Source'] == src and r['Transform'] == transform and r['Term'] == term), None)
    rows = []
    batch_b = batch_a = None   # batch R2 before/after, to express covariates as a share of the non-batch variation
    for term in [batch_column] + covariates:
        b, a = pick('raw', term), (pick(after_source, term) if after_source else None)
        r2b, r2a = value(b['R2']) if b else None, value(a['R2']) if a else None
        pb, pa = value(b['Pr(>F)']) if b else None, value(a['Pr(>F)']) if a else None
        modes = {x['Evaluation'] for x in (b, a) if x}
        ppb = value(b.get('p_within_patient')) if b else None
        ppa = value(a.get('p_within_patient')) if a else None
        row = {'term': term, 'role': 'batch' if term == batch_column else 'covariate',
               'r2_before': r2b, 'r2_after': r2a, 'p_before': pb, 'p_after': pa, 'pp_before': ppb, 'pp_after': ppa,
               'change': None, 'share_before': None, 'share_after': None, 'note': '', 'review': False}
        if b is None:
            row['note'] = 'not evaluated before correction'
        elif a is None:
            row['note'] = 'no corrected or retained table to compare'
        elif len(modes) > 1:
            row['note'] = 'not directly comparable: evaluation mode differs (' + ', '.join(sorted(modes)) + ')'
        elif r2b is None or r2a is None:
            row['note'] = 'R2 not computable; not compared'
        else:
            row['change'] = r2a - r2b
            if row['role'] == 'batch':
                batch_b, batch_a = r2b, r2a
            else:
                # Removing batch variation shrinks the total, so every other R2 grows by itself. Judge a
                # covariate on its share of the non-batch variation, R2 / (1 - batch R2), when batch R2 is known.
                if batch_b is not None and batch_a is not None and batch_b < 1 and batch_a < 1:
                    sb, sa, what = r2b / (1 - batch_b), r2a / (1 - batch_a), 'share of non-batch variation'
                    row['share_before'], row['share_after'] = sb, sa
                else:
                    sb, sa, what = r2b, r2a, 'R2'
                # A fall matters only for a covariate that mattered before; a rise only if it matters after.
                dropped = (sb > 0 and (sb - sa) / sb >= REVIEW_DROP and pb is not None and pb < REVIEW_ALPHA)
                rose = (sa - sb >= REVIEW_RISE_ABS and (sb <= 0 or (sa - sb) / sb >= REVIEW_DROP)
                        and pa is not None and pa < REVIEW_ALPHA)
                lost = pb is not None and pa is not None and pb < REVIEW_ALPHA <= pa
                row['review'] = dropped or rose or lost
                reasons = [f'{what} fell by {100 * (sb - sa) / sb:.0f}% (provisional 20% rule)'] if dropped else []
                reasons += [f'{what} rose from {sb:.3g} to {sa:.3g} (provisional rule: 20% and 0.01)'] if rose else []
                reasons += [f'p moved from {pb:.3g} to {pa:.3g}'] if lost else []
                row['note'] = '; '.join(reasons)
        rows.append(row)
    return rows, after_source


def summary_line(m, quality, status, settings):
    q, st = quality[0], status[0]
    state = st['status']
    if state == 'corrected':
        what = f"Corrected with ConQuR ({esc(q.get('method', 'NA'))})."
    elif state.startswith('uncorrected_'):
        what = f"Correction not applied ({esc(state)}): {esc(q.get('correction_decision_reason', st.get('reason', 'NA')))}. The input table is used downstream."
    else:
        what = f"Status {esc(state)}: {esc(st.get('reason', 'NA'))}."
    parts = [what, f"Samples {esc(q['n_samples'])}, taxa {esc(q['n_taxa'])}.",
             f"Covariates: {esc(settings.get('covariates', 'NA'))}.",
             f"Decontamination controls: {esc(settings.get('decontamination_controls', 'not recorded'))}.",
             f"ConQuR {esc(settings.get('ConQuR_version', 'NA'))}, vegan {esc(settings.get('vegan_version', 'NA'))}."]
    # The correction decision uses relative-abundance Bray; the table below stays on counts.
    if settings.get('selection_distance'):
        parts.append(f"Correction decision: {esc(settings['selection_distance'])}.")
    # The numbers that decision used (relative-abundance Bray), separate from the count-based table below.
    if q.get('baseline_p', 'NA') not in ('NA', '', None):
        after = (f"{fmt(value(q.get('final_R2')))} after correction" if state == 'corrected' else 'not corrected')
        parts.append(f"Decision metric (Bray–Curtis on relative abundances): batch R2 {fmt(value(q.get('baseline_R2')))} "
                     f"(p {fmt(value(q.get('baseline_p')), 3)}) before; {after}.")
    # Depth by batch, so the user can decide whether to add depth as a covariate.
    if settings.get('depth_by_batch'):
        added = settings.get('depth_covariate', '').startswith('added')
        parts.append(f"Sequencing depth (median total count by batch): {esc(settings['depth_by_batch'])}. "
                     + ("Depth (log10 total count) was added as a covariate." if added else DEPTH_NOTE))
    conf = q.get('confounding_note', 'NA')
    if conf not in ('NA', None):
        parts.append(f"Review needed: {esc(conf)} (provisional threshold 0.8).")
    note = q.get('sample_type_note', 'NA')
    if note != 'NA':
        warn = f' {MIXED_WARNING}' if 'not a protected covariate' in note else ''
        parts.append(f"Note: {esc(note)}; the run continued as configured.{warn}")
    return ' '.join(parts)


def term_table(rows):
    patient = any(r.get('pp_before') is not None or r.get('pp_after') is not None for r in rows)
    cols = ['Term', 'R² before', 'R² after', 'Change', 'Non-batch share (before → after)', 'p before', 'p after'] + (
        ['p before (within patient)', 'p after (within patient)'] if patient else []) + ['Look', 'Note']
    head = ''.join(f'<th>{h}</th>' for h in cols)
    def cells(r):
        c = [esc(r['term'] + (' (batch)' if r['role'] == 'batch' else '')), fmt(r['r2_before']), fmt(r['r2_after']),
             fmt(r['change'], 4) if r['change'] is not None else 'NA',
             (f"{fmt(r['share_before'])} → {fmt(r['share_after'])}" if r.get('share_before') is not None else ''),
             fmt(r['p_before'], 3), fmt(r['p_after'], 3)]
        if patient:
            c += [fmt(r.get('pp_before'), 3), fmt(r.get('pp_after'), 3)]
        return c + ['Review needed' if r['review'] else '', esc(r['note'])]
    body = ''.join('<tr>' + ''.join(f'<td>{x}</td>' for x in cells(r)) + '</tr>' for r in rows)
    return f'<div class="table-wrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def effect_summary(rows):
    """Sample-type comparison before/after. Taxa with |log2 ratio| >= 1 before correction count as strong."""
    if not rows or rows[0]['status'] != 'computed':
        return None, (rows[0]['status'] if rows else 'not computed: no effect-size file')
    a, b, col = rows[0]['group_a'], rows[0]['group_b'], rows[0]['group_column']
    pairs = [(r['taxon'], value(r['log2_ratio_before']), value(r['log2_ratio_after'])) for r in rows]
    strong = [(t, x, y) for t, x, y in pairs if x is not None and y is not None and abs(x) >= STRONG_LOG2]
    if not strong:
        return {'a': a, 'b': b, 'col': col, 'strong': [], 'kept': None, 'flips': 0, 'review': False}, None
    kept = sum(abs(y) for _, _, y in strong) / sum(abs(x) for _, x, _ in strong)
    flips = sum((x > 0) != (y > 0) for _, x, y in strong)
    return {'a': a, 'b': b, 'col': col, 'strong': strong, 'kept': kept, 'flips': flips,
            'review': kept <= 1 - REVIEW_DROP or flips > 0}, None


def protection_line(rows):
    """For a user-named comparison: whether its column was among the covariates used for correction."""
    col = esc(rows[0]['group_column'])
    if rows[0].get('protected_covariate') == 'yes':
        return f'<p>{col} is a protected covariate (listed in --batch_corr_covariates).</p>'
    return (f'<p>{col} is not a protected covariate: correction can shrink this difference; '
            f'list it in --batch_corr_covariates to protect it.</p>')


def effect_block(rows, title='Sample-type comparison'):
    e, why = effect_summary(rows)
    named = title != 'Sample-type comparison'
    if e is None:
        return f'<h3>{title}' + (f' ({esc(rows[0]["group_column"])})' if named else '') + f'</h3><p>{esc(why)}.</p>'
    extra = protection_line(rows) if named else ''
    head = (f'<h3>{title}: {esc(e["a"])} vs {esc(e["b"])} ({esc(e["col"])})</h3>') + extra
    if not e['strong']:
        return head + '<p>No taxon differed at least two-fold before correction; nothing to compare.</p>'
    look = ('Review needed: ' if e['review'] else 'No review flag: ') + (
        f'{len(e["strong"])} taxa differed at least two-fold before correction; after correction their summed '
        f'absolute log2 ratio is {100 * e["kept"]:.0f}% of the value before, and {e["flips"]} changed direction.')
    rule = ('“Review needed” here marks a fall of at least 20% in this summed size or any change of direction '
            '(provisional rule, not validated).')
    top = sorted(e['strong'], key=lambda t: -abs(t[1]))[:15]
    table = ('<div class="table-wrap"><table><thead><tr><th>Taxon</th><th>log2 ratio before</th><th>log2 ratio after</th>'
             '</tr></thead><tbody>' + ''.join(f'<tr><td>{esc(t)}</td><td>{fmt(x, 2)}</td><td>{fmt(y, 2)}</td></tr>'
                                             for t, x, y in top) + '</tbody></table></div>')
    return head + f'<p>{look} {rule}</p><details><summary>Largest differences before correction (up to 15 taxa)</summary>{table}</details>'


def patient_line(raw_counts):
    """Within-patient permutations, shown only when a patient was declared."""
    if not raw_counts or raw_counts[0].get('Patient_note', 'NA') in ('NA', 'no patient column declared'):
        return ''
    r = raw_counts[0]
    if r['Patient_note'] != 'computed':
        return f'<p>Within-patient permutation: {esc(r["Patient_note"])}.</p>'
    return (f'<p>Within-patient permutation: {esc(r["Patients"])} patients, {esc(r["Patients_multi"])} with more than one sample. '
            'For terms that vary within patients, the p-value is also computed with permutations restricted within patients; '
            'terms constant within patients show NA there. R² does not depend on the permutation scheme.</p>')


def dispersion_line(rows, after_source):
    """Betadisper on the same Bray distances, batch (and protected sample type), before and after."""
    out = []
    for grp in dict.fromkeys(r['Group'] for r in rows):
        b = next((r for r in rows if r['Group'] == grp and r['Source'] == 'raw'), None)
        a = next((r for r in rows if r['Group'] == grp and r['Source'] == after_source), None)
        if not b or not a or b['Status'] != 'computed' or a['Status'] != 'computed':
            why = next((r['Status'] for r in (b, a) if r and r['Status'] != 'computed'), 'not computed')
            out.append(f'{esc(grp)}: {esc(why)}')
            continue
        out.append(f'{esc(grp)}: F {fmt(value(b["F"]), 2)} → {fmt(value(a["F"]), 2)}, p {fmt(value(b["p"]), 3)} → {fmt(value(a["p"]), 3)}')
    if not out:
        return ''
    return ('<p>Dispersion check (betadisper, Bray): ' + '; '.join(out) + '. A small p means the groups differ in spread, '
            'so the PERMANOVA R² above may partly reflect spread rather than location.</p>')


def render_bundle(path):
    m, data, settings = load(path)
    quality, status, perm = data['batch_correction_quality.tsv'], data['BATCH_CORRECTION_STATUS.tsv'], data['permanova_summary.tsv']
    batch_column = settings.get('batch_column', 'NA')
    covariates = [c for c in settings.get('covariates', '').split(',') if c]
    rows, after_source = term_rows(perm, batch_column, covariates)
    counts = [r for r in perm if r['Transform'] == 'counts' and r['Source'] == 'raw']
    cond = (f"Same {esc(counts[0]['N'])} samples, {esc(counts[0]['Distance'])} distance on counts, "
            f"{esc(counts[0]['Evaluation'])} evaluation, model {esc(counts[0]['Model'])}.") if counts else 'No PERMANOVA rows.'
    transforms = sorted({r['Transform'] for r in perm} - {'counts'})
    extra = ''.join(f'<h4>{esc(t)}</h4>' + term_table(term_rows(perm, batch_column, covariates, t)[0]) for t in transforms)
    # Rows of a user-named comparison carry protected_covariate (yes/no); older bundles have sample-type rows only.
    all_effects = data['batch_correction_effect_sizes.tsv']
    effects = [r for r in all_effects if r.get('protected_covariate') in (None, 'NA')]
    named = [r for r in all_effects if r.get('protected_covariate') not in (None, 'NA')]
    computed = bool(effects) and effects[0]['status'] == 'computed'
    named_computed = bool(named) and named[0]['status'] == 'computed'
    raw_counts = [r for r in perm if r['Transform'] == 'counts' and r['Source'] == 'raw']
    pnote = raw_counts[0].get('Patient_note', 'NA') if raw_counts else 'NA'
    patient_done = pnote == 'computed'
    not_assessed = [x for x in NOT_ASSESSED if not ((computed or named) and x == NAMED_COMPARISON)
                    and not (patient_done and x == 'Patient-restricted permutations.')]
    if named_computed:
        not_assessed.append(f'Clinical comparisons other than {"sample type and " if computed else ""}{named[0]["group_column"]}.')
    elif computed:
        not_assessed.append('Clinical comparisons other than sample type.')
    flagged = [r['term'] for r in rows if r['review']]
    look = ('Review needed: ' + ', '.join(esc(t) for t in flagged) + '.') if flagged else 'No covariate met the review rule.'
    return (f'<section><h2>Batch correction — {esc(m["rank"])}</h2>'
            f'<p>{summary_line(m, quality, status, settings)}</p>'
            f'<p>{cond} After-correction table: {esc(after_source or "none")}. {look} '
            f'“Review needed” marks a covariate whose share of the non-batch variation (R² ÷ (1 − batch R²)) fell by '
            f'at least 20% (if p < 0.05 before) or rose by at least 20% and 0.01 (if p < 0.05 after), or whose p-value '
            f'moved from below 0.05 to 0.05 or above '
            f'(provisional rule, not validated). Removing batch variation makes every other R² grow by itself, so the '
            f'share is what is compared. It is a pointer, not a pass or fail.</p>'
            + term_table(rows) + patient_line(raw_counts) + dispersion_line(data['batch_correction_dispersion.tsv'], after_source) +
            (f'<details><summary>Same comparison after each normalization (sensitivity; RCLR with Euclidean distance is not RPCA)</summary>{extra}</details>' if extra else '') +
            effect_block(effects) + (effect_block(named, 'Named comparison') if named else '') +
            '<h3>Not assessed in this run</h3><ul>' + ''.join(f'<li>{esc(x)}</li>' for x in not_assessed) + '</ul></section>')


def render(paths):
    if not paths:
        return ''
    return '\n'.join(render_bundle(p) for p in paths)
