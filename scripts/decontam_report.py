"""Render explicit current-run decontamination evidence; no results-tree discovery."""
import csv
import hashlib
import html
import math
from pathlib import Path

STAGES = ['input', 'declared_removal', 'selected_cohort', 'abundance', 'prevalence', 'final']
LABELS = ['Input table', 'After declared taxon removal', 'Selected cohort',
          'After count filter', 'After prevalence filter', 'After decontam']
REQUIRED = ['stage_flow.tsv', 'taxon_ledger.tsv', 'batch_type_summary.tsv', 'taxon_detection.tsv',
            'sample_validation.tsv', 'batch_eligibility.tsv', 'contaminants_by_batch.tsv',
            'contaminants_final.tsv', 'filter_report.tsv', 'decontaminated.csv']


def read(path, delimiter='\t'):
    with path.open() as f:
        return [{k: ('NA' if v == '' else v) for k,v in row.items()}
                for row in csv.DictReader(f, delimiter=delimiter)]


def number(x):
    n = float(x)
    if not math.isfinite(n) or n < 0:
        raise ValueError(f'Invalid nonnegative evidence count: {x}')
    return n


def load(path):
    path = Path(path)
    pairs = read(path / 'manifest.tsv')
    m = {r['key']: r['value'] for r in pairs}
    if len(m) != len(pairs) or m.get('report_schema') != '1' or m.get('status') != 'completed':
        raise ValueError('Invalid or incomplete decontamination evidence manifest')
    data = {}
    for name in REQUIRED:
        f = path / name
        if hashlib.sha256(f.read_bytes()).hexdigest() != m.get('report_sha256.' + name):
            raise ValueError(f'Decontamination evidence hash mismatch: {name}')
        data[name] = read(f, ',' if name.endswith('.csv') else '\t')
    flow = data['stage_flow.tsv']
    if [r['stage'] for r in flow] != STAGES:
        raise ValueError('Unexpected decontamination stages')
    for i, r in enumerate(flow):
        for k in ('samples', 'taxon_rows', 'detected_taxa', 'count_sum'):
            number(r[k])
        if number(r['detected_taxa']) > number(r['taxon_rows']):
            raise ValueError('Detected taxa exceed table rows')
        if i and number(flow[i-1]['count_sum']) - number(r['count_sum']) != number(r['lost_count']):
            raise ValueError('Stage count loss is inconsistent')
        if r['known_count_sum'] != 'NA' and number(r['known_count_sum']) > number(r['count_sum']):
            raise ValueError('Blacklist numerator exceeds denominator')
    ledger = data['taxon_ledger.tsv']
    if len({r['taxon_id'] for r in ledger}) != len(ledger):
        raise ValueError('Duplicate taxon ledger rows')
    for key, stage in [('input_all_samples',0), ('after_abundance',3), ('final_count',5)]:
        if sum(number(r[key]) for r in ledger) != number(flow[stage]['count_sum']):
            raise ValueError(f'Taxon ledger and stage counts disagree: {key}')
    sample = data['sample_validation.tsv']
    included = [r for r in sample if r['final_status'] == 'included']
    if len(included) != number(flow[2]['samples']) or len(included) != number(m['prevalence_denominator']):
        raise ValueError('Selected cohort denominator disagrees with sample ledger')
    if sum(number(r['final_total_count']) for r in included) != number(flow[5]['count_sum']):
        raise ValueError('Sample totals disagree with final matrix')
    final = data['decontaminated.csv']
    taxcol = m['taxon_column']
    if len(final) != number(flow[5]['taxon_rows']):
        raise ValueError('Final matrix taxon count disagrees')
    if sum(number(v) for row in final for k,v in row.items() if k != taxcol) != number(flow[5]['count_sum']):
        raise ValueError('Final matrix counts disagree')
    matrix_totals = {r[taxcol]: sum(number(v) for k,v in r.items() if k != taxcol) for r in final}
    if len(matrix_totals) != len(final):
        raise ValueError('Duplicate taxa in final matrix')
    for r in ledger:
        if number(r['final_count']) != matrix_totals.get(r['taxon_id'], 0):
            raise ValueError('Taxon final count disagrees with final matrix')
        if (r['final_action'] == 'retain') != (r['taxon_id'] in matrix_totals):
            raise ValueError('Taxon action disagrees with final matrix')
    reasons = {r['taxon_id']: r for r in data['filter_report.tsv']}
    for r in ledger:
        if r['taxon_id'] in reasons:
            other = reasons[r['taxon_id']]
            if any(r[k] != other[k] for k in ['terminal_reason', 'final_action', 'known_contaminant']):
                raise ValueError('Taxon reason/flag disagrees with filter report')
        elif r['terminal_reason'] != 'declared_taxon_removal':
            raise ValueError('Missing taxon removal reason')
    for stage in STAGES[2:]:
        groups = [r for r in data['batch_type_summary.tsv'] if r['stage'] == stage]
        f = flow[STAGES.index(stage)]
        if sum(number(r['count_sum']) for r in groups) != number(f['count_sum']) or sum(number(r['samples']) for r in groups) != number(f['samples']):
            raise ValueError('Batch/type totals disagree with stage totals')
    return m, data


def esc(value):
    return html.escape(str(value))


def fmt(value):
    if str(value) == 'NA':
        return 'Not applicable'
    return f'{float(value):,.0f}'


def grid(rows, columns):
    return ('<div class="table-wrap" style="overflow:auto;max-height:480px"><table><thead><tr>' +
            ''.join('<th>'+esc(label)+'</th>' for key,label in columns) + '</tr></thead><tbody>' +
            ''.join('<tr>'+''.join('<td>'+esc(r.get(key, 'Not available'))+'</td>' for key,label in columns)+'</tr>' for r in rows) +
            '</tbody></table></div>' if rows else '<p>No rows at this stage.</p>')


def details(title, rows, columns, note=''):
    return '<details><summary>'+esc(title)+'</summary><p>'+esc(note)+'</p>'+grid(rows, columns)+'</details>'


def chart(title, labels, values, suffix='', maximum=None):
    vals = [float(x) for x in values]
    scale = maximum if maximum is not None else max(vals + [1])
    scale = scale or 1
    h = 44 * len(labels) + 25
    body = []
    for i,(label,val) in enumerate(zip(labels,vals)):
        y = 30 + 44*i
        body.append(f'<text x="5" y="{y}" font-size="14">{esc(label)}</text>'
                    f'<rect x="260" y="{y-16}" width="{320*val/scale:.2f}" height="22" fill="#286d8c"/>'
                    f'<text x="590" y="{y}" font-size="14">{val:,.1f}{esc(suffix)}</text>')
    return f'<h3>{esc(title)}</h3><svg role="img" aria-label="{esc(title)}" viewBox="0 0 760 {h}" style="width:100%;max-width:960px">'+''.join(body)+'</svg>'


def render_bundle(path, qc_rows=None):
    m,d = load(path)
    if qc_rows is not None:
        qc = {r['sample']:r for r in qc_rows}
        for r in d['sample_validation.tsv']:
            if r['in_otu'] != 'TRUE':
                continue  # Metadata-only records do not enter the QC sample cohort.
            q = qc.get(r['sample_id'])
            fields = {'sample_type':'decontam_sample_type', 'role':'decontam_role',
                      'final_status':'decontam_final_status', 'result_status':'decontam_result_status',
                      'result_reason':'decontam_result_reason'}
            if q is None or any((q.get(v) or 'NA') != r[k] for k,v in fields.items()):
                raise ValueError('QC summary disagrees with decontamination sample evidence')
    flow, ledger = d['stage_flow.tsv'], d['taxon_ledger.tsv']
    final = d['contaminants_final.tsv']
    eligible = d['batch_eligibility.tsv']
    tests = d['contaminants_by_batch.tsv']
    removed = sum(r['final_contaminant']=='TRUE' for r in final)
    callable_n = sum(r['callable']=='TRUE' for r in final)
    excluded_n = sum(r['final_status']=='excluded' and r['in_otu']=='TRUE' for r in d['sample_validation.tsv'])
    metadata_only_n = sum(r['reason']=='metadata_only' for r in d['sample_validation.tsv'])
    title = f"Decontamination evidence — {m['prefix']} ({m['taxonomic_rank']})"
    out = ['<section class="decontam-evidence" id="decontamination-evidence"><h2>'+esc(title)+'</h2>']
    out += [f'<p><b>{fmt(flow[-1]["detected_taxa"])} taxa retained; {removed} removed by decontam.</b> '
            f'{fmt(flow[2]["samples"])} selected samples; {excluded_n} input samples excluded; {metadata_only_n} metadata-only records excluded. '
            f'Result: <b>{esc(m["result_status"])}</b> ({esc(m["result_reason"])}).</p>',
            '<p>Other losses below are filtering or cohort selection, not decontam contaminant calls. '
            'Retained taxa are not certified contamination-free. A valid-empty output is not a biological negative.</p>']
    out += [f'<p><b>Input:</b> {esc(m["input_source"])} — {esc(m["input_filename"])}. '
            f'<b>Units:</b> {esc(m["count_unit"])}.</p>']
    if m['input_source'] == 'consensus_or_bracken_passthrough_species':
        out += ['<p>The Consensus section above identifies whether MetaPhlAn consensus or Bracken pass-through actually ran.</p>']
    ranknote = ('Species-level decontamination. Any downstream genus table is summed from retained species; this panel is not a direct genus-level test.'
                if m['taxonomic_rank']=='species' else 'Direct genus-level decontamination, not a sum of species-level decisions.')
    out += ['<p>'+ranknote+'</p>', '<h3>Applied rules and controls</h3>',
            f'<p>Counts &lt; <b>{esc(m["min_abundance"])}</b> become zero. Prevalence ≥ <b>{100*float(m["min_prevalence"]):g}%</b> '
            f'is required among <b>{esc(m["prevalence_denominator"])}</b> selected samples (positives and controls, including retained zero-count samples). '
            f'A taxon is removed when its decontam prevalence score is below <b>{esc(m["threshold"])}</b> in at least <b>{esc(m["min_batches"])}</b> eligible batches.</p>',
            f'<p><b>Declared control role:</b> {esc(m["control_role"])}. Positive labels: <b>{esc(m["positive_values"])}</b>; '
            f'control labels: <b>{esc(m["control_values"])}</b> (exact metadata labels). '
            'Biological surrogate controls are not experimental blanks; comparison with blood/normal tissue alone does not establish reagent contamination.</p>',
            '<p>Prevalence compares detection in controls versus positives. Requiring recurrence across batches limits single-batch calls but can miss contamination confined to one batch. '
            'These are configured rules, not thresholds calibrated for this run. Scores are not FDR-adjusted probabilities of contamination.</p>']
    out += [chart('A. How much of the input count remains?', LABELS, [r['count_sum'] for r in flow]),
            chart('B. How many taxa still have nonzero counts?', LABELS, [r['detected_taxa'] for r in flow]),
            '<p>Each bar shows the remaining total at that step. The selected-cohort step removes excluded samples; '
            'later losses use the same selected sample denominator. Taxa with zero counts are excluded from the detected-taxon bars, even if their rows remain in the table.</p>',
            details('Stage totals and losses',flow,[('stage','Stage'),('samples','Samples'),('taxon_rows','Table rows'),('detected_taxa','Nonzero taxa'),('count_sum','Count sum'),('lost_count','Loss from previous stage'),('lost_detected_taxa','Nonzero taxa lost')])]
    out += ['<h3>C. Which batches could support a decision?</h3>',
            f'<p>{sum(r["eligible"]=="TRUE" for r in eligible)} / {len(eligible)} selected-cohort batches eligible; '
            f'{callable_n} / {len(final)} evaluated taxa testable in at least {esc(m["min_batches"])} batches. '
            f'Per batch, at least {esc(m["min_total_per_batch"])} samples, {esc(m["min_positive_per_batch"])} positives and '
            f'{esc(m["min_control_per_batch"])} controls are required. Ineligible, absent and unevaluable results are distinct; none means “no contamination”.</p>',
            grid(eligible,[('batch','Batch'),('n_positive','Positive'),('n_control','Control'),('n_total','Total'),('eligible','Eligible'),('reason','Reason')])]
    diagnostic = []
    for b in eligible:
        ts = [r for r in tests if r['batch']==b['batch']]
        diagnostic.append(dict(batch=b['batch'], evaluated=sum(r['evaluation_status']=='evaluated' for r in ts),
          absent=sum(r['evaluation_status']=='absent_in_batch' for r in ts), unevaluable=sum(r['evaluation_status']=='unevaluable' for r in ts),
          calls=sum(r['contaminant']=='TRUE' for r in ts), status=('not tested — no taxa after filters' if m['analysis_status']=='not_run_no_taxa_after_filters' else 'eligible' if b['eligible']=='TRUE' else 'not tested — ineligible')))
    out += [grid(diagnostic,[('batch','Batch'),('status','Status'),('evaluated','Taxa tested'),('absent','Absent taxa'),('unevaluable','Unevaluable taxa'),('calls','Calls before recurrence rule')])]
    scores=[float(r['p']) for r in tests if r['evaluation_status']=='evaluated' and r['p'] not in ('NA','')]
    bounds=[0,.01,.05,.1,.5,1.00000001]
    bins=[sum(a<=v<b for v in scores) for a,b in zip(bounds,bounds[1:])]
    out += [chart('D. Distribution of available decontam scores', ['0–<0.01','0.01–<0.05','0.05–<0.10','0.10–<0.50','0.50–1.00'],bins),
            f'<p>{len(scores)} available scores from evaluated taxon–batch tests. These fixed display bins do not change the configured threshold. '
            'This distribution is descriptive: it does not measure false-discovery rate or stability.</p>']
    for r in d['batch_type_summary.tsv']:
        r['cell_detection'] = f'{r["positive_cells"]} / {r["possible_cells"]}' if number(r['possible_cells']) else 'Not applicable (no taxa)'
        r['sample_detection'] = f'{r["samples_with_counts"]} / {r["samples"]}'
    out += [details('Batch × sample type: counts and detection at each stage',d['batch_type_summary.tsv'],
            [('stage','Stage'),('batch','Batch'),('sample_type','Sample type'),('role','Role'),('eligible','Eligible'),('samples','Samples'),('count_sum','Count sum'),('detected_taxa','Nonzero taxa'),('cell_detection','Positive taxon/sample cells / possible cells'),('sample_detection','Samples with any count / samples')],
            'Selected cohort only. Cell detection uses the taxon rows retained at each stage; its denominator changes after prevalence/decontam. Compare count totals separately from detection frequency.'),
            details('Taxon × batch × sample type: evidence entering decontam',d['taxon_detection.tsv'],
            [('taxon_id','Taxon'),('batch','Batch'),('sample_type','Sample type'),('role','Role'),('eligible','Eligible batch'),('detected_samples','Detected samples'),('samples','Sample denominator'),('count_sum','Count sum')], 'After abundance and prevalence filters; detection is count > 0. Ineligible batches are shown but are not tested.')]
    out += ['<h3>E. Blacklist overlap is a flag, not a removal rule</h3>',
            f'<p><b>Status:</b> {esc(m["blacklist_status"])}. File: {esc(m["blacklist_filename"])}; '
            f'{esc(m["blacklist_species"])} normalized species names. Content version (SHA-256): <code style="overflow-wrap:anywhere">{esc(m["blacklist_sha256"])}</code>.</p>',
            '<p>'+esc(m['blacklist_matching'])+'. A listed species may remain; an unlisted species may be removed. Overlap is not proof of contamination.</p>']
    if m['blacklist_status']=='applied_as_flag_only':
        shares = []
        for i,r in enumerate(flow):
            denominator=number(r['count_sum']); numerator=number(r['known_count_sum'])
            shares.append(dict(stage=LABELS[i], numerator=fmt(numerator),denominator=fmt(denominator),
                share=f'{100*numerator/denominator:.2f}%' if denominator else 'Not applicable (zero count sum)'))
        out += [grid(shares,[('stage','Computation stage'),('numerator','Listed-species count'),('denominator','All-taxon count'),('share','Pooled share')]),
                '<p>Each share divides the listed-species count sum by the all-taxon count sum at that exact stage, not the average of sample percentages.</p>']
        nonzero = [(LABELS[i], 100*number(r['known_count_sum'])/number(r['count_sum'])) for i,r in enumerate(flow) if number(r['count_sum'])]
        out += [chart('Listed-species share of the count total', [r[0] for r in nonzero], [r[1] for r in nonzero], '%', maximum=100)]
    else:
        out += ['<p>Blacklist fractions are not computed here. Species lists cannot label an entire genus; a missing/empty list is not a zero-overlap result.</p>']
    out += [details('Complete taxon decisions and removal reasons',ledger,
            [('taxon_id','Taxon'),('input_all_samples','Input count: all samples'),('input_selected_samples','Input count: selected samples'),('after_abundance','After count filter'),('final_count','Final count'),('known_contaminant','Blacklist flag'),('final_action','Decision'),('terminal_reason','Reason'),('n_evaluable_batches','Evaluable batches'),('detected_batch_count','Batches called')],
            'Declared removals are separate from decontam. retained_not_callable means insufficient evaluable batches; retained_noncontaminant means not called by this rule, not proven clean.'),
            details('Sample inclusion, exclusions and original blacklist shares',d['sample_validation.tsv'],
            [('sample_id','Sample'),('sample_type','Sample type'),('role','Role'),('final_status','Status'),('reason','Reason'),('final_total_count','Final count'),('result_status','Result'),('result_reason','Result reason'),('known_contaminant_reads','Listed count before abundance filter'),('pct_known_contaminant','Listed count %')],m['sample_blacklist_denominator']),
            details('Per-taxon decontam scores and batch decisions',tests,
            [('taxon_id','Taxon'),('batch','Batch'),('n_positive','Positives'),('n_control','Controls'),('p','Score'),('evaluation_status','Test status'),('contaminant','Call')]),
            '<h3>What has not been assessed</h3><p>Setting sensitivity, repeated depth-matched resampling, alternative-mode comparison (including isNotContaminant), '
            'and false-discovery calibration are not assessed in this run. No claim of decision stability is made.</p>',
            details('All settings, software versions and evidence hashes',[dict(setting=k,value=v) for k,v in m.items()], [('setting','Setting / source'),('value','Value')]), '</section>']
    return '\n'.join(out)


def render(paths, qc_rows=None):
    if not paths:
        return '<section><h2>Decontamination evidence</h2><p>No current-run evidence bundle was supplied. This is not a zero-contamination result; legacy reports need an explicit compatible bundle.</p></section>'
    return '\n'.join(render_bundle(path, qc_rows) for path in paths)
