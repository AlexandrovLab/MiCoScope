"""Actual decontamination.R boundary regressions; not parameter-optimality tests."""
import csv,subprocess,tempfile,unittest,os,json,shutil
from contextlib import nullcontext
from pathlib import Path
from test_decontam_zero_total import decontam_available,_decontam_rscript,ROOT

@unittest.skipUnless(decontam_available(),"R/decontam/phyloseq/optparse/digest unavailable: full-script boundary tests NOT RUN")
class FilterBoundaries(unittest.TestCase):
 def test_selected_denominator_and_final_exclusion(self):
  saved = os.environ.get("CMP_FILTER_BOUNDARY_OUTPUT_DIR")
  if saved: Path(saved).mkdir(parents=True, exist_ok=True)
  with (nullcontext(saved) if saved else tempfile.TemporaryDirectory()) as tmp:
   d=Path(tmp); samples=[f'S{i}' for i in range(1,42)]
   def write(name,rows):
    with (d/name).open('w') as f:csv.writer(f,delimiter='\t',lineterminator='\n').writerows(rows)
   taxa={'Always':[20]*41,'Count4':[4]*40+[999],'Count5':[5]*40+[999],
         'Below':[5]+[0]*39+[999],'At':[5]+[0]*19+[5]+[0]*19+[999],
         'Above':[5,5]+[0]*18+[5]+[0]*19+[999]}
   write('otu.tsv',[['clade_name']+samples]+[[f'g__G|s__{k}']+v for k,v in taxa.items()])
   write('meta.tsv',[['sampleid','sample_type','shipment_batch']]+[[s,'Other' if i==40 else 'Tumor' if i%20<10 else 'Blood',f'B{1+i//20}'] for i,s in enumerate(samples)])
   write('lib.tsv',[['sample','bracken_genus_total','bracken_species_total']]+[[s,5000,5000] for s in samples])
   for threshold,expected in [('.05',{'Always','Count5','At','Above'}),('.0500001',{'Always','Count5','Above'})]:
    prefix='edge'+threshold.replace('.','_')
    case=d/prefix; case.mkdir()
    for name in ['otu.tsv','meta.tsv','lib.tsv']: shutil.copy2(d/name,case/name)
    cmd=[_decontam_rscript(),str(ROOT/'scripts/decontamination.R'),'--otu_table','otu.tsv','--metadata','meta.tsv',
     '--prefix',prefix,'--start_mode','decontam','--environment_spec',str(ROOT/'conda_envs/decontam_env.yml'),
     '--nextflow_version','boundary-test','--taxon_column','clade_name','--taxonomic_rank','species',
     '--sample_id_column','sampleid','--batch_column','shipment_batch','--type_column','sample_type',
     '--control_role','surrogate_biological','--positive_values','Tumor','--control_values','Blood',
     '--unselected_type_policy','exclude','--threshold','.1','--min_prevalence',threshold,'--min_abundance','5',
     '--min_batches','2','--min_total_per_batch','5','--min_positive_per_batch','1','--min_control_per_batch','1',
     '--library_table','lib.tsv','--remove_taxa','']
    r=subprocess.run(cmd,cwd=case,capture_output=True,text=True)
    (case/f"{prefix}.command.json").write_text(json.dumps(cmd,indent=2))
    (case/f"{prefix}.stdout.log").write_text(r.stdout)
    (case/f"{prefix}.stderr.log").write_text(r.stderr)
    self.assertEqual(r.returncode,0,r.stderr)
    with (case/f'{prefix}.decontam_pkg_decontaminated.csv').open() as f:
     reader=csv.reader(f);header=next(reader);rows=list(reader)
    self.assertEqual(len(header)-1,40)
    self.assertEqual({row[0].split('s__')[-1] for row in rows},expected)
    with (case/f'{prefix}.filter_report.tsv').open() as f:ledger={r['taxon_id'].split('s__')[-1]:r for r in csv.DictReader(f,delimiter='\t')}
    self.assertEqual(float(ledger['Count4']['after_abundance_sum']),0)
    self.assertEqual(float(ledger['At']['after_abundance_prevalence']),.05)
    self.assertEqual(ledger['Below']['terminal_reason'],'below_min_prevalence_after_abundance_filter')
if __name__=='__main__':unittest.main()
