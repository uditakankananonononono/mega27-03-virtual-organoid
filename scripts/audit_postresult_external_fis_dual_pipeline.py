"""Compare published demo outputs from two pipelines without treating tracks as truth."""
import collections,csv,hashlib,json
from pathlib import Path
root=Path('/tmp/vo3_external_fis'); prior=json.loads(Path('results/postresult_external_fis_demo_source.json').read_text())
files=[('CellProfiler',root/'FIS_rawdata.csv','TrackObjects_Label_4',prior['summary_csv_sha256']),('Fiji/ImageJ',root/'FIS_rawdata_ij.csv','TrackObjects_Label','072554a0bf397adcaa53b83670626d8ed8e4e01d6617911534d0c5dc5bf10194')]
out={'status':'post-result source schema audit of two pipelines on same demo images','protocol':'notes/postresult_external_fis_dual_pipeline_protocol.md','source_url':prior['source_url'],'pipelines':{},'limits':'Same source images and one plate, not independent cohorts. Duplicate well/time/label groups in both pipelines; no ground-truth tracking. Original VO3 sign/specificity FAIL and benchmark win withdrawn.'}
for name,p,col,sha in files:
 assert p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==sha
 meta={};frames=collections.defaultdict(set);labeltimes=collections.defaultdict(set);seen=collections.Counter();n=0;missing=0;example=None
 for r in csv.DictReader(p.open(newline='')):
  n+=1;w=r['Metadata_wellNum'];t=int(r['Metadata_timeNum']);lab=r[col];a=r['Metadata_compound'];dose=r['Metadata_concentration'];area=float(r['Math_area_micronsq']);assert area>0 and t in range(7)
  if w in meta: assert meta[w]==(a,dose)
  else:meta[w]=(a,dose)
  frames[w].add(t)
  if lab in ('NA','NaN','nan',''):
   missing+=1;continue
  labeltimes[w,lab].add(t);key=(w,t,lab);seen[key]+=1
  if seen[key]==2 and example is None: example={'well':w,'time':t,'label':lab}
 assert len(meta)==64 and all(f==set(range(7)) for f in frames.values())
 arms=collections.Counter(a for a,dose in meta.values());doses=collections.Counter(dose for a,dose in meta.values())
 assert set(arms)=={'fsk','fsk_809','fsk_770','fsk_770_809'} and '0' not in doses
 out['pipelines'][name]={'sha256':sha,'rows':n,'wells':len(meta),'times_per_well':7,'missing_label_rows':missing,'n_distinct_within_well_labels':len(labeltimes),'n_labels_with_baseline_and_final':sum(0 in ts and 6 in ts for ts in labeltimes.values()),'n_duplicate_well_time_label_keys':sum(v>1 for v in seen.values()),'duplicate_excess_object_rows':sum(v-1 for v in seen.values() if v>1),'max_objects_per_well_time_label':max(seen.values()),'example_duplicate_key':example,'arm_wells':dict(sorted(arms.items())),'dose_wells':dict(sorted(doses.items(),key=lambda x:float(x[0])))}
assert out['pipelines']['CellProfiler']['n_duplicate_well_time_label_keys']==prior['n_duplicated_well_time_label_keys']==929
assert out['pipelines']['CellProfiler']['duplicate_excess_object_rows']==prior['n_duplicate_well_time_label_excess_rows']==965
Path('results/postresult_external_fis_dual_pipeline.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
for k,d in out['pipelines'].items():print(k,{z:d[z] for z in ('rows','missing_label_rows','n_distinct_within_well_labels','n_labels_with_baseline_and_final','n_duplicate_well_time_label_keys','duplicate_excess_object_rows')})
