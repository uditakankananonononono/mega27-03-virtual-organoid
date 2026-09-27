"""External FIS demo source/track integrity only; no biological validation."""
import collections,csv,hashlib,json
from pathlib import Path
root=Path('/tmp/vo3_external_fis');p=root/'FIS_rawdata.csv';q=root/'demoplate_01.txt'
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
assert sha(p)=='1cb65d37fe78ad232c9829634c4ed7be8afc52e9e0308a519673eb3a9d91f807' and sha(q)=='70195928bdca463f287451bb3f056203a355893d8a6c37713c166461298d13f8'
rows=list(csv.DictReader(p.open(newline='')));assert len(rows)==31900
meta={};frames=collections.defaultdict(set);labels=collections.defaultdict(lambda:collections.defaultdict(list));NA=collections.Counter();dups=[];seen=set()
for r in rows:
 w=r['Metadata_wellNum'];t=int(r['Metadata_timeNum']);a=r['Metadata_compound'];dose=r['Metadata_concentration'];area=float(r['Math_area_micronsq']);lab=r['TrackObjects_Label_4']
 assert area>0 and t in range(7)
 frames[w].add(t)
 if w in meta:assert meta[w]==(a,dose)
 else:meta[w]=(a,dose)
 if lab in ('NA',''):
  NA[w]+=1;continue
 key=(w,t,lab)
 if key in seen:dups.append(key)
 seen.add(key)
 labels[w][lab].append((t,area))
assert len(meta)==64 and all(x==set(range(7)) for x in frames.values())
full=[];baseline=[];final=[];label_total=0
for w,d in labels.items():
 label_total+=len(d)
 for lab,vals in d.items():
  times={x for x,y in vals}
  if 0 in times:baseline.append((w,lab))
  if 6 in times:final.append((w,lab))
  if {0,6}<=times:full.append((w,lab))
arm_wells=collections.Counter(a for a,dose in meta.values());dose_wells=collections.Counter(dose for a,dose in meta.values())
assert set(arm_wells)=={'fsk','fsk_809','fsk_770','fsk_770_809'} and '0' not in dose_wells
out={'source_url':'https://github.com/hmbotelho/FIS_analysis','readme_url':'https://raw.githubusercontent.com/hmbotelho/FIS_analysis/master/README.md','summary_csv_url':'https://raw.githubusercontent.com/hmbotelho/FIS_analysis/master/demo_dataset/05-images_analysis/demoplate_01--cellprofiler--analysis/FIS_rawdata.csv','summary_csv_sha256':sha(p),'infile_url':'https://raw.githubusercontent.com/hmbotelho/FIS_analysis/master/demo_dataset/02-microscope_infile/demoplate_01.txt','infile_sha256':sha(q),'source_tree_snapshot':'GitHub API repos/hmbotelho/FIS_analysis/git/trees/master?recursive=1, inspected on 2026-09-27','n_raw_exported_tiff':448,'n_cellprofiler_mask_png':448,'n_object_csv':64,'n_wells':len(meta),'frames_per_well':7,'n_object_rows':len(rows),'n_na_track_rows':sum(NA.values()),'wells_with_na_track':sum(v>0 for v in NA.values()),'n_distinct_well_track_labels':label_total,'n_labels_with_baseline':len(baseline),'n_labels_with_endpoint':len(final),'n_labels_with_both':len(full),'n_duplicate_well_time_label_excess_rows':len(dups),'n_duplicated_well_time_label_keys':len(set(dups)),'max_rows_per_well_time_label':max(collections.Counter((w,t,lab) for w,d in labels.items() for lab,vals in d.items() for t,area in vals).values()),'arm_wells':dict(sorted(arm_wells.items())),'dose_wells':dict(sorted(dose_wells.items(),key=lambda z:float(z[0]))),'n_zero_forskolin_wells':0,'n_cftr_inhibitor_wells':0,'assay_limits':'TrackObjects_Label_4 is NOT one-to-one with well x time: multiple objects carry the same label in 929 well/time/label groups; cannot silently treat each label as one unique organoid track. One demonstration plate; README says homozygous for a class II CFTR mutation, but no independent donor IDs or multi-donor cohort. Fluorescent calcein time series, VX-809/VX-770, different from original VO3 brightfield VX445/VX661/VX770. Processed masks are not hand-annotated instance ground truth; track labels algorithm-derived. Not an independent phenotype validation or original benchmark gate.'}
Path('results/postresult_external_fis_demo_source.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print({k:out[k] for k in ('n_wells','n_object_rows','n_na_track_rows','n_distinct_well_track_labels','n_labels_with_baseline','n_labels_with_endpoint','n_labels_with_both','n_duplicate_well_time_label_excess_rows','n_duplicated_well_time_label_keys','arm_wells','dose_wells')})
