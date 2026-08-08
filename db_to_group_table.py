#!/usr/bin/env python3
"""
Bridge TurtleWave/wonambi per-subject output -> group-analysis tables.

For each subject folder under DATA (containing wonambi/neural_events.db and a
wonambi/*_clean_rebuilt.xml staging file), this:
  - reads stage durations from the XML hypnogram (sum of epoch lengths per stage),
  - queries neural_events.db for per-channel COUNT, mean peak2peak amplitude and
    mean duration of each event x frequency band,
  - computes density = count / NREM(2+3) minutes,
and writes the canonical group-analysis layout to OUT:
  eventStat_<param>_<event_band>.csv   rows=subjects, cols=channels  (param in density/amplitude/duration)
  subjects.csv                         Subject, group, session
  channels.csv                         channel  (union across subjects)

No dependency on the parameter/density CSVs (DB + XML only), so it survives the
planned TurtleWave outlier-removal / CSV-removal. Uses only the Python stdlib.

Usage:  python3 db_to_group_table.py <DATA_dir> <OUT_dir> [--group LABEL]

--group is the cohort label written into subjects.csv. It defaults to 'group1',
which is only ever correct for a single-cohort run -- pass it explicitly when you
are building a table with more than one group, or every subject lands in the same
level and the group contrast has nothing to test.
"""
import os, sys, re, csv, sqlite3, glob, argparse
import xml.etree.ElementTree as ET
from collections import defaultdict

EVENT_BANDS = [   # (event_type, freq_lower, freq_upper, label)
    ('spindle',   9.0, 12.0, 'spindle_9-12'),
    ('spindle',  12.0, 15.0, 'spindle_12-15'),
    ('slow_wave', 0.5, 1.25, 'sw_0.5-1.25'),
]
PARAMS = ['density', 'amplitude', 'duration']
NREM_STAGES = {'NREM2', 'NREM3'}     # events are tagged stage 'NREM2NREM3'

def stage_minutes(xml_path, stages=NREM_STAGES):
    root = ET.parse(xml_path).getroot()
    secs = 0.0
    for ep in root.findall('.//epoch'):
        st = ep.find('stage')
        if st is not None and st.text in stages:
            secs += float(ep.find('epoch_end').text) - float(ep.find('epoch_start').text)
    return secs / 60.0

def subject_event_stats(db_path, nrem_min):
    """Return {label: {channel: {density, amplitude, duration}}}."""
    con = sqlite3.connect(db_path); cur = con.cursor()
    out = {}
    for etype, flo, fhi, label in EVENT_BANDS:
        cur.execute(
            "SELECT channel, COUNT(*), AVG(peak2peak_amp), AVG(duration) "
            "FROM events WHERE event_type=? AND freq_lower=? AND freq_upper=? "
            "GROUP BY channel", (etype, flo, fhi))
        d = {}
        for ch, n, amp, dur in cur.fetchall():
            d[ch] = {'density': n / nrem_min if nrem_min else float('nan'),
                     'amplitude': amp, 'duration': dur}
        out[label] = d
    con.close()
    return out

def find_subjects(data_dir):
    subs = []
    for name in sorted(os.listdir(data_dir)):
        sdir = os.path.join(data_dir, name, 'wonambi')
        db = os.path.join(sdir, 'neural_events.db')
        xmls = glob.glob(os.path.join(sdir, '*_clean_rebuilt.xml'))
        if os.path.isfile(db) and xmls:
            subs.append((name, db, xmls[0]))
    return subs

def natkey(ch):  # E1,E2,...,E10,...,Cz  natural order
    m = re.match(r'E(\d+)$', ch)
    return (0, int(m.group(1))) if m else (1, ch)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('data_dir'); ap.add_argument('out_dir')
    ap.add_argument('--group', default='group1',
                    help="cohort label written into subjects.csv (default: group1)")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    subs = find_subjects(a.data_dir)
    assert subs, f'No subjects (wonambi/neural_events.db + xml) under {a.data_dir}'
    print(f'Subjects: {[s[0] for s in subs]}')

    # stats[label][param][subject][channel] = value ; collect channel union
    stats = {lab: {p: {} for p in PARAMS} for *_ , lab in EVENT_BANDS}
    channels = set()
    meta = []
    for name, db, xml in subs:
        nrem = stage_minutes(xml)
        ev = subject_event_stats(db, nrem)
        n_tot = sum(len(ev[l]) for l in ev)
        print(f'  {name}: NREM2+3 = {nrem:.1f} min; channels with events = '
              + ', '.join(f'{l}:{len(ev[l])}' for l in ev))
        for lab in ev:
            for ch, vals in ev[lab].items():
                channels.add(ch)
                for p in PARAMS:
                    stats[lab][p].setdefault(name, {})[ch] = vals[p]
        sess = name.split('_')[-1] if '_' in name else ''
        meta.append((name, a.group, sess))

    chans = sorted(channels, key=natkey)
    subjects = [s[0] for s in subs]

    # write per (param, event-band) subject x channel matrices
    for *_, lab in EVENT_BANDS:
        for p in PARAMS:
            path = os.path.join(a.out_dir, f'eventStat_{p}_{lab}.csv')
            with open(path, 'w', newline='') as f:
                w = csv.writer(f); w.writerow(['Subject'] + chans)
                for s in subjects:
                    row = [s] + [stats[lab][p].get(s, {}).get(ch, '') for ch in chans]
                    w.writerow(row)
    # subjects.csv + channels.csv
    with open(os.path.join(a.out_dir, 'subjects.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['Subject', 'group', 'session']); w.writerows(meta)
    with open(os.path.join(a.out_dir, 'channels.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['channel']); w.writerows([[c] for c in chans])

    print(f'Wrote {len(EVENT_BANDS)*len(PARAMS)} matrices + subjects.csv + channels.csv '
          f'({len(subjects)} subjects x {len(chans)} channels) to {a.out_dir}')

if __name__ == '__main__':
    main()
