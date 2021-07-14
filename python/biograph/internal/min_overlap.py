import sys
import gzip
all_event_count = 0
min_overlap_removed = 0

filtered = open("no_min_overlap80_disc.vcf", 'w')

with gzip.GzipFile(sys.argv[1], 'r') as fh:
    for line in fh:
        if line.startswith("#"):
            filtered.write(line)
            sys.stdout.write(line)
            continue
        data = line.strip().split('\t')
        fmt = data[8]
        fmt_idx = fmt.split(':').index("OV")
        has_min = False
        all_event_count += 1
        for samp in data[9:]:
            try:
                samp_ov = int(samp.split(':')[fmt_idx])
            except Exception:
                continue
            if samp_ov >= 80:
                has_min = True
                break
        if has_min: 
            sys.stdout.write(line)
        else:
            min_overlap_removed += 1
            filtered.write(line)
sys.stderr.write("Removed %d of %d variants %.2f\n" % (min_overlap_removed, \
                all_event_count, min_overlap_removed/float(all_event_count)))

                

