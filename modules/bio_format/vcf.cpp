
#include "modules/bio_format/vcf.h"
#include "modules/io/log.h"
#include "modules/bio_base/dna_base_set.h"
#include "modules/io/utils.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/io/config.h"
#include "modules/bio_base/struct_var.h"
#include "modules/bio_format/struct_var.h"
#include "modules/io/registry.h"
#include "modules/io/file_io.h"
#include "modules/io/version.h"
#include "modules/web/urlencode.h"
#include "tools/build_stamp.h"

#include <time.h>
#include <tr1/array>
#include <boost/foreach.hpp>
#include <numeric>

static const char* VCFversion = "VCFv4.1";
REGISTER_3(exporter, vcf, writable&, bool, const std::string&);

char* todaysDate() {
  // format like this: 20110705
  static char date[10];
  time_t rawtime;
  struct tm* timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);

  strftime(date, 10, "%Y%m%d", timeinfo);
  return date;
}

#define HEADER(hdr, id, num, type, desc)                            \
  m_sink.print("##" #hdr "=<ID=" #id ",Number=" #num ",Type=" #type \
               ",Description=\"" #desc "\">\n");
#define FILTER(id, desc) \
  m_sink.print("##FILTER=<ID=" #id ",Description=\"" #desc "\">\n");
#define ALT(id, desc) \
  m_sink.print("##ALT=<ID=" #id ",Description=\"" #desc "\">\n");
#define INFO(id, num, type, desc) HEADER(INFO, id, num, type, desc)
#define FORMAT(id, num, type, desc) HEADER(FORMAT, id, num, type, desc)

void printInfoString(writable& sink, const char* id,
                     const std::string& description) {
  std::string desc = urldecode(description);
  for (size_t i = 0; i < desc.size(); i++) {
    // mutate ';', ',', ' ', '\n' and  '\t' chars into ' ' to avoid conflicts
    // with the VCF format field separators
    if (desc[i] == ';' || desc[i] == ',' || desc[i] == '\t' ||
        desc[i] == '\n' || desc[i] == ' ') {
      desc[i] = '_';
    }
  }
  sink.print(";%s=%s", id, desc.c_str());
}

vcf_exporter::vcf_exporter(
    writable& byte_sink, const std::string& ref_name,
    const std::map<std::string, std::string>& extra_headers, bool use_events,
    const std::string& sample_name)
    : exporter(byte_sink),
      m_use_events(use_events),
      m_reference(ref_name),
      m_reference_assembly(m_reference.get_assembly()),
      m_extra_headers(extra_headers),
      m_sample_name(sample_name) {}

vcf_exporter::vcf_exporter(writable& byte_sink, const std::string& ref_name,
                           bool use_events, const std::string& sample_name)
    : exporter(byte_sink),
      m_use_events(use_events),
      m_reference(ref_name),
      m_reference_assembly(m_reference.get_assembly()),
      m_sample_name(sample_name) {}

vcf_exporter::vcf_exporter(writable& byte_sink, bool /*unused*/,
                           const std::string& serialized_exporter_data,
                           const std::string& sample_name)
    : exporter(byte_sink),
      m_use_events(true),
      m_reference(""),
      m_reference_assembly(m_reference.get_assembly()),
      m_sample_name(sample_name) {
  // throw io_exception("Not implemented");
}

void vcf_exporter::write_header() {
  SPLOG("vcf_exporter::write_header> Exporting variants");
  //_____ VCF Header ________________________________________________
  m_sink.print("##fileformat=%s\n", VCFversion);
  m_sink.print("##fileDate=%s\n", todaysDate());
  m_sink.print("##BioGraph.Variants=\"source=\"Spiral Genetics BioGraph\"");
  m_sink.print(",version=\"%s\"", biograph_current_version.make_string().c_str()); 
  m_sink.print(",build-revision=\"%s%s\"", get_build_scm_revision().c_str(),
               build_is_clean() ? "" : " (unclean workspace)");
  time_t build_timestamp = get_build_timestamp();
  std::string build_timestamp_text = ctime(&build_timestamp);
  build_timestamp_text.pop_back();  // Remove \n
  m_sink.print(",build-time=\"%s\"", build_timestamp_text.c_str());
  for (const auto& attrib : m_extra_headers) {
    m_sink.print(",%s=\"%s\"", attrib.first.c_str(), attrib.second.c_str());
  }
  m_sink.print("\"\n");
  INFO(DP, 1, Integer, Total Depth)
  INFO(NS, 1, Integer, Number of Samples)
  INFO(SVTYPE, 1, String, Structural Variant Type)
  INFO(MATEID, 1, String, ID of mate breakends)
  INFO(AID, ., Integer, Assembly IDs used in constructing this variant)
  INFO(AMBCOUNT, 1, Integer, Count of alternate locations for this end of an ambiguous breakend)
  INFO(AMBMATES, 1, Integer,
        Count of possible mate locations of an ambiguous breakend)
  INFO(ENTROPYALT, A, Float,
        Shannon entropy of alt allele if longer than 100 bp)
  INFO(SVLEN, 1, Integer,
        Difference in length between REF and ALT alleles)
  INFO(END, 1, Integer,
        End position of the variant described in this record)
  INFO(IMPRECISE, 0, Flag, Imprecise structural variation)
  INFO(CIPOS, 2, Integer, Confidence interval around POS for imprecise variants)
  INFO(CIEND, 2, Integer, Confidence interval around END for imprecise variants)
  INFO(TRANSPOSE, 1, String,
        Transposon FASTA sequence ID that this breakpoint anchor matches)
  INFO(SAS, 1, Float,
        Simple alignment score. Likelihood a breakend is not structural but
            rather aligns to reference simply)
  INFO(FW, A, Float, Percent of forward reads)
  INFO(BQ, A, Integer, Average base quality at this position)
  FORMAT(GT, 1, String, Genotype)
  FORMAT(DP, 1, Integer, Sample Depth)
  FORMAT(AD, ., Integer, Allelic depths for the ref and alt alleles in the order listed)
  FORMAT(ED, 1, Integer, Edit distance)
  FORMAT(OV, A, Integer, Minimum read overlap in assembly)
  FILTER(homologous_breakends,
          The edit_distance between sides of breakpoints was below the
              minimum allowed threshold)
  FILTER(too_many_alleles, The set of possible alleles was too large /
                                supported to be called)
  FILTER(dust_mask, At least 45 bp were considered masked out by DUST)
  FILTER(missing_assembly, No hits were reported by BLAST for this assembly)
  FILTER(non_structural_alignment, BLAST alignment indicates probable SNP)
  FILTER(missing_anchor, One or more anchors were not found by BLAST)
  FILTER(no_unique_anchor, BLAST could not uniquely identify at least one of the anchors for this assembly)
  FILTER(ambiguous_anchor, A BLAST query for this anchor reported multiple ambiguous hits)
  ALT(INS, Insertion)
  ALT(DEL, Deletion)
  for (size_t i = 0; i < m_reference_assembly.scaffold_order.size();
        i++) {
    const scaffold& sc = m_reference_assembly.get_scaffold(
        m_reference_assembly.scaffold_order[i]);
    m_sink.print("##contig=<ID=%s,length=%ld>\n", sc.name.c_str(), sc.len);
  }
  std::string head = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + m_sample_name + "\n";
  m_sink.print("%s", head.c_str());
  //_____ VCF End of Header
  //________________________________________________
}

static std::string check_dot(const dna_sequence& seq,
                             const std::string& prefix) {
  std::string out = prefix + seq.as_string();
  if (out.size() == 0) {
    out = ".";
  }
  return out;
}

void vcf_exporter::nonstruct(const sv_call& call) {
  // Compute depth data
  std::vector<double> depths(call.alleles.size());
  std::vector<int> overlaps(call.alleles.size());
  std::vector<double> fwd_prec(call.alleles.size());
  std::vector<double> avg_qual(call.alleles.size());
  double tot_depth = 0.0;
  bool empty_allele = false;
  for (size_t i = 0; i < call.alleles.size(); i++) {
    double depth_tot = 0.0;
    double fwd_tot = 0.0;
    double qual_tot = 0.0;
    if (call.alleles[i].seq.size() == 0) {
      empty_allele = true;
    }
    for (size_t j = 0; j < call.alleles[i].depth.size(); j++) {
      depth_tot += call.alleles[i].depth[j];
      if (i != 0) {
        fwd_tot += call.alleles[i].fwd[j];
        qual_tot += call.alleles[i].tot_qual[j];
      }
    }
    int overlap = 0;
    for (size_t j = 0; j < call.alleles[i].sub_ids.size(); j++) {
      for (size_t k = 0; k < call.sources.size(); k++) {
        if (call.sources[k].sub_id != call.alleles[i].sub_ids[j]) {
          continue;
        }
        overlap = std::max(overlap, (int)call.sources[k].min_overlap);
      }
    }
    depths[i] = depth_tot / call.alleles[i].depth.size();
    ;
    overlaps[i] = overlap;
    if (i != 0) {
      fwd_prec[i] = fwd_tot / depth_tot;
      avg_qual[i] = qual_tot / depth_tot;
    } else {
      fwd_prec[i] = 0;
      avg_qual[i] = 0;
    }
    tot_depth += depths[i];
  }

  // Set up some helper variables
  std::string sc =
      m_reference_assembly.scaffold_order[call.position.scaffold_id];
  unsigned long pos = call.position.position + 1;
  std::string prefix = "";

  // Deal with the annoying requirement that ref seq can't be empty
  if (call.alleles[0].seq.size() == 0 ||
      (empty_allele && call.alleles.size() > 1)) {
    size_t ref_loc = m_reference.flatten(
        seq_position(call.position.scaffold_id, call.position.position - 1));
    // Pad at front
    pos--;
    prefix += (char)*m_reference.get_dna(ref_loc);
  }

  std::string filter = "PASS";
  std::string gt = "./.";
  bool a0 = (depths[0] / tot_depth > .1);
  bool a1 = (depths.size() > 1 ? depths[1] / tot_depth > .1 : false);
  bool a2 = (depths.size() > 2 ? depths[2] / tot_depth > .1 : false);
  if (call.alleles.size() > 3) {
    filter = "too_many_alleles";
  } else if (a0 && a1 && a2) {
    filter = "too_many_alleles";
  } else if (a0 && a1 && !a2) {
    gt = "0/1";
  } else if (a0 && !a1 && a2) {
    gt = "0/2";
  } else if (a0 && !a1 && !a2) {
    gt = "0/0";
  } else if (!a0 && a1 && a2) {
    gt = "1/2";
  } else if (!a0 && a1 && !a2) {
    gt = "1/1";
  } else if (!a0 && !a1 && a2) {
    gt = "2/2";
  }

  //__________________________________________________________________
  //				VCF Mandatory fixed fields
  // #CHROM
  m_sink.print("%s\t", sc.c_str());
  // POS
  m_sink.print("%ld\t", pos);
  // ID
  m_sink.print(".\t");
  // REF
  m_sink.print("%s\t", check_dot(call.alleles[0].seq, prefix).c_str());
  // ALT
  for (size_t i = 1; i < call.alleles.size(); i++) {
    m_sink.print("%s", check_dot(call.alleles[i].seq, prefix).c_str());
    m_sink.print("%c", i + 1 == call.alleles.size() ? '\t' : ',');
  }
  // QUAL
  m_sink.print("%d\t", 100);
  // FILTER
  m_sink.print("%s\t", filter.c_str());

  //__________________________________________________________________
  //				INFO
  // Number of samples
  m_sink.print("NS=1");
  // Coverage Depth
  m_sink.print(";DP=%d", (int)tot_depth);
  // Assembly ID's
  m_sink.print(";AID=");
  for (size_t i = 0; i < call.sources.size(); i++) {
    m_sink.print("%d", call.sources[i].var_id);
    if (i + 1 != call.sources.size()) {
      m_sink.print(",");
    }
  }
  // Fwd %
  m_sink.print(";FW=");
  for (size_t i = 1; i < fwd_prec.size(); i++) {
    m_sink.print("%2.2f", fwd_prec[i]);
    if (i + 1 != fwd_prec.size()) {
      m_sink.print(",");
    }
  }
  m_sink.print(";BQ=");
  // BQ: Avg Base Qual
  for (size_t i = 1; i < avg_qual.size(); i++) {
    m_sink.print("%d", (int)avg_qual[i]);
    if (i + 1 != avg_qual.size()) {
      m_sink.print(",");
    }
  }
  m_sink.print("\t");

  //__________________________________________________________________
  //				FORMAT
  m_sink.print("GT:DP:AD:OV\t");

  //__________________________________________________________________
  //				SAMPLE
  // GT: Genotype
  m_sink.print("%s:", gt.c_str());
  // DP: Total depth
  m_sink.print("%d:", (int)tot_depth);
  // AD: Allelic depth
  for (size_t i = 0; i < depths.size(); i++) {
    m_sink.print("%d", (int)depths[i]);
    if (i + 1 != depths.size()) {
      m_sink.print(",");
    }
  }
  m_sink.print(":");
  // OV: Overlap
  for (size_t i = 1; i < overlaps.size(); i++) {
    m_sink.print("%d", (int)overlaps[i]);
    if (i + 1 != overlaps.size()) {
      m_sink.print(",");
    }
  }

  //__________________________________________________________________
  //                              NEWLINE
  m_sink.print("\n");
}

void vcf_exporter::struct_event(const struct_var& sv, double ref_depth,
                                bool imprecise) {
  // Only process canonical version
  if (sv.ref_end <= sv.ref_start) return;
  //English added the <= instead of < because some variant have same start/end - wrong!
  // Only output ++ and -- variants
  if (sv.rev_start != sv.rev_end) return;

  //Where the SV starts, 1-based, with anchor
  size_t my_pos = sv.ref_start.position + 1;
  //The anchor base's position
  size_t ref_loc = m_reference.flatten(sv.ref_start);
  seq_position n(sv.ref_end.scaffold_id, sv.ref_end.position);
  size_t end_loc = m_reference.flatten(sv.ref_end);

  //This is the anchor base - still 0 based
  char ref_base = (char)*m_reference.get_dna(ref_loc);
  //Include the anchor base with ref_seq
  dna_sequence ref_seq = dna_sequence(m_reference.get_dna(ref_loc), m_reference.get_dna(end_loc));
  //Doesn't have the anchor base
  dna_sequence alt_seq = sv.assembled.subseq(sv.var_start, sv.var_end - sv.var_start);
  //No idea.. won't use it.
  std::string prefix = "";
  if (!imprecise && (ref_seq.size() == 0 || alt_seq.size() == 0)) {
    // Pad at front
    //offset = 0;ew... the plus one...
    prefix += (char)*m_reference.get_dna(ref_loc + 1);
  }

  std::string svtype;

  //why? -1 here?
  size_t ref_len = end_loc - ref_loc - 1;
  size_t ins_len = sv.var_end - sv.var_start;
  int diff = int(ins_len) - int(ref_len);
  svtype = (diff < 0) ? "DEL" : "INS";

  int tot_depth = (int)sv.depth + (int)ref_depth;
  // need left-alignment...this is bigger than just rolling it over..
  // need to loop - I can't get the comparision correct
  // wtf are these dna_const_iterators
  //if (*ref_base == alt_seq.subseq(alt_seq.size()-1, 1)) {
    //ref_base -= 1;
    //alt_seq = dna_sequence((char)*ref_base) + alt_seq.subseq(0, alt_seq.size() - 1);
    //my_pos -= 1;
  //}
  //Also - whatabout del..?
  // that'd be rollleft while ref_base == end_loc+1
  std::string filter = "PASS";
  std::string gt = "./.";
  // Dumb GT - allele coverage is less than 90% of total coverage is het else hom
  if (tot_depth > 0) {
    bool a0 = (sv.depth / (tot_depth) < .90);
    gt = a0 ? "0/1" : "1/1";
  }

  //__________________________________________________________________
  //				VCF Mandatory fixed fields
  // #CHROM
  std::string sc =
      m_reference_assembly.scaffold_order[sv.ref_start.scaffold_id];
  m_sink.print("%s\t", sc.c_str());
  // POS - which is 0-based including anchor base so no reason to correct
  //     Again, except for point insertions, apparently
  //    OLD - m_sink.print("%ld\t", sv.ref_start.position + offset);
  m_sink.print("%ld\t", my_pos);  // +1 to anchor base
  // ID
  m_sink.print("sv_%d\t", sv.var_id);
  // REF
  if (imprecise) {
    m_sink.print("%c\t", ref_base);
    //m_sink.print("%c\t", (char)*anch_base);
  } else {
    m_sink.print("%s\t", ref_seq.as_string().c_str());
    //m_sink.print("%s%s\t", prefix.c_str(), ref_seq.as_string().c_str());
    //m_sink.print("%c%s\t", (char)*anch_base, ref_seq.as_string().c_str());
  }
  // ALT
  if (imprecise) {
    m_sink.print("<%s>\t", svtype.c_str());
  } else {
    m_sink.print("%c%s\t", ref_base, alt_seq.as_string().c_str());
    //m_sink.print("%c%s\t", (char)*anch_base, alt_seq.as_string().c_str());
  }
  // QUAL
  m_sink.print("%d\t", 100);
  // FILTER
  m_sink.print("%s\t", sv.filter != "" ? sv.filter.c_str() : filter.c_str());

  //__________________________________________________________________
  //				INFO
  // Coverage Depth
  m_sink.print("NS=1");
  m_sink.print(";DP=%d", tot_depth);
  m_sink.print(";SVTYPE=%s", svtype.c_str());
  m_sink.print(";END=%ld", sv.ref_end.position);
  m_sink.print(";SVLEN=%d", diff);
  m_sink.print(";AID=%d", sv.var_id);
  if (imprecise) {
    m_sink.print(";IMPRECISE;CIPOS=0,0;CIEND=0,0");
  }
  m_sink.print("\t");

  //__________________________________________________________________
  //				FORMAT
  m_sink.print("GT:DP:AD:ED:OV\t");

  //__________________________________________________________________
  //				SAMPLE
  // GT: Genotype
  m_sink.print(gt);
  // DP: Total depth
  m_sink.print(":%d", tot_depth);
  // AD: Allelic depth
  if (ref_depth < 0) {
    m_sink.print(":.,%d", (int)sv.depth);
  } else {
    m_sink.print(":%d,%d", (int)ref_depth, (int)sv.depth);
  }
  // ED: Edit distance (. if < 0)
  int ref_diff = sv_compute_edit_distance(sv, m_reference);
  if (ref_diff < 0) {
    m_sink.print(":.");
  } else {
    m_sink.print(":%d", ref_diff);
  }
  // OV: Min overlap
  m_sink.print(":%d", sv.min_overlap);

  //__________________________________________________________________
  //                              NEWLINE
  m_sink.print("\n");
}

void vcf_exporter::breakend(const struct_var& sv, double ref_depth) {
  size_t ref_loc = m_reference.flatten(sv.ref_start);
  char ref_base = (char)*m_reference.get_dna(ref_loc);

  bool is_left = (sv.ref_start < sv.ref_end);
  std::string filter;
  if (sv.filter.empty()) {
    filter = "PASS";
  } else {
    filter = sv.filter;
  }
  std::string amb_field;
  if (sv.filter == "left_transposable_element") {
    return;
  }
  if (sv.ambiguous_count) {
    if ((sv.ambiguous_side == struct_var::amb_left && is_left) ||
        (sv.ambiguous_side == struct_var::amb_right && !is_left)) {
      filter = "ambiguous_anchor";
      amb_field = "AMBCOUNT";
    } else {
      amb_field = "AMBMATES";
    }
  }

  int tot_depth = (int)sv.depth + (int)ref_depth;

  //__________________________________________________________________
  //				VCF Mandatory fixed fields
  // #CHROM
  std::string sc =
      m_reference_assembly.scaffold_order[sv.ref_start.scaffold_id];
  m_sink.print("%s\t", sc.c_str());
  // POS
  m_sink.print("%ld\t", sv.ref_start.position + 1);
  // ID
  m_sink.print("bnd_%d\t", sv.var_id * 2 + (is_left ? 0 : 1));
  // REF
  m_sink.print("%c\t", ref_base);
  // ALT
  dna_sequence middle =
      sv.assembled.subseq(sv.var_start, sv.var_end - sv.var_start);
  seq_position position_other = sv.ref_end;
  std::string scaffold_other =
      m_reference_assembly.scaffold_order[position_other.scaffold_id];
  std::string other =
      printstring("%s:%ld", scaffold_other.c_str(), position_other.position + 1);
  if (sv.rev_start) {
    if (sv.rev_end) {
      // --
      m_sink.print("]%s]%s%c\t", other.c_str(),
                   middle.rev_comp().as_string().c_str(), ref_base);
    } else {
      // -+
      m_sink.print("[%s[%s%c\t", other.c_str(),
                   middle.rev_comp().as_string().c_str(), ref_base);
    }
  } else {
    if (sv.rev_end) {
      // +-
      m_sink.print("%c%s]%s]\t", ref_base, middle.as_string().c_str(),
                   other.c_str());
    } else {
      // ++
      m_sink.print("%c%s[%s[\t", ref_base, middle.as_string().c_str(),
                   other.c_str());
    }
  }
  // QUAL
  m_sink.print("%d\t", 100);
  // FILTER
  m_sink.print("%s\t", filter.c_str());

  //__________________________________________________________________
  //				INFO
  // Coverage Depth
  m_sink.print("NS=1");
  m_sink.print(";DP=%d", tot_depth);
  m_sink.print(";SVTYPE=BND");
  m_sink.print(";AID=%d", sv.var_id);
  m_sink.print(";MATEID=bnd_%d", sv.var_id * 2 + (is_left ? 1 : 0));
  // Assume the assembly must have some identity, so 0 means blast did not run,
  // e.g. homologous_breakends
  if (sv.simple_alignment_score >= 0.01) {
    m_sink.print(";SAS=%4.2f", sv.simple_alignment_score);
  }
  if (!amb_field.empty()) {
    m_sink.print(";%s=%zu", amb_field.c_str(), sv.ambiguous_count);
  }
  if (!sv.transpose.empty()) {
    m_sink.print(";TRANSPOSE=%s", sv.transpose.c_str());
  }
  if (middle.size() > 100) {
    m_sink.print(";ENTROPYALT=%.3f", compute_entropy(middle.as_string()));
  }
  m_sink.print("\t");

  //__________________________________________________________________
  //				FORMAT
  m_sink.print("GT:DP:AD:ED:OV\t");

  //__________________________________________________________________
  //				SAMPLE
  // GT: Genotype
  m_sink.print("./.");
  // DP: Total depth
  m_sink.print(":%d", (int)sv.depth);
  // AD: Allelic depth
  if (ref_depth < 0) {
    m_sink.print(":.,%d", (int)sv.depth);
  } else {
    m_sink.print(":%d,%d", (int)ref_depth, (int)sv.depth);
  }
  // ED: Edit distance (. if < 0)
  int ref_diff = sv_compute_edit_distance(sv, m_reference);
  if (ref_diff < 0) {
    m_sink.print(":.");
  } else {
    m_sink.print(":%d", ref_diff);
  }
  // OV: Min overlap
  m_sink.print(":%d", sv.min_overlap);

  //__________________________________________________________________
  //                              NEWLINE
  m_sink.print("\n");
}

void vcf_exporter::write(const std::string& key, const std::string& value) {
  // SPLOG("Getting a value");
  //_____ VCF Body ______________________________
  sv_call svc;
  msgpack_deserialize(svc, value);

  if (svc.alleles.size()) {
    nonstruct(svc);
  } else {
    struct_var sv = svc.sources[0];
    if (m_use_events && svc.sv_ref_depth != -1 && sv.ambiguous_count == 0 &&
        sv.rev_start == sv.rev_end) {
      struct_event(sv, svc.sv_ref_depth, false);
    } else {
      breakend(sv, svc.sv_ref_depth);
    }
  }
  //_____ VCF End of Body ______________________________
}

void vcf_exporter::write_footer() {
  SPLOG("vcf_exporter::write_footer> VCF export complete.");
}

double vcf_exporter::compute_entropy(const std::string& a_string) {
  std::map<char, unsigned> char_freq_map;
  typedef std::map<char, unsigned>::value_type char_freq_value_t;
  std::for_each(a_string.cbegin(), a_string.cend(),
                [&char_freq_map](char c) { char_freq_map[c]++; });

  double entropy = std::accumulate(
      char_freq_map.cbegin(), char_freq_map.cend(), 0.0,
      [&a_string](double cumulative_entropy, char_freq_value_t freq_map_entry) {
        double normalized_freq =
            static_cast<double>(freq_map_entry.second) / a_string.size();
        return cumulative_entropy +
               normalized_freq * log(normalized_freq) / log(2.0);
      });

  SPLOG_P(LOG_DEBUG, "vcf_exporter::compute_entropy> Computed entropy of %f",
          -entropy);
  return -entropy;
}
