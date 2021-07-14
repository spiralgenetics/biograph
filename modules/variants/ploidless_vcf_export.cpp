#include "modules/variants/ploidless_vcf_export.h"
#include "modules/variants/scaffold.h"
#include "tools/build_stamp.h"

namespace variants {

namespace {

const char* remove_beginning_newline(const char* str) {
  CHECK_EQ(str[0], '\n');
  return str + 1;
}

}

void ploidless_vcf_export::on_assembly(assembly_ptr a) {
  if (a->matches_reference) {
    return;
  }

  int min_depth = 0;
  int min_pair_depth = 0;
  if (!a->coverage.empty()) {
    min_depth = container_min(a->coverage);
  }
  if (!a->pair_coverage.empty()) {
    min_pair_depth = container_min(a->pair_coverage);
  }

  std::string var_seq = a->seq.as_string();
  std::string ref_seq =
      m_options.scaffold->subscaffold_str(a->left_offset, a->right_offset - a->left_offset);

  aoffset_t left_offset = a->left_offset;

  if (var_seq.empty() || ref_seq.empty()) {
    if (left_offset == 0) {
      std::string base_to_add = m_options.scaffold->subscaffold_str(a->right_offset, 1);
      var_seq += base_to_add;
      ref_seq += base_to_add;
    } else {
      std::string base_to_add = m_options.scaffold->subscaffold_str(a->left_offset - 1, 1);
      var_seq = base_to_add + var_seq;
      ref_seq = base_to_add + ref_seq;
      --left_offset;
    }
  }

  std::string unphased_genotype;
  std::string phased_genotype;
  std::string filter_str;
  switch (a->strand_count) {
    case 0:
      filter_str = "FILTERED";
      unphased_genotype = "0/1";
      phased_genotype = "0|1";
      break;
    case 1:
      filter_str = "PASS";
      unphased_genotype = "0/1";
      phased_genotype = "0|1";
      break;
    case 2:
      filter_str = "PASS";
      unphased_genotype = "1/1";
      phased_genotype = "1|1";
      break;
    default:
      CHECK_LE(a->strand_count, 2) << "Outputting VCF with more than 2 strands not supported";
  }
  mem_io line("", track_alloc("ploidless_vcf_export:line"));

  //__________________________________________________________________
  //				VCF Mandatory fixed fields
  // #CHROM
  line.print("%s\t", m_scaffold_name.c_str());
  // POS
  line.print("%d\t", left_offset + 1 /* VCFs are 1-indexed */);
  // ID
  line.print(".\t");
  // REF
  line.print("%s\t", ref_seq.c_str());
  // ALT
  line.print("%s\t", var_seq.c_str());
  // QUAL
  line.print("%d\t", 100);
  // FILTER
  line.print("PASS\t");

  //__________________________________________________________________
  //				INFO

  // Always output NS to avoid a potentially empty INFO column.
  line.print("NS=1");

  // Optional AIDs
  if (m_options.output_assembly_ids) {
    line.print(";AID=");
    line.print("%ld", a->assembly_id);
    for (size_t i = 0; i < a->merged_assembly_ids.size(); i++) {
      line.print(",%ld", a->merged_assembly_ids[i]);
    }
  }
  // SVLEN and SVTYPE
  if (var_seq.size() >= m_options.vcf_sv_size_threshold ||
      ref_seq.size() >= m_options.vcf_sv_size_threshold) {
    aoffset_t right_offset = left_offset + ref_seq.size();
    aoffset_t rightmost_ref_base = right_offset - 1;
    line.print(";END=%d", rightmost_ref_base + 1 /* VCFs are 1-indexed */);
    int svlen = int(var_seq.size()) - int(ref_seq.size());
    line.print(";SVLEN=%d", svlen);
    if (svlen < 0) {
      line.print(";SVTYPE=DEL");
    } else if (svlen > 0) {
      line.print(";SVTYPE=INS");
    } else {
      line.print(";SVTYPE=CPX");
    }
  }
  if (m_options.use_bidir_tracer) {
    line.print(";GENBY=%s", a->tags.to_string_short().c_str());
  } else {
    if (a->tags.contains("POP")) {
      CHECK(m_options.pop_trace_anchor_drop || m_options.use_pop_tracer);
      line.print(";POP");
    }
  }
  line.print("\t");

  //__________________________________________________________________
  //				FORMAT
  //
  // NOTE: Per the VCF spec, GT *must* be the first field if it is present. Because reasons.
  // https://samtools.github.io/hts-specs/VCFv4.2.pdf
  //
  line.print("GT:PG:GQ:PI:OV:DP:AD:PDP:PAD");
  if (m_options.output_ml_features) {
    line.print(":LASCORE:LAREFSPAN:LARANCH:LALANCH:LAREFGC:LAALTGC:LAALTSEQLEN:NUMASM");
  }
  line.print("\t");

  //__________________________________________________________________
  //				SAMPLE
  // GT: Genotype
  line.print("%s:", unphased_genotype.c_str());
  // PG: Phased Genotype
  line.print("%s:", phased_genotype.c_str());
  // GQ: Genotype quality
  line.print("%d:", int(a->genotype_quality * 100));
  // PI: Phase ID
  line.print("%ld:", a->assembly_id);
  // OV: Overlap
  line.print("%d:", a->min_overlap);
  // DP: total depth across all alleles
  line.print("%d:", a->other_depth + min_depth + a->ref_depth);
  // AD: allelic depth
  line.print("%d,%d:", a->ref_depth, min_depth);
  // PDP: total pair-confirmed-read depth across all alleles
  // TODO(nils): Add ref pair coverage.
  line.print("%d:", a->other_pair_depth + min_pair_depth);
  // PAD: pair-confirmed-read allelic depth for this and reference
  // TODO(nils): Add ref pair coverage.
  line.print("%d,%d", 0, min_pair_depth);

  if (m_options.output_ml_features) {
    CHECK(a->ml_features) << "Missing ML features on assembly? " << *a;
    const assembly_ml_features& f = *a->ml_features;
    line.print(":%d", f.score);
    line.print(":%d", f.refspan);
    line.print(":%d", f.lanch);
    line.print(":%d", f.ranch);
    line.print(":%f", f.refgc);
    line.print(":%f", f.altgc);
    line.print(":%lu", f.alt_seq.size());
    line.print(":%lu", a->merged_assembly_ids.size() + 1);
  }
  //__________________________________________________________________
  //                              NEWLINE
  line.print("\n");

  m_output(line.str());
}

std::string ploidless_vcf_export::header(const assemble_options& options,
                                         const std::map<std::string, std::string>& extra_headers,
                                         const std::string sample_name) {
  // VCF version
  std::string header = "##fileformat=VCFv4.1\n";

  // Today's date as 20180726
  static char date[10];
  time_t rawtime;
  struct tm* timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(date, 10, "%Y%m%d", timeinfo);
  header += printstring("##fileDate=%s\n", date);

  // BioGraph version and runtime options
  header += R"VER(##source="Spiral Genetics BioGraph")VER";
  header += printstring(",version=\"%s\"", biograph_current_version.make_string().c_str());
  header += printstring(",description=\"build-revision='%s%s'", get_build_scm_revision().c_str(),
                        build_is_clean() ? "" : " (unclean workspace)");
  time_t build_timestamp = get_build_timestamp();
  std::string build_timestamp_text = ctime(&build_timestamp);
  build_timestamp_text.pop_back();  // Remove \n
  header += printstring(",build-time='%s'", build_timestamp_text.c_str());
  for (const auto& attrib : extra_headers) {
    header += printstring(",%s='%s'", attrib.first.c_str(), attrib.second.c_str());
  }
  header += "\"\n";

  // Reference
  header += printstring("##reference=%s\n", options.ref->path().c_str());

  // Optional Assembly IDs
  if (options.output_assembly_ids) {
    header += remove_beginning_newline(R"AID(
##INFO=<ID=AID,Number=.,Type=Integer,Description="Assembly IDs used in constructing this variant">
)AID");
  }

  if (options.use_bidir_tracer) {
    header += remove_beginning_newline(R"POP(
##INFO=<ID=GENBY,Number=1,Type=String,Description="Type of tracer used that discovered this variant">
)POP");

  } else if (options.use_pop_tracer || options.pop_trace_anchor_drop) {
    header += remove_beginning_newline(R"POP(
##INFO=<ID=POP,Number=0,Type=Flag,Description="Found using pop tracer">
)POP");
  }

  // Other static fields
  header += remove_beginning_newline(R"HEADER(
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural Variant Type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample Depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=PDP,Number=1,Type=Integer,Description="Sample Pair Depth">
##FORMAT=<ID=PAD,Number=.,Type=Integer,Description="Allelic pair depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=OV,Number=1,Type=Integer,Description="Minimum read overlap in assembly">
##FORMAT=<ID=PG,Number=1,Type=String,Description="Phased genotype">
##FORMAT=<ID=PI,Number=1,Type=Integer,Description="Phase group">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
)HEADER");

  if (options.output_ml_features) {
    header += remove_beginning_newline(R"HEADER(
##FORMAT=<ID=LASCORE,Number=1,Type=Integer,Description="Score of longest assembly">
##FORMAT=<ID=LAREFSPAN,Number=1,Type=Integer,Description="Ref span length of longest assembly">
##FORMAT=<ID=LARANCH,Number=1,Type=Integer,Description="Right anchor length of longest assembly">
##FORMAT=<ID=LALANCH,Number=1,Type=Integer,Description="Leftanchor length of longest assembly">
##FORMAT=<ID=LAREFGC,Number=1,Type=Float,Description="Portion of G/C bases in ref span of longest assembly">
##FORMAT=<ID=LAALTGC,Number=1,Type=Float,Description="Portion of G/C bases in sequence of longest assembly">
##FORMAT=<ID=LAALTSEQLEN,Number=1,Type=Integer,Description="Sequence length of longest assembly">
##FORMAT=<ID=NUMASM,Number=1,Type=Integer,Description="Number of assemblies that independently produced this variant">
)HEADER");
  }

header += remove_beginning_newline(R"HEADER(
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DEL,Description="Deletion">
)HEADER");

  // Chromosomes
  const auto& refasm = options.ref->get_assembly();
  for (size_t i = 0; i < refasm.scaffold_order.size(); i++) {
    const auto& sc = refasm.get_scaffold(refasm.scaffold_order[i]);
    header += printstring("##contig=<ID=%s,length=%ld>\n", sc.name.c_str(), sc.len);
  }

  // Header + Sample ID
  header += remove_beginning_newline(R"HEADER(
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	)HEADER");
  header += printstring("%s\n", sample_name.c_str());
  return header;
}

}  // namespace variants
