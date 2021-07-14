#pragma once

#include "modules/bio_base/reference.h"
#include "modules/bio_base/sv_call.h"
#include "modules/bio_format/exporter.h"

class vcf_exporter : public exporter {
 public:
  vcf_exporter(writable& byte_sink, const std::string& ref_name,
               const std::map<std::string, std::string>& extra_headers,
               bool use_events = false,
               const std::string& sample_name = "SAMPLE");

  vcf_exporter(writable& byte_sink, const std::string& ref_name,
               bool use_events = false,
               const std::string& sample_name = "SAMPLE");

  vcf_exporter(writable& byte_sink, bool /*unused*/,
               const std::string& serialized_ref_name,
               const std::string& sample_name = "SAMPLE");

  void write(const std::string& key, const std::string& value) override;
  void write_header() override;
  void write_footer() override;

  static double compute_entropy(const std::string& a_string);

 private:
  void breakend(const struct_var& sv, double ref_depth);
  void struct_event(const struct_var& sv, double ref_depth, bool imprecise);
  void nonstruct(const sv_call& sv);

  bool m_use_events;
  reference m_reference;
  const reference_assembly& m_reference_assembly;
  std::map<std::string, std::string> m_extra_headers;
  std::string m_sample_name;
};
