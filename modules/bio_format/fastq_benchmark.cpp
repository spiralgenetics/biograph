#include "benchmark/benchmark.h"
#include "modules/io/io.h"
#include "modules/io/file_io.h"
#include "modules/io/zip.h"
#include "modules/bio_format/fastq.h"
#include "base/base.h"

#include <string.h>

namespace {

class fastq_generator : public read_wrapper {
 public:
  fastq_generator(std::string pattern)
      : m_data(pattern) {
    m_data_ptr = m_data.data();
    m_data_end = m_data.data() + m_data.size();
  }

  int base_read(char* buf, size_t len) override;

  size_t get_tot_processed() const { return m_tot_processed; }

 private:
  std::string m_data;
  size_t m_tot_processed = 0;

  const char* m_data_ptr = nullptr;
  const char* m_data_end = nullptr;
};

int fastq_generator::base_read(char* buf, size_t len)
{
  len = std::min<size_t>(len, m_data_end - m_data_ptr);
  DCHECK_GT(len, 0);
  memcpy(buf, m_data_ptr, len);
  m_data_ptr += len;
  
  if (m_data_ptr == m_data_end) {
    m_data_ptr = m_data.data();
  }
  m_tot_processed += len;
  return len;
}

void read_single(fastq_reader& fq) {
  read_id id;
  unaligned_reads value;
  fq.read(id, value);
}

void read_paired(fastq_reader& fq1, fastq_reader& fq2) {
  read_id id1, id2;
  unaligned_reads value;
  value.emplace_back();
  fq1.read(id1, value.back());
  value.emplace_back();
  fq2.read(id2, value.back());
}

}  // namespace

static void BM_read_fastq(benchmark::State& state) {
  fastq_generator gen(
      "@6000:1:1101:1049:2117/1\n"
      "GAAACCGTTGCAGGAAACGTAACCGCGGCAGCGTCAGACACAGCCAGTTGTGTCGATTGCGGTTCCACAGGC"
      "GCTTCCACTGTGCGGCTTTTTATATATA\n"
      "+\n"
      "@<@D:==DHHF>FHIG92A<+C@DEAFHAHABG;C//=ACEE?6;;>.;>;=(-9,5@?CB@272443:<??"
      ";@&55@AC@C##################\n");

  fastq_reader fq(gen);
  
  while (state.KeepRunning()) {
    read_single(fq);
  }
  state.SetBytesProcessed(gen.get_tot_processed());
}

BENCHMARK(BM_read_fastq);

static void BM_read_fastq_paired(benchmark::State& state) {
  fastq_generator gen1(
      "@6000:1:1101:1049:2117/1\n"
      "GAAACCGTTGCAGGAAACGTAACCGCGGCAGCGTCAGACACAGCCAGTTGTGTCGATTGCGGTTCCACAGGC"
      "GCTTCCACTGTGCGGCTTTTTATATATA\n"
      "+\n"
      "@<@D:==DHHF>FHIG92A<+C@DEAFHAHABG;C//=ACEE?6;;>.;>;=(-9,5@?CB@272443:<??"
      ";@&55@AC@C##################\n");
  fastq_generator gen2(
      "@6000:1:1101:1049:2117/2\n"
      "GAAACCGTTGCAGGAAACGTAACCGCGGCAGCGTCAGACACAGCCAGTTGTGTCGATTGCGGTTCCACAGGC"
      "GCTTCCACTGTGCGGCTTTTTATATATA\n"
      "+\n"
      "@<@D:==DHHF>FHIG92A<+C@DEAFHAHABG;C//=ACEE?6;;>.;>;=(-9,5@?CB@272443:<??"
      ";@&55@AC@C##################\n");

  fastq_reader fq1(gen1);
  fastq_reader fq2(gen2);
  
  while (state.KeepRunning()) {
    read_paired(fq1, fq2);
  }
  state.SetBytesProcessed(gen1.get_tot_processed() + gen2.get_tot_processed());
}

BENCHMARK(BM_read_fastq_paired);

static void BM_read_fastq_gz(benchmark::State& state) {
  fastq_generator gen(slurp_file(
      "/share/datasets/panels/SRR081224/100/SRR081224_100Genes_r1.fastq.gz"));

  zip_reader zr(gen);
  fastq_reader fq(zr);
  
  while (state.KeepRunning()) {
    read_single(fq);
  }
  state.SetBytesProcessed(gen.get_tot_processed());
}

BENCHMARK(BM_read_fastq_gz);

static void BM_read_fastq_gz_paired(benchmark::State& state) {
  fastq_generator gen1(slurp_file(
      "/share/datasets/panels/SRR081224/100/SRR081224_100Genes_r1.fastq.gz"));
  fastq_generator gen2(slurp_file(
      "/share/datasets/panels/SRR081224/100/SRR081224_100Genes_r2.fastq.gz"));

  zip_reader zr1(gen1);
  zip_reader zr2(gen2);
  fastq_reader fq1(zr1);
  fastq_reader fq2(zr2);
  
  while (state.KeepRunning()) {
    read_paired(fq1, fq2);
  }
  state.SetBytesProcessed(gen1.get_tot_processed() + gen2.get_tot_processed());
}

BENCHMARK(BM_read_fastq_gz_paired);

BENCHMARK_MAIN();
