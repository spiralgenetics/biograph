#include "modules/bio_base/seqset_export.h"
#include "modules/bio_base/seqset.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_format/fastq.h"
#include "modules/io/file_io.h"
#include "modules/io/log.h"
#include "modules/io/make_unique.h"
#include "modules/io/zip.h"
#include "modules/main/main.h"

#include <stdexcept>

// helper to only update progress when the delta is > 0.01%
static void update_progress(const float& new_progress) {
  static float prev_progress = 0;
  if (fabs(new_progress - prev_progress) > 0.0001) {
    prev_progress = new_progress;
    print_progress(new_progress);
  }
}

class SEQSETExportFastqMain : public Main {
 public:
  SEQSETExportFastqMain() {
    m_usage =
        "%1% version %2%\n\n"
        "Usage: %1% [OPTIONS] --in <file.seqset> --readmap <file.readmap> "
        "--out <file.fastq1> --pair <file.fastq2>\n\n"
        "Write out all the reads from a seqset + redmap.\n";
  }

 protected:
  const product_version& get_version() override { return biograph_current_version; }
  void add_args() override;
  int run(po::variables_map vars) override;

 private:
  std::string m_seqset_file;
  std::string m_readmap_file;
  std::string m_fastq_out1;
  std::string m_fastq_out2;

  void do_export();
};

void SEQSETExportFastqMain::add_args() {
  m_options.add_options()("in", po::value(&m_seqset_file)->required(),
                          "Seqset file")(
      "readmap", po::value(&m_readmap_file)->required(),
      "Readmap to get reads from")("out", po::value(&m_fastq_out1)->required(),
                                   "Destination fastq for first part of pairs")(
      "pair", po::value(&m_fastq_out2)->required(),
      "Destination fastq for second part of pairs and unpaired entries");
}

int SEQSETExportFastqMain::run(po::variables_map vars) {
  initialize_app("");
  do_export();

  return 0;
}

namespace {

class export_chunk {
 public:
  export_chunk(uint32_t start, uint32_t limit)
      : m_start(start), m_limit(limit), m_outbuf("", track_alloc("export_fastq:chunk")) {
  }

  ~export_chunk() {
    CHECK(!m_exporter);
    CHECK_EQ(0, m_outbuf.size());
  }

  void flush() {
    if (!m_exporter) {
      CHECK_EQ(0, m_outbuf.size());
      return;
    }
    m_exporter.reset();
    m_zipper->close();
    m_zipper.reset();
  }

  void write_to(file_writer& out) {
    CHECK(!m_exporter);
    CHECK(!m_zipper);
    out.write(m_outbuf.buffer(), m_outbuf.size());
    m_outbuf.clear();
  }

  fastq_exporter& get_exporter() {
    if (!m_exporter) {
      m_zipper.emplace(m_outbuf);
      m_exporter.emplace(*m_zipper);
    }
    return *m_exporter;
  }

  uint32_t start_read_id() const { return m_start; }
  uint32_t limit_read_id() const { return m_limit; }

 private:
  const uint32_t m_start, m_limit;
  mem_io m_outbuf;
  boost::optional<zip_writer> m_zipper;
  boost::optional<fastq_exporter> m_exporter;
};

struct export_state {
  export_state(file_writer& o1, file_writer& o2) : out1(o1), out2(o2) {}

  std::mutex mu;

  file_writer& out1;
  file_writer& out2;

  uint32_t next_read_id_1 = 0;
  uint32_t next_read_id_2 = 0;

  std::map<uint32_t /* start read id */, std::unique_ptr<export_chunk>> chunks1;
  std::map<uint32_t /* start read id */, std::unique_ptr<export_chunk>> chunks2;
};

void flush_next_chunks_part(
    std::mutex& mu, uint32_t& next_read_id,
    std::map<uint32_t /* start read id */, std::unique_ptr<export_chunk>>& chunks,
    file_writer& out) {
  std::unique_lock<std::mutex> l(mu);
  while (!chunks.empty()) {
    auto chunk_it = chunks.begin();
    if (chunk_it->first != next_read_id) {
      return;
    }
    std::unique_ptr<export_chunk> chunk = std::move(chunk_it->second);
    CHECK_EQ(chunk_it->first, chunk->start_read_id());
    chunks.erase(chunk_it);

    l.unlock();
    chunk->write_to(out);
    l.lock();

    CHECK_EQ(chunk->start_read_id(), next_read_id);
    next_read_id = chunk->limit_read_id();
  }
}

void flush_next_chunks(export_state& state) {
  flush_next_chunks_part(state.mu, state.next_read_id_1, state.chunks1, state.out1);
  flush_next_chunks_part(state.mu, state.next_read_id_2, state.chunks2, state.out2);
}

class fastq_export_worker : public seqset_export_worker {
 public:
  fastq_export_worker(export_state& state) : m_state(state) {}
  void start_chunk(uint32_t start, uint32_t limit) override {
    CHECK(!m_chunk1);
    CHECK(!m_chunk2);
    m_chunk1 = make_unique<export_chunk>(start, limit);
    m_chunk2 = make_unique<export_chunk>(start, limit);
  }

  void output_paired(uint32_t this_read_id, dna_slice r1, dna_slice r2) override {
    unaligned_reads reads1;
    reads1.emplace_back();
    reads1.back().sequence = r1.as_string();
    reads1.back().quality = std::string(r1.size(), '"');

    unaligned_reads reads2;
    reads2.emplace_back();
    reads2.back().sequence = r2.as_string();
    reads2.back().quality = std::string(r2.size(), '"');

    read_id id;
    id.pair_name = std::to_string(this_read_id);

    CHECK(m_chunk1);
    CHECK(m_chunk2);
    m_chunk1->get_exporter().write(id, reads1);
    m_chunk2->get_exporter().write(id, reads2);
  }

  void output_unpaired(uint32_t this_read_id, dna_slice r) override {
    unaligned_reads reads;
    reads.emplace_back();
    reads.back().sequence = r.as_string();
    reads.back().quality = std::string(r.size(), '!');

    read_id id;
    id.pair_name = std::to_string(this_read_id);
    m_chunk2->get_exporter().write(id, reads);
  }

  void done_chunk() override {
    CHECK(m_chunk1);
    m_chunk1->flush();
    CHECK(m_chunk2);
    m_chunk2->flush();

    {
      std::lock_guard<std::mutex> l(m_state.mu);
      bool did_insert;
      did_insert = m_state.chunks1.emplace(m_chunk1->start_read_id(), std::move(m_chunk1)).second;
      CHECK(did_insert);
      did_insert = m_state.chunks2.emplace(m_chunk2->start_read_id(), std::move(m_chunk2)).second;
      CHECK(did_insert);
    }
    m_chunk1.reset();
    m_chunk2.reset();
    flush_next_chunks(m_state);
    CHECK(!m_chunk1);
    CHECK(!m_chunk2);
  }

 private:
  export_state& m_state;
  std::unique_ptr<export_chunk> m_chunk1;
  std::unique_ptr<export_chunk> m_chunk2;
};

}  // namespace

void SEQSETExportFastqMain::do_export() {
  auto ss_f = std::make_shared<seqset>(m_seqset_file);
  readmap m_readmap(ss_f, m_readmap_file);

  std::cout << "\nLoading seqset\n";
  ss_f->membufs().cache_in_memory(update_progress);

  file_writer out1(m_fastq_out1);
  file_writer out2(m_fastq_out2);

  seqset_export e(ss_f.get(), &m_readmap, m_tmp_dir);

  std::cout << "\nPreparing for export\n";
  e.prepare(update_progress);

  std::cout << "\nExporting paired data\n";
  {
    export_state paired_state(out1, out2);
    e.write_paired(
        [&]() -> std::unique_ptr<seqset_export_worker> {
          return make_unique<fastq_export_worker>(paired_state);
        },
        update_progress);

    CHECK_EQ(paired_state.chunks1.size(), 0);
    CHECK_EQ(paired_state.next_read_id_1, m_readmap.size());
    CHECK_EQ(paired_state.chunks2.size(), 0);
    CHECK_EQ(paired_state.next_read_id_2, m_readmap.size());
  }

  std::cout << "\nExporting unpaired data\n";
  {
    export_state unpaired_state(out1, out2);
    e.write_unpaired(
        [&]() -> std::unique_ptr<seqset_export_worker> {
          return make_unique<fastq_export_worker>(unpaired_state);
        },
        update_progress);
    CHECK_EQ(unpaired_state.chunks1.size(), 0);
    CHECK_EQ(unpaired_state.next_read_id_1, m_readmap.size());
    CHECK_EQ(unpaired_state.chunks1.size(), 0);
    CHECK_EQ(unpaired_state.next_read_id_2, m_readmap.size());
  }

  std::cout << "\nExport complete\n";
}

std::unique_ptr<Main> export_fastq_main() {
  return std::unique_ptr<Main>(new SEQSETExportFastqMain);
}
