#include "modules/bio_base/seqset_export.h"
#include "modules/io/config.h"
#include "modules/io/parallel.h"

seqset_export_worker::~seqset_export_worker() = default;

void seqset_export::prepare(progress_handler_t progress) {
  CHECK(!m_flat);
  CHECK(!m_tmp_dir.empty());

  static std::atomic<size_t> g_counter{0};
  boost::filesystem::create_directories(m_tmp_dir);
  std::string flat_path = m_tmp_dir + "/export_flat" + std::to_string(g_counter.fetch_add(1));
  {
    spiral_file_create_mmap c(flat_path);
    seqset_flat_builder b(&m_seqset->get_seqset());
    b.build(c.create(), progress);
  }

  spiral_file_open_mmap o(flat_path);
  m_flat.emplace(o.open(), &m_seqset->get_seqset());
}

namespace {

class seqset_export_local : public parallel_local {
 public:
  seqset_export_local(const readmap* rm, const seqset_flat* flat, bool paired,
                      const std::function<std::unique_ptr<seqset_export_worker>()>& worker_f)
      : m_readmap(rm), m_flat(flat), m_paired(paired), m_worker(worker_f()) {}

  void process_range(size_t start, size_t limit) {
    m_worker->start_chunk(start, limit);
    for (size_t read_id = start; read_id != limit; ++read_id) {
      if (!m_readmap->get_is_forward(read_id)) {
        continue;
      }

      uint32_t mate_read_id = std::numeric_limits<uint32_t>::max();
      if (m_paired) {
        if (!m_readmap->has_mate(read_id)) {
          continue;
        }

        mate_read_id = m_readmap->get_mate(read_id);
        if (read_id > mate_read_id) {
          // Only output each pair once.
          continue;
        }
      } else {
        if (m_readmap->has_mate(read_id)) {
          continue;
        }
      }

      uint64_t read_entry = m_readmap->index_to_entry(read_id);
      int read_len = m_readmap->get_readlength(read_id);
      dna_slice read_sequence = m_flat->get(read_entry);

      if (m_paired) {
        uint64_t mate_entry = m_readmap->index_to_entry(mate_read_id);
        int mate_len = m_readmap->get_readlength(mate_read_id);
        dna_slice mate_sequence = m_flat->get(mate_entry);

        m_worker->output_paired(read_id, read_sequence.subseq(0, read_len),
                                mate_sequence.subseq(0, mate_len));
      } else {
        m_worker->output_unpaired(read_id, read_sequence.subseq(0, read_len));
      }
    }
    m_worker->done_chunk();
  }

  void flush() override {
    m_worker.reset();
  };

 private:
  const readmap* m_readmap;
  const seqset_flat* m_flat;
  const bool m_paired;
  std::unique_ptr<seqset_export_worker> m_worker;
};

}  // namespace

void seqset_export::write_paired(
    const std::function<std::unique_ptr<seqset_export_worker>()>& worker_f,
    progress_handler_t progress) {
  CHECK(m_flat);

  size_t num_reads = m_readmap->size();
  parallel_for(  //
      0, num_reads,
      [&](size_t start, size_t end, parallel_state& ps) {
        seqset_export_local* local =
            ps.get_local<seqset_export_local>(m_readmap, &*m_flat, true /* paired */, worker_f);
        local->process_range(start, end);
      },
      progress);
}

void seqset_export::write_unpaired(
    const std::function<std::unique_ptr<seqset_export_worker>()>& worker_f,
    progress_handler_t progress) {
  CHECK(m_flat);

  // Now, output unpaired entries.
  size_t num_reads = m_readmap->size();
  parallel_for(  ///
      0, num_reads,
      [&](size_t start, size_t end, parallel_state& ps) {
        seqset_export_local* local = ps.get_local<seqset_export_local>(
            m_readmap, &*m_flat, false /* not paired */, worker_f);
        local->process_range(start, end);
      },
      progress);
}
