#pragma once

#include "modules/bio_base/seqset.h"
#include "modules/bio_base/seqset_flat.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_base/unaligned_read.h"

class seqset_export_worker {
public:
  virtual void start_chunk(uint32_t start, uint32_t limit) {}
  virtual void output_paired(uint32_t read_id, dna_slice read_seq, dna_slice pair_seq) = 0;
  virtual void output_unpaired(uint32_t read_id, dna_slice read_seq) = 0;
  virtual void done_chunk() {}
  virtual ~seqset_export_worker();

protected:
  seqset_export_worker() = default;
};

class seqset_export {
 public:
  seqset_export(const seqset_file* the_seqset,
                const readmap* the_readmap, std::string tmp_dir)
      : m_seqset(the_seqset), m_readmap(the_readmap), m_tmp_dir(tmp_dir) {}

  void prepare(progress_handler_t progress = null_progress_handler);

  // Outputs all reads. Paired reads will return two slices, and
  // unpaired slices will only return one.  The provided output function must be
  // threadsafe.
  //
  // It is guaranteed that all paired entries will be output before all unpaired
  // entries.
  void write_paired(const std::function<std::unique_ptr<seqset_export_worker>()>& worker_f,
                    progress_handler_t progress = null_progress_handler);
  void write_unpaired(const std::function<std::unique_ptr<seqset_export_worker>()>& worker_f,
                      progress_handler_t progress = null_progress_handler);

 private:
  const seqset_file* m_seqset = nullptr;
  const readmap* m_readmap = nullptr;
  std::string m_tmp_dir;

  boost::optional<seqset_flat> m_flat;
};
