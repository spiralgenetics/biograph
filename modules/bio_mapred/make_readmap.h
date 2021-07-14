#pragma once

#include "modules/io/version.h"
#include "modules/io/file_io.h"
#include "modules/io/packed_vector.h"
#include "modules/mapred/manifest.h"
#include "modules/bio_base/seqset_mergemap.h"
#include "modules/bio_base/seqset_bitmap.h"
#include "modules/bio_base/corrected_read.h"
#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"

class make_readmap
{
public:
	// Create a disk file with all the readmap data that can be loaded by
	// the constructor.
	static void do_make(
		const std::string& readmap_file_path
		, const seqset_file& the_seqset_file
		, manifest corrected_reads
		, bool is_paired
        , unsigned max_read_len
        , progress_handler_t progress = null_progress_handler
	);

	// Migrate the current readmap to point to a new seqset that presumably is
	// a superset of the current seqset.  If it's not a superset, the set the throw
	// on mismatch to false.
	static void migrate(
		const seqset_file& old_seqset
		, const readmap& old_readmap
		, const seqset_file& new_seqset
		, const std::string& new_readmap_path
		, bool throw_on_read_not_in_new_seqset = true
	);

  static void fast_migrate(
      const readmap& old_readmap,
      const seqset_mergemap& mergemap,
      const spiral_file_create_state& new_readmap,
      progress_handler_t progress = null_progress_handler);

  static void upgrade(
      const readmap& old_readmap,
      const seqset_file& the_seqset_file, const std::string& new_readmap_path,
      const std::function<dna_sequence(uint64_t /* seqset_id */, unsigned)>&
          lookup_seq,
      progress_handler_t progress = null_progress_handler);

 private:
  // This does mate loops in the same way as "Readmap" does, except
  // that MATE_RC does not point back to LOOP_START.
  enum loop_entry_type { LOOP_START, RC, MATE, MATE_RC };

  static constexpr unsigned k_entry_id_bits = 37;
  static constexpr unsigned k_read_length_bits = 10;
  static constexpr unsigned k_type_bits = 2;
  static constexpr uint64_t k_no_loop_entry = (1UL << k_entry_id_bits) - 1;
  struct mate_loop_table_entry {
    static_assert((k_entry_id_bits * 2 + k_read_length_bits * 2 + k_type_bits) == 8 * 12,
                  "mate loop table entry should take 12 bytes");
    static constexpr unsigned k_max_read_len = (1ULL << k_read_length_bits) - 1;

    uint64_t entry_id : k_entry_id_bits;
    loop_entry_type type : k_type_bits;
    uint32_t read_length : k_read_length_bits;
    uint32_t mate_read_length : k_read_length_bits;
    uint64_t loop_entry_id : k_entry_id_bits;

    constexpr mate_loop_table_entry()
        : entry_id(k_no_loop_entry),
          type(LOOP_START),
          read_length(0),
          mate_read_length(0),
          loop_entry_id(k_no_loop_entry) {}
    mate_loop_table_entry(loop_entry_type the_entry_type, uint64_t the_entry_id,
                          unsigned the_read_length, uint64_t the_loop_entry_id,
                          unsigned the_mate_read_length)
        : entry_id(the_entry_id),
          type(the_entry_type),
          read_length(the_read_length),
          mate_read_length(the_mate_read_length),
          loop_entry_id(the_loop_entry_id) {
      CHECK_EQ(sizeof(mate_loop_table_entry), 12);
      CHECK_LT(the_entry_id, k_no_loop_entry)
          << "Entry id too long to fit in mate loop table entry";
      CHECK_LE(the_read_length, k_max_read_len) << "Read too long to fit in mate loop table entry";
    }
  } __attribute__((packed));
  static_assert(sizeof(mate_loop_table_entry) == 12,
                "Mate loop table entry has an unexpected size");

  friend bool operator<(const mate_loop_table_entry& lhs, const mate_loop_table_entry& rhs);

  static std::string make_error_string(const seqset_file& old_seqset_file,
                                       const seqset_file& new_seqset_file,
                                       const dna_sequence& entry_seq);

  struct parallel_mate_loop_table_builder {
    parallel_mate_loop_table_builder(
        tracked_vector<make_readmap::mate_loop_table_entry>& the_mate_loop_table,
        const seqset* the_seqset, size_t file_info_count, size_t manifest_record_count,
        bool is_paired)
        : m_seqset(the_seqset), m_mate_loop_table(the_mate_loop_table),
          m_paired_counts(track_alloc("mate_loop_table:paired_counts")),
          m_unpaired_counts(track_alloc("mate_loop_table:unpaired_counts")), m_is_paired(is_paired) {
      SPLOG(
          "Constructing parallel_mate_loop_table_builder, file_info_count = %lu, "
          "manifest_record_count = %lu",
          file_info_count, manifest_record_count);
      m_paired_counts.resize(file_info_count);
      m_unpaired_counts.resize(file_info_count);
    }

    void operator()(const std::string&, const corrected_reads&, size_t file_info_id,
                    size_t record_id);

    std::pair<uint64_t, uint64_t> get_counts() const;  // Paired first, unpaired second.

    void add_read_to_mate_loop_table(
        loop_entry_type type, size_t record_id, uint64_t entry_id, unsigned read_len,
        uint64_t loop_entry_id, unsigned mate_read_len);

private:
		const seqset* m_seqset = nullptr;
		tracked_vector<make_readmap::mate_loop_table_entry>& m_mate_loop_table;
		tracked_vector<uint64_t> m_paired_counts;
		tracked_vector<uint64_t> m_unpaired_counts;
		bool m_is_paired;
		
		size_t calc_mate_loop_table_id(size_t record_id, bool is_paired) const
		{
			if (is_paired) {
				return 4 * record_id;
			} else {
				return 2 * record_id;
			}
		}
	};

    make_readmap(const seqset_file* the_seqset_file)
        : m_seqset_file(the_seqset_file),
          m_seqset(the_seqset_file?&the_seqset_file->get_seqset():nullptr),
          m_mate_loop_table(track_alloc("mate_loop_table")) {}

    void create_from_reads(manifest corrected_reads, bool is_paired,
                           unsigned max_read_len,
                           const spiral_file_create_state& state,
                           progress_handler_t progress = null_progress_handler);
    void create_from_migrate(const seqset_file& old_seqset_file,
                      const readmap& old_readmap,
                      const spiral_file_create_state& state,
                      bool throw_on_read_not_in_new_seqset);
    void create_from_fast_migrate(
        const readmap& old_readmap, const seqset_mergemap& mergemap,
        const spiral_file_create_state& new_readmap,
        progress_handler_t progress = null_progress_handler);

    void create_from_upgrade(
        const readmap& old_readmap, const seqset_file& the_seqset_file,
        const spiral_file_create_state& new_readmap,
        const std::function<dna_sequence(uint64_t /* seqset_id */, unsigned)>&
            lookup_seq,
        progress_handler_t progress = null_progress_handler);

   private:
	void create_common(const spiral_file_create_state& state,
                       const std::string& seqset_uuid, size_t seqset_size,
                       size_t num_reads, unsigned max_read_len);

	void import_reads_from(manifest corrected_reads_manifest, bool is_paired, progress_handler_t progress = null_progress_handler);

    const seqset_file* m_seqset_file = nullptr;
    const seqset* m_seqset = nullptr;
    std::unique_ptr<sparse_multi_builder> m_sparse_multi;
	std::unique_ptr<mutable_packed_varbit_vector> m_read_lengths;
    bool m_pairing_data_present = false;
    std::unique_ptr<mutable_packed_vector<uint32_t, 32>> m_mate_pair_ptr;
    std::unique_ptr<mutable_packed_varbit_vector> m_mate_loop_ptr;
    std::unique_ptr<mutable_packed_vector<unsigned, 1>> m_is_forward;

    // For creating from reads:
	tracked_vector<mate_loop_table_entry> m_mate_loop_table;
};

inline bool operator<(const make_readmap::mate_loop_table_entry& lhs,
                      const make_readmap::mate_loop_table_entry& rhs) {
  if (lhs.entry_id != rhs.entry_id) {
    return lhs.entry_id < rhs.entry_id;
  }
  if (lhs.type != rhs.type) {
    return lhs.type < rhs.type;
  }
  if (lhs.read_length != rhs.read_length) {
    return lhs.read_length < rhs.read_length;
  }
  if (lhs.mate_read_length != rhs.mate_read_length) {
    return lhs.mate_read_length < rhs.mate_read_length;
  }
  if (lhs.loop_entry_id != rhs.loop_entry_id) {
    return lhs.loop_entry_id < rhs.loop_entry_id;
  }

  // Change when https://gcc.gnu.org/bugzilla/show_bug.cgi?id=69285 is
  // fixed.
  constexpr bool k_buggy_parallel_partition = true;
  if (!k_buggy_parallel_partition) {
    DCHECK_EQ(0, memcmp(&lhs, &rhs, sizeof(make_readmap::mate_loop_table_entry)))
        << "Non-total ordering detected in readmap creation; result may not be deterministic";
  }
  return false;
}
