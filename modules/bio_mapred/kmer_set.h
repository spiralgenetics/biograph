#pragma once

#include "modules/bio_base/dna_sequence.h"
#include "modules/mapred/manifest.h"
#include "modules/io/keyvalue.h"
#include "modules/io/mmap_buffer.h"
#include "modules/io/progress.h"
#include "modules/io/packed_vector.h"
#include "modules/bio_base/kmer.h"

#include <algorithm>
#include <boost/iterator/iterator_facade.hpp>

class kmer_set
{
public:
  // Flags associated with each kmer set entry.  0 means no flags are used.
  static constexpr unsigned k_flag_bits = 2;
  // Flags per entry
  static constexpr unsigned k_fwd_starts_read = 0b01;
  static constexpr unsigned k_rev_starts_read = 0b10;

  // Returned from find_table_index when the kmer was not found in the set.
  static constexpr size_t k_not_present = std::numeric_limits<size_t>::max();

	typedef std::function<void(
		size_t index,
		const kmer_t& kmer,
		size_t kmer_size,
		const std::string& value)
	> callback_t;

	inline static
	void null_callback(size_t index, const kmer_t& k, size_t ks, const std::string& v) {}

	typedef kmer_t key_type;
	typedef kmer_t value_type;
	typedef size_t size_type;
	typedef uint32_t lookup_t;

	class const_iterator : public boost::iterator_facade<
		const_iterator,
		kmer_t const,
		std::random_access_iterator_tag,
		kmer_t>
	{
	public:
		inline
		const_iterator(const kmer_set& self, size_t lookup_index, size_t table_index)
			: m_self(&self)
			, m_lookup_index(lookup_index)
			, m_table_index(table_index)
		{
			fixup();
		}

		inline bool equal(const_iterator const& other) const
		{
			return m_table_index == other.m_table_index;
		}

		kmer_t dereference() const;
    unsigned get_flags() const;

		void increment()
		{
			m_table_index++;
			fixup();
		}

    void decrement()
    {
      m_table_index--;
      seek_fixup();
    }

    void advance(difference_type distance) {
      m_table_index += distance;
      seek_fixup();
    }

		difference_type distance_to(const_iterator const& other) const
		{
			return difference_type(other.m_table_index) - difference_type(m_table_index);
		}

		inline size_t index() { return m_table_index; }

	private:
		void fixup()
		{
			while (m_self->m_lookup[m_lookup_index + size_t(1)] == m_table_index) {
				m_lookup_index++;
			}
		}

    void seek_fixup();

	private:
		const kmer_set* m_self;
		size_t m_lookup_index;
		size_t m_table_index;
	};

	// Build from a kv_source
	kmer_set(
		kv_source& source,
		size_t count,
		size_t kmer_size,
		const callback_t& progress = null_callback
	);

	// Build from a previous serialized kmer_set
	kmer_set(
		const std::string& serialized,
		const progress_handler_t& progress = null_progress_handler
	);

  // Build from a source of kmers.  "max_count" must be at least as
  // large as the number of kmers present, but should not be too much
  // larger as to not waste space.
  //
  // get_kmers must be able to make multiple passes through the list
  // of kmers.  kmers need not be sorted, and the kmer_output_f may be
  // called from multiple threads at once.
  using kmer_output_f = std::function<void(kmer_t, unsigned /* flags */)>;
  using kmer_source_f = std::function<void(const kmer_output_f&, progress_handler_t)>;
  kmer_set(size_t max_count, size_t kmer_size, size_t max_ram, const kmer_source_f& get_kmers,
           progress_handler_t progress = null_progress_handler);

	// Generate resource
	std::string save(
		const path& root,
		const progress_handler_t& progress = null_progress_handler
	);

	inline bool empty() const { return size() == 0; }
	inline size_t size() const { return m_size; }
	inline size_t kmer_size() const { return m_kmer_size; }
	const_iterator begin() const { return const_iterator(*this, 0, 0); }
	inline const_iterator end() const { return const_iterator(*this, m_lookup_size, m_size); }
	const_iterator find(kmer_t x) const;  // Get the index of a kmer, or npos
	size_t count(kmer_t x) const { return find_table_index(x) == k_not_present ? 0 : 1; }
  size_t find_table_index(kmer_t x) const;

	struct kmer_serialized;

  void copy_into_ram();
  unsigned get_flags(size_t index) const;

private:
  template <size_t /* tail bytes */>
  size_t sized_find_internal(kmer_t x) const;
	void create_sizes(size_t size, size_t kmer_size);
	void create_sizes_with_ram(size_t size, size_t kmer_size, size_t max_ram_bytes);
	void alloc_tables();
	void alloc_tables_in_memory();
	void save_memory_tables();
	kmer_t kmer_tail(size_t index) const;

  kmer_t lookup_for_kmer(kmer_t kmer) const { return kmer >> (2 * m_kmer_size - m_lookup_bits); }
  kmer_t tail_for_kmer(kmer_t kmer) const { return right(kmer, m_tail_bases); }
  kmer_t kmer_from_parts(size_t lookup, kmer_t tail) const;

private:
  class kmer_tail_reference;
  struct kmer_tail_and_flags_reference;
  struct kmer_tail_and_flags;
  class kmer_tail_iterator;
  using flags_table_t = mutable_packed_vector<unsigned, k_flag_bits>;
  using flags_reference = flags_table_t::mutable_accessor;
  size_t m_orig_size = 0;
	size_t m_size = 0;
	size_t m_kmer_size = 0;
	size_t m_lookup_bits = 0;
	size_t m_lookup_size = 0;
	size_t m_tail_bases = 0;
	size_t m_tail_bytes = 0;
	mmap_buffer m_lookup_buf;
	mmap_buffer m_table_buf;
	mutable_membuf m_lookup_membuf;
	mutable_membuf m_table_membuf;
	lookup_t* m_lookup = nullptr;
	unsigned char* m_table = nullptr;
  std::unique_ptr<flags_table_t> m_flags_table;
};

// Made this public for the GC, TODO: fix this
struct kmer_set::kmer_serialized
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(size);
		FIELD(kmer_size);
		FIELD(table);
		FIELD(lookup);
    FIELD(orig_size);
	}

	size_t size;
	size_t kmer_size;
	manifest table;
	manifest lookup;
  size_t orig_size;

	void validate();
};

