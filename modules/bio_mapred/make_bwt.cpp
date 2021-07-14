#include <algorithm>
#include <parallel/algorithm>

#include "base/base.h"
#include "modules/bio_base/corrected_read.h"
#include "modules/bio_mapred/make_bwt.h"
#include "modules/io/bitcount.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/mmap_vector.h"
#include "modules/io/parallel.h"
#include "modules/io/progress.h"
#include "modules/io/uuid.h"
#include "modules/mapred/manifest_parallel.h"
#include "modules/mapred/resource_manager.h"
#include "modules/mapred/temp_file.h"
#include "modules/io/track_mem.h"

REGISTER_TASK(make_bwt_task);

void make_bwt_task::validate()
{
	SPLOG_P(LOG_DEBUG,"make_bwt_task::validate> cent_mod: %lu", cent_mod);
}

void make_bwt_task::run()
{
	validate();

	subprogress sp_compute([this](double p) { update_progress(p); }, 0.1, 1.0);

	SPLOG("make_bwt_task::run> Loading %s", input_ref.c_str());
	m_ref = make_unique<flat_ref>(input_ref);
	const flat_ref::index_t& idx = m_ref->get_index();

	size_t entries = 0;
	size_t cent_count = 0;
	for(const auto& ext : idx.extents) {
		entries += 1 + ext.size;
		cent_count += (ext.size + cent_mod - 1) / cent_mod;
	}

	SPLOG("make_bwt_task::run> total entries = %lu", entries);
	SPLOG("make_bwt_task::run> century entries = %lu", cent_count);

	tracked_vector<bwt_flyweight> flyweights(entries, track_alloc("make_bwt:bwt_flyweight"));
	size_t cur = 0;
	for(size_t i = 0; i < idx.extents.size(); i++) {
		for(size_t j = 0; j <= idx.extents[i].size; j++) {
			flyweights[cur].extent = i;
			flyweights[cur].offset = j;
			cur++;
		}
	}

	lambda_watchdog(sp_compute, 0.01, [&]() {
		dna_const_iterator flat_dna = m_ref->get_dna(0);
		parallel_sort_in_place(flyweights.begin(), flyweights.end(), [&](const bwt_flyweight& a, const bwt_flyweight& b) {
			const flat_ref::extent_t& aext = idx.extents[a.extent];
			const flat_ref::extent_t& bext = idx.extents[b.extent];
			dna_const_iterator ait = flat_dna + aext.flat + a.offset;
			dna_const_iterator bit = flat_dna + bext.flat + b.offset;
			size_t alen = aext.size - a.offset;
			size_t blen = bext.size - b.offset;
			dna_compare_result cmp = subseq_compare(ait, bit, alen, blen);
			if (cmp != dna_compare_result::EQUAL) {
              return cmp == dna_compare_result::FIRST_IS_LESS ||
                  cmp == dna_compare_result::FIRST_IS_PREFIX;
			}
			// Make sort stable by using original location as tie breaker
			if (a.extent != b.extent) {
				return a.extent < b.extent;
			}
			return a.offset < b.offset;
		});
		size_t bits_size = bitcount::compute_size(entries);
		// pointer to header + five bitmaps(A,C,T,G,Century) + century entries
		size_t file_size = 2*sizeof(uint64_t) + 5 * bits_size + sizeof(uint32_t) * cent_count;
		mmap_buffer bwt_file(output_bwt, file_size);

		char* base = bwt_file.buffer() + 2*sizeof(uint64_t);
		bitcount base_bits[4] = {
			{ base + (0 * bits_size), entries},
			{ base + (1 * bits_size), entries},
			{ base + (2 * bits_size), entries},
			{ base + (3 * bits_size), entries},
		};
		for(size_t i = 0; i < 4; i++) { base_bits[i].init(); }

		bitcount century_bits(base + (4 * bits_size), entries);
		century_bits.init();

		uint32_t* century_table = reinterpret_cast<uint32_t *>(base + (5 * bits_size));

		bwt_header header;
		size_t cur_century = 0;
		int cur_base = -1;
		for(size_t i = 0; i < entries; i++) {
			const flat_ref::extent_t& extent = idx.extents[flyweights[i].extent];
			uint64_t flat_pos = extent.flat + flyweights[i].offset;
			// Save C(a) table
			if (flyweights[i].offset != extent.size && int(*(flat_dna + flat_pos)) != cur_base && cur_base < 4) {
				cur_base = int(*(flat_dna + flat_pos));
				header.ca_table.push_back(i);
			}
			if (flyweights[i].offset != 0) {
				dna_base b = *(flat_dna + flat_pos - 1);
				base_bits[(int) b].set(i, true);
			}
			if (flyweights[i].offset != extent.size &&
				flyweights[i].offset % cent_mod == 0) {
				century_bits.set(i, true);
				century_table[cur_century++] = flat_pos;
			}
		}
		// C(a) table is now A, C, G, T, total entries, length of century table
		header.ca_table.push_back(entries);
		header.ca_table.push_back(cur_century);

		CHECK_EQ(cur_century, cent_count);
		CHECK(header.ca_table.size() == 6);

		for(size_t i = 0; i < 4; i++) { base_bits[i].finalize(); }
		century_bits.finalize();

		// Set first 8 bytes of file to k_magic
		*reinterpret_cast<uint64_t*>(bwt_file.buffer()) = bwt_file::k_magic;
		// Set next 8 bytes of file to point to position of msg_packed footer
		*reinterpret_cast<uint64_t*>(bwt_file.buffer() + sizeof(uint64_t)) = file_size;
		bwt_file.sync();
		bwt_file.close();

		file_writer fw(output_bwt, true);
		std::string str_header = msgpack_serialize(header);
		fw.write(str_header.data(), str_header.size());

		fw.close();
	});

	set_output(output_bwt);
}
