#include <cmath>
#include "modules/bio_base/kmer.h"
#include "modules/bio_mapred/kmerize_reads_mapper.h"
#include "modules/io/int_seq.h"
#include "modules/io/json_transfer.h"
#include "modules/io/log.h"

REGISTER_1(mapper, kmerize_reads, const std::string&);

namespace {

// We multiply all the logs by k_big_multiplier so that we can
// use integer arithmetic in the calculations.
const int64_t k_big_multiplier = 1000000000000000ULL; // 10^15
const double k_phred_recal = log(.1)/10.0;

} // anonymous namespace

constexpr int64_t calc_prob_inner(int qual)
{
	return qual > 33 ? static_cast<int64_t>(
		k_big_multiplier * log(1.0 - exp(
			k_phred_recal * static_cast<double>(qual - 33)
		))) : 0.0;
}

constexpr int64_t calc_prob(int qual)
{
	return qual == 33 ? calc_prob_inner(34) : calc_prob_inner(qual);
}

template<unsigned... Ints>
constexpr std::array<int64_t, 127> generate_array(int_seq<Ints...>)
{
	return {{ calc_prob(Ints)... }};
}

constexpr std::array<int64_t, 127> generate_array()
{
	return generate_array(generate_int_seq<127>{});
}

std::array<int64_t, 127> kmerize_reads_mapper::mg_log_lookup_table(generate_array());

void kmerize_reads_params::validate()
{
	SPLOG_P(LOG_DEBUG, "kmerize_reads_params::validate> kmer_size: %lu, trim: %lu, use_score: %s",
		kmer_size, trim, use_score ? "true" : "false"
	);

	if (kmer_size < 20 || kmer_size > 32) {
		SPLOG("Invalid kmer_size %lu (must be 20-32)", kmer_size);
		throw io_exception("Invalid kmer_size (must be 20-32)");
	}
	if ((trim > 0) && (trim <= kmer_size)) {
		SPLOG("trim (%lu) must be larger than kmer_size (%lu)", trim, kmer_size);
		throw io_exception("trim must be larger than kmer_size");
	}
	if (use_score != true && use_score != false) {
		SPLOG("use_score unspecified");
		throw io_exception("use_score unspecified");
	}
}

kmerize_reads_mapper::kmerize_reads_mapper(const std::string& params)
{
	json_deserialize(m_params, params);
	m_params.validate();
}

task_requirements kmerize_reads_mapper::get_requirements()
{
	return task_requirements {
		.profile = "normal",
		.cpu_minutes = 3,
	};
}

void kmerize_reads_mapper::typed_map(const read_id& key, const unaligned_reads& value)
{
	for (size_t i = 0; i < value.size(); i++) {
		map_one_read(key, value[i]);
	}
}

void kmerize_reads_mapper::map_one_read(const read_id& read_id, const unaligned_read& original_read)
{
	unaligned_read r = original_read;
	if (m_params.trim) {
		if (r.sequence.size() < m_params.trim) {
			return;
		}
		r.trim3(r.sequence.size() - m_params.trim);
	}

	if (r.sequence.size() < m_params.kmer_size) {
		return;
	}

	dna_sequence read_sequence(r.sequence);
	kmer_t kmer = make_kmer(read_sequence.begin(), m_params.kmer_size);
	int64_t log_prob = 0;

	uint32_t score = 1;
	if (m_params.use_score) {
		log_prob = check_qual(read_id, r.quality, 0);
		score = static_cast<unsigned>(log(1.0 - exp(
			static_cast<double>(log_prob) / k_big_multiplier)) / k_phred_recal);
	}
	kmer_t ckmer = canonicalize(kmer, m_params.kmer_size);
	if (ckmer == kmer) {
		output(ckmer, kcount_pair(score, 0));
	} else {
		output(ckmer, kcount_pair(0, score));
	}

	for (size_t offset = 1; offset <= r.sequence.size() - m_params.kmer_size; offset++) {
		int base = (int) read_sequence[m_params.kmer_size + offset - 1];
		rotate_left(kmer, m_params.kmer_size, base);

		score = 1;
		if (m_params.use_score) {
			log_prob = rotate_qual(read_id, r.quality, offset, log_prob);
			auto score = static_cast<unsigned>(log(1.0 - exp(
				static_cast<double>(log_prob) / k_big_multiplier)) / k_phred_recal);
			if (score == 0) {
				continue;
			}
		}
		kmer_t ckmer = canonicalize(kmer, m_params.kmer_size);
		if (ckmer == kmer) {
			output(ckmer, kcount_pair(score, 0));
		} else {
			output(ckmer, kcount_pair(0, score));
		}
	}
}

int64_t kmerize_reads_mapper::check_qual(
	const read_id& read_id,
	const std::string& qual,
	size_t start)
{
	int64_t log_prob = 0;
	for (size_t i = start; i < start + m_params.kmer_size; i++) {
		unsigned char base_qual = qual[i];
		if (base_qual < 33) {
			std::string error_string("Invalid quality score of ");
			error_string += std::to_string(base_qual);
			error_string += " in read ";
			error_string += read_id.pair_name;
			throw io_exception(error_string);
		}

		log_prob += kmerize_reads_mapper::mg_log_lookup_table[base_qual];
	}
	return log_prob;
}

int64_t kmerize_reads_mapper::rotate_qual(
	const read_id& read_id,
	const std::string& qual,
	size_t start,
	int64_t score)
{
	unsigned char base_qual = qual[start + m_params.kmer_size - 1];
	if (base_qual < 33) {
		std::string error_string("Invalid quality score of ");
		error_string += std::to_string(base_qual);
		error_string += " in read ";
		error_string += read_id.pair_name;
		throw io_exception(error_string);
	}

	return score - kmerize_reads_mapper::mg_log_lookup_table[qual[start - 1]]
		+ kmerize_reads_mapper::mg_log_lookup_table[base_qual];
}

