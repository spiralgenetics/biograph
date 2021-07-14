#include "base/base.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_format/dna_io.h"

void dna_writer::write(const dna_sequence& seq_to_write)
{
	m_packed_seq = seq_to_write.as_packed();
	CHECK(m_packed_seq.size() < 256);
	char packed_seq_size = static_cast<char>(m_packed_seq.size());
	m_target.write(&packed_seq_size, sizeof(packed_seq_size));
	m_target.write(m_packed_seq.data(), m_packed_seq.size());
}

dna_sequence dna_reader::read()
{
	dna_sequence seq_just_read;

	if (m_source) {
		char packed_seq_size;
		size_t amount_read = m_source->read(&packed_seq_size, sizeof(packed_seq_size));
		if (amount_read == 1) {
			m_packed_seq.resize(packed_seq_size);
			amount_read = m_source->read(&m_packed_seq[0], packed_seq_size);
			if (amount_read == static_cast<size_t>(packed_seq_size)) {
				seq_just_read = dna_sequence(m_packed_seq, true);
			} else {
				throw io_exception(printstring("dna_reader::read> Expected to read %lu bytes, but got %lu",
						m_packed_seq.size(), amount_read));
			}
		} // else EOF detected, return empty dna_sequence.
	}

	return seq_just_read;
}

void multi_file_dna_buffer::advance()
{
	m_current_sequence = m_dna_reader.read();
	while (at_eof() && m_current_file != m_source_files.end()) {
		m_dna_reader = dna_reader(make_unique<file_reader>((*m_current_file++)->path()));
		m_current_sequence = m_dna_reader.read();
	}
}
