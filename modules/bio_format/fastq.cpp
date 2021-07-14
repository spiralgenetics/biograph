#include "modules/bio_format/fastq.h"
#include "modules/io/simple_metadata.h"
#include "modules/io/utils.h"
#include "modules/io/log.h"
#include "modules/io/msgpack_transfer.h"
#include "modules/io/registry.h"

static const int k_maxline = 65536;

REGISTER_3(importer, fastq, readable&, bool, std::string const&);
REGISTER_3(exporter, fastq, writable&, bool, std::string const&);

fastq_reader::fastq_reader(readable& source, bool keep_quality)
    : m_source(source), m_keep_quality(keep_quality) {}

bool fastq_reader::read(std::string& key, std::string& value)
{
	read_id id;
	unaligned_reads reads;
	bool result = read(id, reads);
	if (!result) {
		return false;
	}
	key = msgpack_serialize(id);
	value = msgpack_serialize(reads);
	return true;
}

bool fastq_reader::read(read_id& id, unaligned_reads& reads) {
  reads.emplace_back();
  if (read(id, reads.back())) {
    return true;
  } else {
    reads.pop_back();
    return false;
  }
}

bool fastq_reader::read(read_id& id, unaligned_read& out) {
  const char* line;
  size_t len;

  // Read until non blank line or EOF
  while (true) {
    m_linenum++;
    if (!m_source.readline_no_copy(line, len, k_maxline)) {
      if (len > 0) {
        throw io_exception(
            printstring("line %d: Partial line in fastq file", m_linenum));
      }
      return false;
    }
    if (len > 0) {
      break;
    }
  }

  // Validate seq id, etc
  if (len < 2) {
    throw io_exception(
        printstring("line %d: Sequence id too short", m_linenum));
  }
  if (line[0] != '@') {
    throw io_exception(
        printstring("line %d: Sequence id missing @", m_linenum));
  }

  parse_read_name(std::string(line + 1, len - 1), id.pair_name, out);

  if (!m_source.readline_no_copy(line, len, k_maxline)) {
    throw io_exception(printstring(
        "line %d: End of file while reading sequence line", m_linenum));
  }

  m_linenum++;
  if (len == 0) {
    throw io_exception(printstring(
        "line %d: Expecting sequence, found empty line", m_linenum));
  }

  out.sequence = std::string(line, len);
  if (out.sequence.find_first_not_of("ACGTN") != std::string::npos) {
    throw io_exception(printstring(
        "line %d: Sequence contains unexpected characters", m_linenum));
  }

  if (!m_source.readline_no_copy(line, len, k_maxline)) {
    throw io_exception(
        printstring("line %d: End of file while reading + line", m_linenum));
  }

  m_linenum++;
  if (len == 0) {
    throw io_exception(
        printstring("line %d: Expecting +, found empty line", m_linenum));
  }
  if (line[0] != '+') {
    throw io_exception(
        printstring("line %d: Expecting + as first char of line", m_linenum));
  }

  if (!m_source.readline_no_copy(line, len, k_maxline)) {
    throw io_exception(printstring(
        "line %d: End of file while reading quality line", m_linenum));
  }

  if (m_keep_quality) {
    out.quality = std::string(line, len);

    for (char c : out.quality) {
      if (c < 33 or c > 126) {
        throw io_exception(printstring(
            "line %d: Quality line contains unexpected characters", m_linenum));
      }
    }
  }
  m_linenum++;

  if (len != out.sequence.size()) {
    throw io_exception(printstring(
        "line %d: Quality line not same length as sequence", m_linenum));
  }

  m_bases += out.sequence.size();

  return true;
}

void fastq_importer::import(kv_sink& sink, simple_metadata& meta)
{
	SPLOG("fastq_importer::import>");
	std::string key;
	std::string value;
	while (m_reader.read(key, value)) {
		sink.write(key, value);
	}
	SPLOG("fastq_importer::import> done");

	meta.set_simple("sample_bases", m_reader.get_bases());
}

void fastq_exporter::write(const std::string& key, const std::string& value)
{
	read_id id;
	unaligned_reads reads;

	msgpack_deserialize(id, key);
	msgpack_deserialize(reads, value);

	write(id, reads);
}

void fastq_exporter::write(const read_id& id, const unaligned_reads& reads)
{
	for (size_t i = 0; i < reads.size(); i++) {
		std::string name = build_read_name(id.pair_name, reads[i]);
		// Line 1: @<read name>
		m_sink.write("@", 1);
		m_sink.write(name);
		m_sink.write("\n", 1);
		// Line 2: Sequence
		m_sink.write(reads[i].sequence);
		m_sink.write("\n", 1);
		// Line 3: +<optional read name>
		m_sink.write("+\n", 2);
		// Line 4: Read quality
		m_sink.write(reads[i].quality);
		m_sink.write("\n", 1);
	}
}
