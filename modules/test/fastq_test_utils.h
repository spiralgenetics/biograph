#include <string> 

// imports fastq_file as MSGPACK-serialized key-values into kv_file
void make_fastq_kv(const std::string& fastq_file, const std::string& kv_file);
void make_zipped_fastq_kv(const std::string& zipped_fastq_file, const std::string& kv_file);

