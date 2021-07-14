// message digest computation from OpenSSL

#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <openssl/evp.h>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

// string digest
std::string mdsum(const std::string& in_str, const std::string& method) {
  OpenSSL_add_all_digests();

  const EVP_MD* md = EVP_get_digestbyname(method.c_str());
  if (!md) {
    throw std::runtime_error("Unknown message digest: " + method);
  }

  EVP_MD_CTX* mdctx = EVP_MD_CTX_create();
  unsigned char md_value[EVP_MAX_MD_SIZE] = {0};
  uint32_t md_len;

  EVP_DigestInit_ex(mdctx, md, NULL);
  EVP_DigestUpdate(mdctx, in_str.c_str(), in_str.length());
  EVP_DigestFinal_ex(mdctx, md_value, &md_len);
  EVP_MD_CTX_destroy(mdctx);

  EVP_cleanup();

  std::ostringstream os;
  for (uint32_t i = 0; i < md_len; ++i) {
    os << std::hex << std::setw(2) << std::setfill('0') << (int)md_value[i];
  }
  return os.str();
}

// file digest
std::string mdsum(const fs::path& in_file, const std::string& method) {
  OpenSSL_add_all_digests();

  const EVP_MD* md = EVP_get_digestbyname(method.c_str());
  if (!md) {
    throw std::runtime_error("Unknown message digest: " + method);
  }

  EVP_MD_CTX* mdctx = EVP_MD_CTX_create();
  unsigned char md_value[EVP_MAX_MD_SIZE] = {0};
  uint32_t md_len;

  EVP_DigestInit_ex(mdctx, md, NULL);

  std::ifstream ifs(in_file.string(), std::ios::binary);

  if (not ifs.good()) {
    throw std::runtime_error("Cannot open " + in_file.string());
  }

  char buf[1024];

  while (ifs.good()) {
    ifs.read(buf, sizeof(buf));
    EVP_DigestUpdate(mdctx, buf, ifs.gcount());
  }
  ifs.close();

  EVP_DigestFinal_ex(mdctx, md_value, &md_len);
  EVP_MD_CTX_destroy(mdctx);

  EVP_cleanup();

  std::ostringstream os;
  for (uint32_t i = 0; i < md_len; ++i) {
    os << std::hex << std::setw(2) << std::setfill('0') << (int)md_value[i];
  }
  return os.str();
}

// Helper for running a sha1 digest
std::string sha1sum(const std::string& in_str) { return mdsum(in_str, "sha1"); };
std::string sha1sum(const boost::filesystem::path& in_file) { return mdsum(in_file, "sha1"); };
