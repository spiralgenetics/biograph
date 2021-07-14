
#pragma once

#include <openssl/aes.h>
#include <openssl/bio.h>
#include <openssl/conf.h>
#include <openssl/crypto.h>
#include <openssl/err.h>
#include <openssl/evp.h>
#include <openssl/pem.h>
#include <openssl/x509_vfy.h>
#include <openssl/x509v3.h>

#include "modules/io/io.h"
#include "modules/io/mem_io.h"

// Gives us 16 bytes of strong randomness
void generate_salt(char* salt);

// AES encryption for SpEC
class crypto_ctx {
 public:
  // Salt size must be 16 bytes long
  crypto_ctx(const char* salt, const char* key, size_t key_size);

  // Raw API
  // Tag must be 16 bytes long, out and in can be the same pointer
  void encrypt(uint64_t iv, char* tag, char* out, const char* in, size_t size) const;
  // Decrypt + validate tag, return false if invalid
  // If invalid, out may be overwritten with junk
  bool decrypt(uint64_t iv, const char* tag, char* out, const char* in, size_t size) const;

  // Block API
  // Encrypts + writes out block, destroys mem_io
  // Returns number of bytes written
  uint64_t encrypt(writable& out, uint64_t iv, mem_io&& block);
  // Reads block and decrypts and verifies, throws if crypto validation error
  void decrypt(mem_io& out, uint64_t iv, readable& in);

 private:
  AES_KEY m_key;
};

// RSA signature verification
class rsa_ctx {
 public:
  void Base64Encode(const unsigned char* buffer, size_t length, char** base64Text);
  void Base64Decode(const char* b64message, unsigned char** buffer, size_t* length);

  RSA* createPublicRSA(std::string key);

  bool RSAVerifySignature(RSA* rsa, unsigned char* MsgHash, size_t MsgHashLen, const char* Msg,
                          size_t MsgLen, bool* Authentic);

  bool verifySignature(std::string publicKey, std::string plainText, const char* signatureBase64);
};
