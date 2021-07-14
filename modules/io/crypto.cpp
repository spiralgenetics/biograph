
#include "modules/io/crypto.h"
#include <endian.h>
#include <openssl/rand.h>
#include <openssl/sha.h>
#include <string.h>
#include "modules/io/hash.h"
#include "modules/io/hexdump.h"
#include "modules/io/log.h"

extern "C" {
// WTF openssl did you forget this one?
#include <openssl/modes.h>
}

constexpr size_t k_max_block_size = 64UL * 1024UL * 1024UL * 1024UL;

void generate_salt(char* salt) { RAND_bytes((unsigned char*)salt, 16); }

struct iv_struct {
  iv_struct(uint64_t _iv) : iv(htole64(_iv)), pad(0) {}
  uint64_t iv;
  uint32_t pad;
} __attribute__((packed));

struct block_header {
  block_header() {}
  block_header(uint32_t _block_size) : check(0), block_size(htole32(_block_size)), reserved(0) {}
  uint32_t get_size() {
    if (check != 0) {
      throw io_exception("Block decrypt failed: Check your decryption key");
    }
    return le32toh(block_size);
  }
  uint64_t check;  // Used to quick check key
  uint32_t block_size;
  uint32_t reserved;
} __attribute__((packed));

crypto_ctx::crypto_ctx(const char* salt, const char* key, size_t size) {
  unsigned char sha_key[32];
  // Hash the key + salt
  SHA256_CTX sha_ctx;
  SHA256_Init(&sha_ctx);
  SHA256_Update(&sha_ctx, salt, 16);
  SHA256_Update(&sha_ctx, key, size);
  SHA256_Final(sha_key, &sha_ctx);
  // Setup key, GCM only uses encrypt side
  AES_set_encrypt_key(sha_key, 256, &m_key);
}

void crypto_ctx::encrypt(uint64_t iv, char* tag, char* out, const char* in, size_t size) const {
  GCM128_CONTEXT* ctx;
  iv_struct ivs(iv);
  ctx = CRYPTO_gcm128_new(const_cast<void*>(static_cast<const void*>(&m_key)),
                          block128_f(AES_encrypt));
  CRYPTO_gcm128_setiv(ctx, (const unsigned char*)&ivs, sizeof(iv_struct));
  CRYPTO_gcm128_encrypt(ctx, (const unsigned char*)in, (unsigned char*)out, size);
  CRYPTO_gcm128_tag(ctx, (unsigned char*)tag, 16);
  // SPLOG("Tag: %s", hexdump(std::string(tag, 16)).c_str());
  CRYPTO_gcm128_release(ctx);
}

bool crypto_ctx::decrypt(uint64_t iv, const char* tag, char* out, const char* in,
                         size_t size) const {
  GCM128_CONTEXT* ctx;
  iv_struct ivs(iv);
  unsigned char vtag[16];
  ctx = CRYPTO_gcm128_new(const_cast<void*>(static_cast<const void*>(&m_key)),
                          block128_f(AES_encrypt));
  CRYPTO_gcm128_setiv(ctx, (const unsigned char*)&ivs, sizeof(iv_struct));
  CRYPTO_gcm128_decrypt(ctx, (const unsigned char*)in, (unsigned char*)out, size);
  CRYPTO_gcm128_tag(ctx, vtag, 16);
  // SPLOG("Tag: %s", hexdump(std::string((char*) vtag, 16)).c_str());
  CRYPTO_gcm128_release(ctx);
  return memcmp(tag, vtag, 16) == 0;
}

uint64_t crypto_ctx::encrypt(writable& out, uint64_t iv, mem_io&& block) {
  if (block.size() > k_max_block_size) {
    throw io_exception(boost::format("Encrypted block size %1% too large") % block.size());
  }

  GCM128_CONTEXT* ctx;
  // Setup context
  iv_struct ivs(iv);
  ctx = CRYPTO_gcm128_new(const_cast<void*>(static_cast<const void*>(&m_key)),
                          block128_f(AES_encrypt));
  CRYPTO_gcm128_setiv(ctx, (const unsigned char*)&ivs, sizeof(iv_struct));
  // Write out the block header (encrypted of course)
  block_header hdr(block.size());
  CRYPTO_gcm128_encrypt(ctx, (unsigned char*)&hdr, (unsigned char*)&hdr, sizeof(block_header));
  out.write((const char*)&hdr, sizeof(block_header));
  // Write out the main block
  CRYPTO_gcm128_encrypt(ctx, (unsigned char*)block.buffer(), (unsigned char*)block.buffer(),
                        block.size());
  out.write(block.buffer(), block.size());
  // Write the tag
  unsigned char tag_buf[16];
  CRYPTO_gcm128_tag(ctx, (unsigned char*)tag_buf, 16);
  out.write((const char*)tag_buf, 16);
  // Goodbye context
  CRYPTO_gcm128_release(ctx);
  return 16 + sizeof(block_header) + block.size();
}

void crypto_ctx::decrypt(mem_io& out, uint64_t iv, readable& in) {
  GCM128_CONTEXT* ctx;
  // Setup IV
  iv_struct ivs(iv);
  ctx = CRYPTO_gcm128_new(const_cast<void*>(static_cast<const void*>(&m_key)),
                          block128_f(AES_encrypt));
  CRYPTO_gcm128_setiv(ctx, (const unsigned char*)&ivs, sizeof(iv_struct));
  // Read the block header
  block_header hdr;
  size_t r = in.read((char*)&hdr, sizeof(block_header));
  if (r != sizeof(block_header)) {
    throw io_exception("EOF encountered while reading block header");
  }
  CRYPTO_gcm128_decrypt(ctx, (unsigned char*)&hdr, (unsigned char*)&hdr, sizeof(block_header));
  size_t size = hdr.get_size();
  ;
  // Check size, implode if too large
  if (size > k_max_block_size) {
    throw io_exception(boost::format("Decrypted block size %1% too large") % size);
  }
  // Read in the main block
  out.resize(size);
  out.reset();
  r = in.read(out.buffer(), size);
  if (r != size) {
    throw io_exception("EOF encountered while reading block data");
  }
  CRYPTO_gcm128_decrypt(ctx, (unsigned char*)out.buffer(), (unsigned char*)out.buffer(),
                        out.size());
  // Validate tag
  char tag_buf[16];
  unsigned char vtag_buf[16];
  r = in.read(tag_buf, 16);
  if (r != 16) {
    throw io_exception("EOF encountered while reading block tag");
  }
  CRYPTO_gcm128_tag(ctx, (unsigned char*)vtag_buf, 16);
  if (memcmp(tag_buf, vtag_buf, 16) != 0) {
    throw io_exception("Cyptographic checksum of block failed, data corruption");
  }
  // Goodbye context
  CRYPTO_gcm128_release(ctx);
}

// Basic RSA signature support inspired by:
// https://eclipsesource.com/blogs/2016/09/07/tutorial-code-signing-and-verification-with-openssl/
// https://gist.github.com/irbull/08339ddcd5686f509e9826964b17bb59

void rsa_ctx::Base64Encode(const unsigned char* buffer, size_t length, char** base64Text) {
  BIO *bio, *b64;
  BUF_MEM* bufferPtr;

  b64 = BIO_new(BIO_f_base64());
  bio = BIO_new(BIO_s_mem());
  bio = BIO_push(b64, bio);

  BIO_write(bio, buffer, length);
  BIO_flush(bio);
  BIO_get_mem_ptr(bio, &bufferPtr);
  BIO_set_close(bio, BIO_NOCLOSE);
  BIO_free_all(bio);

  *base64Text = (*bufferPtr).data;
}

size_t calcDecodeLength(const char* b64input) {
  size_t len = strlen(b64input), padding = 0;

  if (b64input[len - 1] == '=' && b64input[len - 2] == '=')  // last two chars are =
    padding = 2;
  else if (b64input[len - 1] == '=')  // last char is =
    padding = 1;
  return (len * 3) / 4 - padding;
}

void rsa_ctx::Base64Decode(const char* b64message, unsigned char** buffer, size_t* length) {
  BIO *bio, *b64;

  int decodeLen = calcDecodeLength(b64message);
  *buffer = (unsigned char*)malloc(decodeLen + 1);
  (*buffer)[decodeLen] = '\0';

  bio = BIO_new_mem_buf(b64message, -1);
  b64 = BIO_new(BIO_f_base64());
  bio = BIO_push(b64, bio);

  *length = BIO_read(bio, *buffer, strlen(b64message));
  BIO_free_all(bio);
}

// Make an RSA public key from a string
RSA* rsa_ctx::createPublicRSA(std::string key) {
  RSA* rsa = NULL;
  BIO* keybio;
  const char* c_string = key.c_str();
  keybio = BIO_new_mem_buf((void*)c_string, -1);
  if (keybio == NULL) {
    return 0;
  }
  rsa = PEM_read_bio_RSA_PUBKEY(keybio, &rsa, NULL, NULL);
  return rsa;
}

// Returns true if verification ran (ie. the signature could be checked).
// Sets Authentic to true if the signature is verified vs. the public key.
// You probably want crypto_ctx::verifySignature instead.
bool rsa_ctx::RSAVerifySignature(RSA* rsa, unsigned char* MsgHash, size_t MsgHashLen,
                                    const char* Msg, size_t MsgLen, bool* Authentic) {
  // default: fail closed
  *Authentic = false;

  EVP_PKEY* pubKey = EVP_PKEY_new();
  EVP_PKEY_assign_RSA(pubKey, rsa);
  EVP_MD_CTX* m_RSAVerifyCtx = EVP_MD_CTX_create();

  if (EVP_DigestVerifyInit(m_RSAVerifyCtx, NULL, EVP_sha256(), NULL, pubKey) <= 0) {
    return false;
  }
  if (EVP_DigestVerifyUpdate(m_RSAVerifyCtx, Msg, MsgLen) <= 0) {
    return false;
  }

  int AuthStatus = EVP_DigestVerifyFinal(m_RSAVerifyCtx, MsgHash, MsgHashLen);
  if (AuthStatus == 1) {
    *Authentic = true;
    return true;
  } else if (AuthStatus == 0) {
    return true;
  }

  // default: fail closed.
  return false;
}

// Verify that the publicKey used to sign plainText matches the Base64-encoded signature.
// Returns true on success.
bool rsa_ctx::verifySignature(std::string publicKey, std::string plainText,
                                 const char* signatureBase64) {
  RSA* publicRSA = createPublicRSA(publicKey);
  if (!publicRSA) {
    std::cerr << "rsa_ctx::verifySignature: invalid public key, aborting." << std::endl;
    return false;
  }
  unsigned char* encMessage;
  size_t encMessageLength;
  bool authentic;

  Base64Decode(signatureBase64, &encMessage, &encMessageLength);
  bool result = RSAVerifySignature(publicRSA, encMessage, encMessageLength, plainText.c_str(),
                                   plainText.length(), &authentic);
  return result & authentic;
}
