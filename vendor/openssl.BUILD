hdr_files = ["openssl/include/openssl/" + hdr for hdr in [
    "opensslconf.h",
    "__DECC_INCLUDE_EPILOGUE.H",
    "__DECC_INCLUDE_PROLOGUE.H",
    "aes.h",
    "asn1.h",
    "asn1_mac.h",
    "asn1err.h",
    "asn1t.h",
    "async.h",
    "asyncerr.h",
    "bio.h",
    "bioerr.h",
    "blowfish.h",
    "bn.h",
    "bnerr.h",
    "buffer.h",
    "buffererr.h",
    "camellia.h",
    "cast.h",
    "cmac.h",
    "cms.h",
    "cmserr.h",
    "comp.h",
    "comperr.h",
    "conf.h",
    "conf_api.h",
    "conferr.h",
    "crypto.h",
    "cryptoerr.h",
    "ct.h",
    "cterr.h",
    "des.h",
    "dh.h",
    "dherr.h",
    "dsa.h",
    "dsaerr.h",
    "dtls1.h",
    "e_os2.h",
    "ebcdic.h",
    "ec.h",
    "ecdh.h",
    "ecdsa.h",
    "ecerr.h",
    "engine.h",
    "engineerr.h",
    "err.h",
    "evp.h",
    "evperr.h",
    "hmac.h",
    "idea.h",
    "kdf.h",
    "kdferr.h",
    "lhash.h",
    "md2.h",
    "md4.h",
    "md5.h",
    "mdc2.h",
    "modes.h",
    "obj_mac.h",
    "objects.h",
    "objectserr.h",
    "ocsp.h",
    "ocsperr.h",
    "opensslv.h",
    "ossl_typ.h",
    "pem.h",
    "pem2.h",
    "pemerr.h",
    "pkcs12.h",
    "pkcs12err.h",
    "pkcs7.h",
    "pkcs7err.h",
    "rand.h",
    "rand_drbg.h",
    "randerr.h",
    "rc2.h",
    "rc4.h",
    "rc5.h",
    "ripemd.h",
    "rsa.h",
    "rsaerr.h",
    "safestack.h",
    "seed.h",
    "sha.h",
    "srp.h",
    "srtp.h",
    "ssl.h",
    "ssl2.h",
    "ssl3.h",
    "sslerr.h",
    "stack.h",
    "store.h",
    "storeerr.h",
    "symhacks.h",
    "tls1.h",
    "ts.h",
    "tserr.h",
    "txt_db.h",
    "ui.h",
    "uierr.h",
    "whrlpool.h",
    "x509.h",
    "x509_vfy.h",
    "x509err.h",
    "x509v3.h",
    "x509v3err.h",
]]

lib_files = [
    "openssl/libssl.a",
    "openssl/libcrypto.a",
]

cc_library(
    name = "openssl",
    srcs = lib_files,
    hdrs = hdr_files,
    includes = ["openssl/include"],
    linkopts = [
        "-ldl",
    ],
    linkstatic = 1,
    visibility = ["//visibility:public"],
    alwayslink = 1,
)

TCMALLOC_FLAGS = ("-fno-builtin-malloc " +
                  "-fno-builtin-free " +
                  "-fno-builtin-realloc " +
                  "-fno-builtin-calloc " +
                  "-fno-builtin-cfree " +
                  "-fno-builtin-memalign " +
                  "-fno-builtin-posix_memalign " +
                  "-fno-builtin-valloc " +
                  "-fno-builtin-pvalloc")

genrule(
    name = "openssl_build",
    srcs = glob(["openssl-1.1.1f/**/*"]),
    outs = lib_files + hdr_files,
    cmd = """
SRCDIR=$$(dirname $(location openssl-1.1.1f/Configure))
rm -rf $(@D)/openssl-build
cp -pr $${SRCDIR} $(@D)/openssl-build
(cd $(@D)/openssl-build/ &&
./Configure linux-x86_64 -fPIC no-rc2 no-rc4 no-rc5 no-camellia """ +
          #          select({
          #              "//conditions:default": "",
          #              "@spiral//tools:tcmalloc_setting": TCMALLOC_FLAGS,
          #          }) +
          #
          #  It's fine to specify TCMALLOC_FLAGS even if we're not using
          #  TCMALLOC, and it makes compilation take much longer if we have to
          #  recompile openssl every time we change the TCMALLOC setting.
          TCMALLOC_FLAGS +
          """ &&
make depend &&
make
) > openssl-build.log 2>&1 || (cat openssl-build.log; exit 1)
rm -rf $(@D)/openssl
cp -pr $(@D)/openssl-build $(@D)/openssl
    """,
)
