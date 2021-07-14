#include <gtest/gtest.h>
#include <sys/time.h>
#include <boost/filesystem.hpp>
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/log.h"
#include "modules/web/restful.h"

namespace fs = boost::filesystem;

class unittest_config_environment : public ::testing::Environment {
  void SetUp() override {
    int http_port = -1;
    // Try 100 random ports to find one we can use.

    // TODO(nils): We should pass a port number of 0 instead, and
    // let the OS pick an open port for us.  Also, we should
    // figure out why two processes can bind to the same port at
    // the same time during a bazel test, and then connect to each
    // other when they try to connect to themselves.
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    srand((11304120250909662091ULL * tv.tv_sec) ^ (18020468069336417183ULL * tv.tv_usec) ^
          (14238857486369442079ULL * getpid()));
    for (int i = 0; i < 100; i++) {
      int try_port = (rand() % (65536 - 1024)) + 1024;
      bind_info bi;
      bi.ip = "";
      bi.ssl = false;
      bi.port = try_port;
      bind_list_t bi_list;
      bi_list.push_back(bi);
      try {
        run_restful_server(bi_list, "", "", "thread", false /* don't block */);
        // No exception?  Success!
        http_port = try_port;
        break;
      } catch (...) {
        std::cerr << "Unable to bind to local port " << try_port << "\n";
        // Try a different port
      }
    }
    if (http_port == -1) {
      std::cerr << "Unable to find a local port to bind to.\n";
      exit(1);
    }

    std::string textual_port = printstring("%d", http_port);

    const int overwrite = 1;
    setenv("MASTER_PORT_5984_TCP_ADDR", "127.0.0.1", overwrite);
    setenv("MASTER_PORT_5984_TCP_PORT", textual_port.c_str(), overwrite);
    setenv("MASTER_PORT_5985_TCP_ADDR", "127.0.0.1", overwrite);
    setenv("MASTER_PORT_5985_TCP_PORT", textual_port.c_str(), overwrite);
    setenv("SPLOG_STDERR", "1", overwrite);

    std::string tmpdir;
    if (getenv("TEST_TMPDIR")) {
      tmpdir = std::string(getenv("TEST_TMPDIR")) + "/unittest_XXXXXX";
    } else {
      tmpdir = "/tmp/unittest_XXXXXX";
    }
    char buf[tmpdir.length() + 1];
    strcpy(buf, tmpdir.c_str());
    char* r = mkdtemp(buf);
    if (r == NULL) {
      throw io_exception("Unable to make temp directory " + tmpdir);
    }
    tmpdir = fs::canonical(buf).native();

    fs::create_directories(tmpdir);

    file_writer w(tmpdir + "/unittest.json");
    static const char json_config_template[] =
        "{\n"
        "	\"install_root\": \"/src\",\n"
        "	\"storage_root\" : \"file://%s/storage\",\n"
        "	\"path_allow_children\" : [\"s3://spiraleast/test\"],\n"
        "\n"
        "	\"S3_hostname\" : \"s3.amazonaws.com\",\n"
        "	\"access_key\" : \"\",\n"
        "	\"secret_key\" : \"\",\n"
        "	\"from_email\" : \"SpiralGeneticsJobBot@SpiralGenetics.com\",\n"
        "	\"email_id\" : \"\",\n"
        "	\"email_pass\" : \"\",\n"
        "\n"
        "	\"ottoman_bind_list\" : [ { \"port\" : %d } ],\n"
        "	\"taskdb_bind_list\" : [ { \"port\" : %d } ],\n"
        "\n"
        "	\"task_timeout\" : 600,\n"
        "	\"task_max_timeouts\" : 3,\n"
        "	\"taskdb_backup_period_in_seconds\" : 3,\n"
        "\n"
        "	\"temp_root\" : \"%s\",\n"
        "	\"resources_root\" : \"%s/storage/resources\",\n"
        "\n"
        "	\"test_root\" : \"%s\",\n"
        "	\"reference_path\" : \"%s/build_ref\"\n"
        "\n"
        "}\n";
    std::string json_config =
        printstring(json_config_template, tmpdir.c_str(), http_port, http_port, tmpdir.c_str(),
                    tmpdir.c_str(), tmpdir.c_str(), tmpdir.c_str());
    w.write(json_config.data(), json_config.size());
    w.close();
    boost::filesystem::create_directory(tmpdir + "/storage");
    boost::filesystem::create_directory(tmpdir + "/storage/resources");
    boost::filesystem::create_directory(tmpdir + "/storage/bulkdata");
    boost::filesystem::create_directory(tmpdir + "/storage/reference");
    boost::filesystem::create_directory(tmpdir + "/storage/reference/meta");
    boost::filesystem::create_directory(tmpdir + "/storage/build_ref");

    Config::load(tmpdir + "/unittest.json");
    // temp_root gets overridden by the TMPDIR environment variable somewhere I simply cannot
    // find. So override it manually here.
    Config::set("temp_root", tmpdir);

    log_init(nullptr, STDERR_FILENO);
    SPLOG("unittest_config_environment> Started up local http thread on port %d", http_port);
  }
};

static ::testing::Environment* const unittest_config_env __attribute__((unused)) =
    ::testing::AddGlobalTestEnvironment(new unittest_config_environment);
