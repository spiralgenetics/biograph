#include "gtest/gtest.h"
#include "modules/io/digest.h"
#include "modules/io/file_io.h"
#include "modules/io/log.h"
#include "modules/io/runtime_stats.h"
#include "modules/io/utils.h"

TEST(io, slurp) {
  std::string contents = "Test contents";
  {
    file_writer w("slurp_test_file");
    w.write(contents.data(), contents.size());
  }

  EXPECT_EQ(contents, slurp_file("slurp_test_file"));
}

TEST(io, sha1) {
  std::string ohai = "Hello world!";
  std::string sha = "d3486ae9136e7856bc42212385ea797094475802";

  EXPECT_EQ(sha1sum(ohai), sha);
}

TEST(io, sha256) {
  std::string ohai = "Hello world!";
  std::string sha = "c0535e4be2b79ffd93291305436bf889314e4a3faec05ecffcbb7df31ad9e51a";

  EXPECT_EQ(mdsum(ohai, "sha256"), sha);
}

TEST(io, sha1file) {
  std::string in_file = "golden/e_coli_10000snp.fq";
  std::string sha = "0c882474184fec1af00c815088a5c53cb88cff34";

  EXPECT_EQ(sha1sum(boost::filesystem::path(in_file)), sha);
}

TEST(io, runtime_stats) {
  // simple stats
  runtime_stats stats("stats.json");
  stats.add(js::Pair("foo", 42));
  stats.save();
  EXPECT_EQ(slurp_file("stats.json"), "{\"foo\":42,\"runtime_seconds\":0,\"stages\":[]}");

  // write somewhere else
  stats.save_to("another.json");
  // clear the stats
  stats.clear();
  // save
  stats.save();
  EXPECT_EQ(slurp_file("another.json"), "{\"runtime_seconds\":0,\"stages\":[]}");

  // add some stats
  stats.add(js::Pair("bar", 23.4));
  stats.add("another", "string");
  stats.save();
  EXPECT_EQ(
      slurp_file("another.json"),
      "{\"bar\":23.40000000000000,\"another\":\"string\",\"runtime_seconds\":0,\"stages\":[]}");

  // try an array
  js::Array ary;
  ary.push_back(1);
  ary.push_back(2.0);
  ary.push_back("three");
  stats.add(js::Pair("array", ary));
  // don't save yet
  EXPECT_EQ(
      slurp_file("another.json"),
      "{\"bar\":23.40000000000000,\"another\":\"string\",\"runtime_seconds\":0,\"stages\":[]}");

  // now write it out
  stats.save();
  EXPECT_EQ(slurp_file("another.json"),
            "{\"bar\":23.40000000000000,\"another\":\"string\",\"array\":[1,2.000000000000000,"
            "\"three\"],\"runtime_seconds\":0,\"stages\":[]}");

  // stages
  stats.clear();
  stats.save_to("stages.json");
  stats.add("something", "okay");
  stats.start_stage("begin");
  stats.end_stage("begin");
  stats.start_stage("middle");
  stats.end_stage("middle");
  stats.start_stage("end");
  stats.end_stage("end");
  stats.save();

  js::mValue stages;
  EXPECT_TRUE(js::read(slurp_file("stages.json"), stages));

  js::mObject s = stages.get_obj();

  EXPECT_EQ(s.at("something").get_str(), "okay");
  EXPECT_EQ(s.at("runtime_seconds").get_int(), 0);

  auto a = s.find("stages")->second.get_array();
  EXPECT_EQ(a.at(0).get_obj().at("name").get_str(), "begin");
  EXPECT_EQ((int)a.at(0).get_obj().at("wall_seconds").get_real(), 0);
  EXPECT_EQ((int)a.at(0).get_obj().at("cpu_user_seconds").get_real(), 0);
  EXPECT_EQ((int)a.at(0).get_obj().at("cpu_sys_seconds").get_real(), 0);

  EXPECT_EQ(a.at(1).get_obj().at("name").get_str(), "middle");
  EXPECT_EQ((int)a.at(1).get_obj().at("wall_seconds").get_real(), 0);
  EXPECT_EQ((int)a.at(1).get_obj().at("cpu_user_seconds").get_real(), 0);
  EXPECT_EQ((int)a.at(1).get_obj().at("cpu_sys_seconds").get_real(), 0);

  EXPECT_EQ(a.at(2).get_obj().at("name").get_str(), "end");
  EXPECT_EQ((int)a.at(2).get_obj().at("wall_seconds").get_real(), 0);
  EXPECT_EQ((int)a.at(2).get_obj().at("cpu_user_seconds").get_real(), 0);
  EXPECT_EQ((int)a.at(2).get_obj().at("cpu_sys_seconds").get_real(), 0);
}
