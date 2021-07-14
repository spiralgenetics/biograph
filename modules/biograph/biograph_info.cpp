#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include "modules/main/main.h"

#include "modules/bio_base/biograph_dir.h"
#include "modules/bio_base/seqset.h"
#include "modules/io/config.h"
#include "modules/io/file_io.h"
#include "modules/io/json_transfer.h"
#include "modules/io/log.h"
#include "modules/io/spiral_file_mmap.h"
#include "modules/io/version.h"

namespace fs = boost::filesystem;

class BiographInfoMain: public Main
{
public:
    BiographInfoMain()
    {
        m_usage =
            "%1% version %2%\n\n"
            "Usage: %1% <biograph> [biograph ...]\n\n"
            "Show metadata and other information about a BioGraph.\n"
        ;
    }
    size_t get_size(fs::path path);
  std::string get_size_fmt(fs::path path);
  void print_basic_info(biograph_dir bgdir);
protected:
    void add_args() override;
    int run(po::variables_map vars) override;
    const product_version& get_version() override { return biograph_current_version; }

private:
    std::vector<std::string> m_in_files;
    std::vector<std::string> m_sample_ids;
    std::string m_accession_id;

    po::options_description m_editing_options{"Metadata editing options", get_terminal_width()};
};

void BiographInfoMain::add_args()
{
    m_options.add_options()
        ("in", po::value(&m_in_files)->required(), "Input BioGraph")
    ;

    m_editing_options.add_options()
        ("accession-id", po::value(&m_accession_id), "Change the accession ID for the entire BioGraph")
        ("sample-id", po::value(&m_sample_ids)->multitoken(), "Change the accession ID for the specified sample, old=new")
    ;

    m_positional.add("in", -1);

    m_options.add(m_editing_options);
}

size_t BiographInfoMain::get_size(fs::path path) {
    size_t bytes = 4096; // Counting '.'
    fs::recursive_directory_iterator end;
    for (fs::recursive_directory_iterator rdi(path); rdi != end; rdi++) {
        if (fs::is_regular_file(rdi->path())) {
            bytes += fs::file_size(rdi->path());
        } else if (fs::is_directory(rdi->path())) {
            bytes += 4096;
        }
    }
    return bytes;
}

std::string BiographInfoMain::get_size_fmt(fs::path path) {
    size_t bytes = get_size(path);
    std::string suffix = "";

    if (bytes >= 1000 * 1000 * 1000) {
        suffix = " GB";
        bytes = bytes / (1000 * 1000 * 1000);
    } else if (bytes >= 1000 * 1000) {
        suffix = " MB";
        bytes = bytes / (1000 * 1000);
    } else {
        suffix = " KB";
        bytes = bytes / (1000);
    }
    return boost::str(boost::format("%d%s") % bytes % suffix);
}

void BiographInfoMain::print_basic_info(biograph_dir bgdir) {
  //Try to get some basic information out of the biograph
  spiral_file_open_mmap sf(bgdir.seqset());
  std::cout << "        command line: " << sf.file_info().command_line_str();

  std::cout << "\n\ncoverage file\n-------------\n";
  fs::directory_iterator end_itr;
  for(fs::directory_iterator i(fs::path(bgdir.path() / fs::path("coverage"))); i != end_itr; i++) {
      if(i->path().filename().string().find(".readmap") != std::string::npos) {
          std::cout << i->path().filename() << std::endl;
      }
  }
}

int BiographInfoMain::run(po::variables_map vars)
{
    // Sanity check: don't allow mass-changing of accession ids.
    if (m_in_files.size() > 1 and (!m_accession_id.empty() or !m_sample_ids.empty())) {
        throw std::runtime_error("You may only change the accession ID or sample ID on one BioGraph at a time.");
    }

    // Build a map of old -> new sample ids
    samples_t sample_map;
    if (not m_sample_ids.empty()) {
        for (auto s : m_sample_ids) {
            std::vector<std::string> kv;
            boost::split(kv, s, boost::is_any_of("="));
            if (kv.size() != 2) {
                throw std::runtime_error("Specify sample IDs to change as --sample-id old=new");
            }
            sample_map[kv[0]] = kv[1];
        }
    }

    for (const auto& in_file : m_in_files) {
        if (not fs::exists(in_file)) {
            std::cout << "*** Cannot open " << in_file << " (does not exist) ***\n\n\n";
            continue;
        }

        biograph_dir bgdir(in_file, READ_BGDIR);

        biograph_metadata meta = bgdir.get_metadata();
        if (meta.biograph_id.empty() or meta.samples.empty()) {
            std::cout << "*** Invalid metadata for " << in_file << " ***\n\n\n";
            print_basic_info(bgdir);
            continue;
        }

        seqset_file the_seqset_file(bgdir.seqset());
        const seqset& the_seqset = the_seqset_file.get_seqset();

        if (meta.biograph_id != the_seqset.uuid()) {
            std::cout << "*** Invalid metadata for " << in_file << " ***\n"
                << " (biograph_id is " << meta.biograph_id << ", not " << the_seqset.uuid() << ")\n\n\n";
            print_basic_info(bgdir);
            continue;
        }

        // sanity check: do all of the specified sample ids exist?
        auto samples = bgdir.samples();
        for (auto change : sample_map) {
            samples_t::iterator it = samples.find(change.first);
            if (it == samples.end()) {
                throw std::runtime_error("The sample ID '" + change.first + "' does not exist.");
            }
            it = samples.find(change.second);
            if (it != samples.end()) {
                throw std::runtime_error("The sample ID '" + change.second + "' already exists. Choose a unique name.");
            }
        }

        // Make sure we can write metadata if needed
        if ((!m_accession_id.empty()) or (!m_sample_ids.empty())) {
            bgdir.set_metadata(meta);
            bgdir.save_metadata();
        }

        std::cout <<   "       biograph path: " << fs::canonical(fs::path(in_file)) << "\n";

        if (m_accession_id.empty()) {
            std::cout << "        accession_id: " << meta.accession_id << "\n";
        } else {
            std::cout << "        accession_id: " << m_accession_id << " (was " << meta.accession_id << ")\n";
            meta.accession_id = m_accession_id;
        }

        spiral_file_open_mmap sf(bgdir.seqset());
        std::cout << "          created on: " << sf.file_info().create_timestamp_text << "\n";
        std::cout << "   number of samples: " << bgdir.samples().size() << "\n";
        std::cout << "   file size on disk: " << get_size_fmt(fs::path(in_file)) << "\n";
        std::cout << "  seqset entry count: " << the_seqset.size() << "\n";

        std::cout << "\n";
        std::cout << "    biograph version: " << meta.version << "\n";
        std::cout << "         biograph_id: " << meta.biograph_id << "\n";
        std::cout << "        command line: " << sf.file_info().command_line_str() << "\n";

        std::cout << "\n";
        std::cout << "coverage file" << std::string(37, ' ') << "sample id\n";
        std::cout << "-------------" << std::string(37, ' ') << "---------" << "\n";

        for (auto sample : bgdir.samples()) {
            if (fs::exists(bgdir.readmap(sample.second))) {
                std::cout << sample.second << ".readmap  ";
                if (sample_map.empty()) {
                    std::cout << sample.first << "\n";
                } else {
                    samples_t::iterator it = sample_map.find(sample.first);
                    if (it == sample_map.end()) {
                        std::cout << sample.first << "\n";
                    } else {
                        std::cout << it->second << " (was " << sample.first << ")\n";
                        meta.samples.erase(sample.first);
                        meta.samples[it->second] = sample.second;
                    }
                }
            } else {
                std::cout << "(missing)" << std::string(41, ' ') << sample.first << "\n";
            }
        }

        std::cout << "\ncommand history:\n";
        std::cout << "----------------\n";
        for (auto cmd = meta.command_history.rbegin(); cmd != meta.command_history.rend(); ++cmd) {
          std::cout << *cmd << "\n";
        }

        // Save updated metadata as needed
        if ((!m_accession_id.empty()) or (!m_sample_ids.empty())) {
            bgdir.set_metadata(meta);
            bgdir.save_metadata();
        }

        std::cout << "\n\n";
    }

    return 0;
}

std::unique_ptr<Main> biograph_info_main()
{
    return std::unique_ptr<Main>(new BiographInfoMain);
}


