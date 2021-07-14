
#include "modules/main/main.h"
#include <boost/filesystem.hpp>

std::unique_ptr<Main> assemble_main(); // retired
std::unique_ptr<Main> biograph_info_main();
std::unique_ptr<Main> bwt_query_main();
std::unique_ptr<Main> discovery_main();
std::unique_ptr<Main> export_fastq_main();
std::unique_ptr<Main> export_main();
std::unique_ptr<Main> make_ref_main();
std::unique_ptr<Main> merge_seqset_main();
std::unique_ptr<Main> migrate_readmap_main();
std::unique_ptr<Main> ref2bwt_main();
std::unique_ptr<Main> ref2seqset_main();
std::unique_ptr<Main> seqset_dump_main();
std::unique_ptr<Main> seqset_main();
std::unique_ptr<Main> seqset_query_main();
std::unique_ptr<Main> upgrade_readmap_main();

void print_generic_help() {
  std::cerr << "bgbinary version " << biograph_current_version.make_string()
            << "\n\nCreate and manipulate BioGraph files.\n"
               "For more information, run any of the following commands with "
               "--help:\n"
               "\n"
               "  bgbinary create\n"
               "  bgbinary discovery\n"
               "  bgbinary reference\n"
               "  bgbinary metadata\n"
               "\n";
}

int main(int argc, char* argv[]) {
  spiral_init(&argc, &argv);
  if (argc == 1) {
    print_generic_help();
    return 1;
  }

  // Shift all args down by one, dropping the 'biograph'
  argc--;
  char* newargs[argc];
  for (int i = 1; i <= argc; i++) {
    newargs[i - 1] = argv[i];
  }

  srandom(time(0) * 0xffff + getpid());

  std::string program{boost::filesystem::basename(newargs[0])};

  std::map<std::string, main_f> programs = {
      {"create", seqset_main},
      {"metadata", biograph_info_main},
      {"merge", merge_seqset_main},
      {"reference", make_ref_main},
      {"upgrade", upgrade_readmap_main},
      {"variants", assemble_main}, // retired
      {"discovery", discovery_main},

      // Dev commands
      {"bwtquery", bwt_query_main},
      {"dump_flat", seqset_dump_main},
      {"dump_taskdb", dump_taskdb_main},
      {"export", export_main},
      {"migrate", migrate_readmap_main},
      {"query", seqset_query_main},
      {"ref2bwt", ref2bwt_main},
      {"ref2seqset", ref2seqset_main},
      {"rerun", rerun_main},
      {"resurrect", resurrect_main},
      {"export_fastq", export_fastq_main},
  };

  auto it = programs.find(program);
  if (it != programs.end()) {
    return it->second()->main(it->first, argc, newargs);
  }

  // List all known commands
  if (program == "COMMANDS") {
    for (auto it = programs.begin(); it != programs.end(); ++it) {
      std::cerr << it->first << std::endl;
    }
    std::cerr << std::endl;
    return 0;
  }

  print_generic_help();
  return 1;
}
