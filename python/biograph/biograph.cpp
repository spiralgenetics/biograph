#include <pybind11/pybind11.h>

#include "modules/bio_base/biograph.h"
#include "python/biograph/biograph.h"

using namespace pybind11;

void bind_biograph(module& m) {
  enum_<biograph::cache_strategy>(m, "CacheStrategy",
                                  R"DOC(
Access mode to use when opening a BioGraph.

The BioGraph Seqset is a large binary file with a nearly random access
pattern. This works well on local SSD storage and RAIDs, but it can
lead to performance issues on single magnetic hard disks or network
volumes (such as NFS or GPFS) that are optimized for sustained
sequential reads.

If performance issues are suspected, try opening the BioGraph with a
different caching strategy.

 * **MMAPCACHE**: memory-map and pre-cache (default)
 * **MMAP**: memory-map files without pre-caching
 * **RAM**: pull the entire BioGraph into RAM

Example:
  >>> # cache the BioGraph in RAM instead of using mmap
  >>> from biograph import BioGraph, CacheStrategy
  >>> my_bg = BioGraph('datasets/lambdaToyData/benchmark/family_lambda.bg', CacheStrategy.RAM)

Note:
  CacheStrategy.RAM uses significantly more memory than the other
  methods.
)DOC")
      .value("MMAP", biograph::cache_strategy::MMAP)
      .value("MMAPCACHE", biograph::cache_strategy::MMAPCACHE)
      .value("RAM", biograph::cache_strategy::RAM);

  class_<biograph, std::shared_ptr<biograph>>(m,  //
                                              "BioGraph",
                                              R"DOC(
Loads a BioGraph into an object in memory for querying.

Raises a *RuntimeError* if the BioGraph cannot be opened.

Args:
  path (str): The path to the BioGraph
  mode (CacheStrategy): Caching strategy to use for this BioGraph
                        (optional; default is MMAPCACHE)

Returns:
  BioGraph: A BioGraph object ready to query.

Example:
  >>> from biograph import BioGraph
  >>> my_bg = BioGraph('datasets/lambdaToyData/benchmark/family_lambda.bg')

Note:
  Performance may be negatively impacted if the BioGraph resides on
  network storage (such as an NFS or GPFS volume). See
  biograph.CacheStrategy to choose a different access mode.
)DOC")  //
      .def(init<const std::string&, biograph::cache_strategy>(), arg("dirname"),
           arg("strategy") = biograph::cache_strategy::MMAPCACHE)
      .def_property_readonly("seqset", &biograph::get_seqset, "Accessor to the Seqset object")
      .def("open_readmap", &biograph::open_readmap, arg("accession_id") = "",
           R"DOC(
open_readmap(accession_id)

Loads the Readmap of the specified accession ID into memory.
If None, tries to open the single Readmap in the BioGraph.

Args:
  accession_id (str): The accession ID of the Readmap to load

Returns:
  The loaded Readmap object.
)DOC")
      .def_property_readonly("metadata", &biograph::get_metadata,
                             "Accessor to the BioGraph file's metadata");
}
