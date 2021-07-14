#include <pybind11/pybind11.h>
#include "python/common.h"

#include "python/biograph/biograph_metadata.h"

#include "modules/bio_base/biograph_dir.h"

using namespace pybind11;

namespace {

dict get_samples(const biograph_metadata& metadata) {
  dict result;
  for (const auto& elem : metadata.samples) {
    result[elem.first.c_str()] = elem.second;
  }
  return result;
}

}  // namespace

void bind_biograph_metadata(module& m) {
  class_<biograph_metadata>(m, "Metadata",
                            R"DOC(
Accessor for the BioGraph metadata, including accession ID, readmap ID, etc.

Example:
    >>> from biograph import BioGraph
    >>> my_sample = BioGraph('datasets/lambdaToyData/benchmark/proband_lambda.bg')
    >>> bg_metadata = my_sample.metadata
    >>> bg_metadata.accession_id
    'proband'
)DOC")
      .def_readonly("version", &biograph_metadata::version,
                    "The version of BioGraph used to create this file\n")
      .def_readonly("accession_id", &biograph_metadata::accession_id,
                    "The accession ID for this file\n")
      .def_readonly("biograph_id", &biograph_metadata::biograph_id,
                    "The unique BioGraph ID for this file\n")
      .def_property_readonly("samples", get_samples,
                             "Returns dictionary of key=accession_id, value=readmap_id for all "
                             "samples in this file\n");
}
